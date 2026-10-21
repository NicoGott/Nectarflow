#!/usr/bin/env python3
#
# Resolves the reference library once and writes assets/reference_manifest.tsv:
# one row per sequence, pinned to accession.version, with the coordinates it is
# trimmed to, the md5 of the trimmed sequence and a taxonomy frozen against a
# dated NCBI taxdump. Run by hand, not by the pipeline; build_reference.py then
# rebuilds the library from the manifest and searches nothing.
#
# Which plants belong in the library comes from one published list of plants,
# the Plant Atlas 2020 (bin/plant_atlas_taxa.py). The sequences come from one
# archive, GenBank. The authors' library also covered the garden's own living
# collection, which is not ours to redistribute, so it is not a source here. Barcode UK (Jones, Twyford et al. 2021) has precedence
# inside it: a taxon with a Barcode UK barcode gets that and nothing is fetched
# for it. Nothing from the authors' results enters at any point.
#
import argparse, collections, csv, gzip, hashlib, os, re, subprocess, sys, tempfile, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from entrez import Entrez, parse_fasta
from refnames import Taxdump, binomial, genus, normalise

BARCODE_UK = "JN890545:JN896265,KX165423:KX167996,MK924423:MK926404"
MARKER_TERMS = {"rbcL": "rbcL", "ITS2": '"internal transcribed spacer 2" OR ITS2', "matK": "matK"}
# Search window: a record shorter than the shortest query cannot cover it, and one longer than
# this is a plastid genome or a full ITS region we would rather not download. Floors are the
# pipeline's own min_len (450 rbcL, 300 ITS2) less the primer fragments trim_primers.py removes.
SLEN = {"rbcL": (450, 5000), "ITS2": (300, 5000), "matK": (450, 5000)}
FLOOR = {"rbcL": 420, "ITS2": 270, "matK": 420}
# The extent a Barcode UK barcode covers, which is the region every other record is trimmed to
REGION = {"rbcL": 542, "ITS2": 446, "matK": 800}
JUNK = re.compile(r"UNVERIFIED|uncultured|environmental sample|whole genome shotgun|"
                  r"complete genome|complete plastid|chloroplast, complete", re.I)
CMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
COLS = ("marker accession source target target_rank target_name ncbi_taxid ncbi_name "
        "species genus family order strand start end length md5 lineage").split()


#
# Reverse complement
def rc(s):
    return s.translate(CMP)[::-1]


#
# Reads a committed snapshot, skipping the comment header
def read_tsv(path):
    return list(csv.DictReader((l for l in open(path, encoding="utf-8") if not l.startswith("#")), delimiter="\t"))


#
# Turns the two snapshots into targets. A species-rank target is an accepted
# binomial; a genus-rank target is a genus of either list. Names that are not
# binomials, that no Viridiplantae taxid matches, or that match several, are
# written to the rejects file with the reason - that file is the denominator of
# every coverage figure, so nothing is dropped in silence
def read_targets(td, atlas):
    want, rejects, genera = {}, [], {}
    src = [(r["taxon"], "plant_atlas", "species") for r in read_tsv(atlas)]
    for name, source, rank in src:
        g = genus(name)
        if g:
            genera.setdefault(g, source)
        if rank == "genus":
            continue
        b = binomial(name)
        if not b:
            rejects.append((name, source, "not a binomial"))
            continue
        ids = td.lookup(b)
        if not ids:
            rejects.append((name, source, "no Viridiplantae taxid"))
            continue
        ids = [t for t in ids if td.rank.get(t) in ("species", "subspecies", "varietas", "forma")] or ids
        if len(set(ids)) > 1:
            rejects.append((name, source, "ambiguous: " + ",".join(str(i) for i in ids)))
            continue
        t = ids[0]
        key = ("species", t)
        if key not in want:
            want[key] = {"taxid": t, "rank": "species", "name": b, "source": source,
                         "ncbi_name": td.sci.get(t, b)}
    # every genus of either list is a genus-rank candidate; whether it is actually
    # needed is decided per marker, once we know which species targets got a sequence
    for g, source in sorted(genera.items()):
        ids = [t for t in td.lookup(g) if td.rank.get(t) == "genus"]
        if len(set(ids)) != 1:
            rejects.append((g, source, "genus: " + ("ambiguous" if ids else "no Viridiplantae taxid")))
            continue
        want[("genus", ids[0])] = {"taxid": ids[0], "rank": "genus", "name": g, "source": source,
                                   "ncbi_name": td.sci.get(ids[0], g)}
    return want, rejects


#
# Downloads Barcode UK whole, by the accession ranges of its data statement,
# and keeps the records whose taxid really is a plant. These go into the
# library verbatim: they are the barcodes the region is defined by
def fetch_barcode_uk(ez, markers, ranges, td, log):
    out = {}
    acc = " OR ".join(f"{r}[ACCN]" for r in ranges.split(","))
    for m in markers:
        term = f"({acc}) AND ({MARKER_TERMS[m]})"
        n, webenv, key = ez.esearch("nuccore", term)
        meta = {}
        for r in ez.esummary_all("nuccore", webenv, key, n):
            if r["taxid"] not in td.keep:
                log.write(f"barcode_uk\t{m}\t{r['accession']}\tdropped\tnot Viridiplantae (taxid {r['taxid']})\n")
                continue
            meta[r["accession"]] = r
        seqs = ez.fetch_fasta(sorted(meta))
        missing = sorted(set(meta) - set(seqs))
        if missing:
            raise SystemExit(f"barcode_uk {m}: {len(missing)} accessions not returned, first {missing[:3]}")
        for a, r in meta.items():
            r["seq"] = seqs[a]
        out[m] = meta
        log.write(f"barcode_uk\t{m}\tsearch\t{n}\tkept {len(meta)}\n")
        print(f"barcode_uk {m}: {n} hits, {len(meta)} plant records", file=sys.stderr)
    return out


#
# The templates that define the barcode region: one Barcode UK record per
# family, per marker, the one whose length is the modal one for that family,
# smallest accession breaking the tie. A candidate is searched against all of
# them at once and keeps its best hit, so there is no fallback chain to follow
def exemplars(bu, td):
    out = {}
    for m, recs in bu.items():
        pools = collections.defaultdict(list)
        for a, r in sorted(recs.items()):
            if abs(len(r["seq"]) - REGION[m]) > REGION[m] * 0.15:
                continue
            fam = td.at_rank(r["taxid"], "family") or "*"
            pools[fam].append((a, r["seq"]))
        pick = {}
        for v in pools.values():
            modal = collections.Counter(len(s) for _, s in v).most_common(1)[0][0]
            a, seq = min((x for x in v if len(x[1]) == modal), key=lambda x: x[0])
            pick[a] = seq
        out[m] = pick
    return out


#
# Where the barcode region sits inside each candidate, found with vsearch rather
# than by hand. Every candidate of a marker is searched once against the whole
# set of templates and keeps its best hit; the aligned extent of that hit is the
# region. This is what stops a 2.9 kb rbcL gene or a full ITS1-5.8S-ITS2 record
# from competing with a 542 bp barcode: every sequence in the library covers one
# region. Returns {accession: (strand, start, end)}, 0-based, end exclusive
def anchor_all(m, seqs, templates, vsearch, tmpdir, min_id):
    if not seqs:
        return {}
    q, db, out = (os.path.join(tmpdir, f"{m}.{x}") for x in ("query.fasta", "templates.fasta", "hits.tsv"))
    with open(q, "w") as fh:
        for acc, seq in sorted(seqs.items()):
            fh.write(f">{acc}\n{seq}\n")
    with open(db, "w") as fh:
        for acc, seq in sorted(templates.items()):
            fh.write(f">{acc}\n{seq}\n")
    subprocess.run([vsearch, "--usearch_global", q, "--db", db, "--id", str(min_id), "--strand", "both",
                    "--maxaccepts", "1", "--top_hits_only", "--userout", out,
                    "--userfields", "query+qilo+qihi+qstrand", "--quiet"], check=True)
    hits = {}
    with open(out) as fh:
        for line in fh:
            acc, lo, hi, strand = line.rstrip("\n").split("\t")
            hits[acc] = (strand, int(lo) - 1, int(hi))
    print(f"{m}: vsearch placed {len(hits)} of {len(seqs)} candidates in the barcode region",
          file=sys.stderr)
    return hits


#
# One manifest row. Barcode UK records go in whole; everything else is trimmed
# to the template and must still be at least as long as the shortest query the
# pipeline accepts, or it is refused and the next candidate is tried
def row(m, acc, seq, taxid, td, source, target, hit):
    if hit == "whole":                       # a Barcode UK record is already the barcode
        strand, a, b = "+", 0, len(seq)
    elif hit is None:
        return None, "vsearch found no barcode region in this sequence"
    else:
        strand, a, b = hit
        if b - a < FLOOR[m]:
            return None, f"trimmed to {b - a} bp, below the {FLOOR[m]} bp floor"
    s = (rc(seq) if strand == "-" else seq)[a:b]
    sp = td.at_rank(taxid, "species")
    return {"marker": m, "accession": acc, "source": source, "target": target["source"],
            "target_rank": target["rank"], "target_name": target["name"], "ncbi_taxid": taxid,
            "ncbi_name": td.sci.get(taxid, ""),
            "species": binomial(sp) if target["rank"] == "species" else "",
            "genus": td.at_rank(taxid, "genus"), "family": td.at_rank(taxid, "family"),
            "order": td.at_rank(taxid, "order"), "strand": strand, "start": a + 1, "end": b,
            "length": len(s), "md5": hashlib.md5(s.encode()).hexdigest(),
            "lineage": ";".join(td.lineage(taxid)), "_seq": s}, ""


#
# Orders the candidates for one target: a record that can contain the whole
# barcode region first, shortest of those first because it carries least
# besides the barcode; then the partial ones, longest first. Accession breaks
# every tie, so the order does not depend on what NCBI sorted by
def rank_key(c, m):
    long_enough = c["length"] >= REGION[m]
    return (0 if long_enough else 1, c["length"] if long_enough else -c["length"], c["accession"])


#
# Searches GenBank once per genus per marker over the whole hit set - no retmax
# cap, and the order NCBI returns is never used - and returns the candidates
# that pass the title screen and are plants, grouped by genus taxid
def search_genera(ez, marker, gids, td, log):
    out = {}
    for i, g in enumerate(sorted(gids), 1):
        lo, hi = SLEN[marker]
        term = f"txid{g}[Organism:exp] AND ({MARKER_TERMS[marker]}) AND {lo}:{hi}[SLEN]"
        n, webenv, key = ez.esearch("nuccore", term)
        keep = []
        if n:
            for r in ez.esummary_all("nuccore", webenv, key, n):
                if JUNK.search(r["title"]) or r["taxid"] not in td.keep or not r["accession"]:
                    continue
                keep.append(r)
        out[g] = keep
        log.write(f"search\t{marker}\ttxid{g}\t{td.sci.get(g, '')}\t{n}\t{len(keep)}\t{term}\n")
        if i % 50 == 0:
            print(f"  {marker}: {i}/{len(gids)} genera searched", file=sys.stderr)
    return out


#
# Walks a taxid up to the root once, so a candidate can be matched to the
# targets above it without re-walking for every target
def ancestors(td, t):
    out, x = [], t
    while x and x != 1:
        out.append(x)
        x = td.parent.get(x, 0)
    return out


#
# The whole resolve: targets, Barcode UK, the gaps, the GenBank fill, the trim
# and the manifest
def main():
    ap = argparse.ArgumentParser(description="Resolves the reference library and writes the pinned manifest")
    ap.add_argument("--plant-atlas", required=True)
    ap.add_argument("--taxdump", required=True, help="a DATED taxdmp zip from NCBI's taxdump_archive")
    ap.add_argument("--ranges", default=BARCODE_UK)
    ap.add_argument("--markers", default="rbcL,ITS2")
    ap.add_argument("--max-species", type=int, default=3, help="sequences per species-rank target per marker")
    ap.add_argument("--max-genus-species", type=int, default=10, help="species per genus-rank target per marker")
    ap.add_argument("--genus-sample", type=int, default=0, help="resolve a random sample of genera (for costing)")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--email", required=True)
    ap.add_argument("--cache", default=".entrez_cache")
    ap.add_argument("--anchor-id", type=float, default=0.7,
                    help="vsearch --id for placing a candidate against a template")
    ap.add_argument("--vsearch", default="vsearch",
                    help="vsearch binary; used to place each candidate in the barcode region")
    ap.add_argument("--out", required=True)
    ap.add_argument("--rejects", required=True)
    ap.add_argument("--log", required=True)
    ap.add_argument("--vendor", help="directory to write the trimmed <marker>.fasta.gz the build reads")
    a = ap.parse_args()
    t0 = time.time()
    markers = a.markers.split(",")
    td_md5 = hashlib.md5(open(a.taxdump, "rb").read()).hexdigest()
    print(f"taxdump {os.path.basename(a.taxdump)} md5 {td_md5}", file=sys.stderr)
    td = Taxdump(a.taxdump)
    print(f"taxdump loaded: {len(td.keep)} Viridiplantae nodes, {len(td.index)} names", file=sys.stderr)

    log = open(a.log, "w", encoding="utf-8")
    log.write("# resolver log; every search term, hit count and decision\n")
    want, rejects = read_targets(td, a.plant_atlas)
    sp_t = {k: v for k, v in want.items() if v["rank"] == "species"}
    gn_t = {k: v for k, v in want.items() if v["rank"] == "genus"}
    print(f"targets: {len(sp_t)} species-rank, {len(gn_t)} genus-rank, {len(rejects)} rejected", file=sys.stderr)

    ez = Entrez(email=a.email, cache=a.cache)
    bu = fetch_barcode_uk(ez, markers, a.ranges, td, log)
    tmpl = exemplars(bu, td)
    rows, used = [], set()

    for m in markers:
        for acc in sorted(bu[m]):
            r = bu[m][acc]
            tgt = {"source": "barcode_uk", "rank": "species", "name": r["organism"]}
            got, why = row(m, acc, r["seq"], r["taxid"], td, "barcode_uk", tgt, "whole")
            if got:
                rows.append(got)
                used.add((m, acc))

    # a target is satisfied for a marker when some Barcode UK record sits inside its subtree
    covered = {m: set() for m in markers}
    for m in markers:
        for r in bu[m].values():
            covered[m] |= set(ancestors(td, r["taxid"]))

    # which genera to search: those with an unsatisfied species target, plus every
    # genus of either list, because a genus with no representative at all cannot
    # even be called 'Genus spp.' and its reads go to the nearest relative instead
    need = {m: sorted({v["taxid"] for v in sp_t.values() if v["taxid"] not in covered[m]}) for m in markers}
    genus_of = {}
    for v in sp_t.values():
        gid = next((x for x in ancestors(td, v["taxid"]) if td.rank.get(x) == "genus"), 0)
        genus_of[v["taxid"]] = gid or v["taxid"]
    allg = sorted({genus_of[v["taxid"]] for v in sp_t.values()} | {v["taxid"] for v in gn_t.values()})
    if a.genus_sample:
        import random
        random.Random(a.seed).shuffle(allg)
        allg = sorted(allg[:a.genus_sample])
        print(f"COSTING RUN: a random {len(allg)} of the genera, seed {a.seed}", file=sys.stderr)

    for m in markers:
        pool = search_genera(ez, m, allg, td, log)
        by_target = collections.defaultdict(list)
        by_genus_species = collections.defaultdict(lambda: collections.defaultdict(list))
        for g, cands in pool.items():
            for c in cands:
                anc = set(ancestors(td, c["taxid"]))
                for t in need[m]:
                    if t in anc:
                        by_target[t].append(c)
                spid = next((x for x in ancestors(td, c["taxid"]) if td.rank.get(x) == "species"), 0)
                if spid:
                    by_genus_species[g][spid].append(c)
        want_seq = []
        for t in need[m]:
            tgt = sp_t[("species", t)]
            want_seq.append((tgt, sorted(by_target.get(t, []), key=lambda c: rank_key(c, m)), a.max_species))
        # genus-rank filling, for the genera that still have nothing in this marker
        filled = {g for g in allg if any(by_target.get(t) for t in need[m] if genus_of.get(t) == g)}
        filled |= {g for g in allg if g in covered[m]}
        for g in allg:
            if g in filled or not by_genus_species.get(g):
                continue
            tgt = gn_t.get(("genus", g)) or {"taxid": g, "rank": "genus", "name": td.sci.get(g, str(g)),
                                             "source": "plant_atlas", "ncbi_name": td.sci.get(g, "")}
            tgt = dict(tgt, rank="genus")
            for spid in sorted(by_genus_species[g])[:a.max_genus_species]:
                want_seq.append((tgt, sorted(by_genus_species[g][spid], key=lambda c: rank_key(c, m)), 1))
        # fetch a few more than needed, because a candidate that will not trim is refused
        need_acc = []
        for tgt, cands, cap in want_seq:
            for c in cands[:cap + 3]:
                if (m, c["accession"]) not in used:
                    need_acc.append(c["accession"])
        need_acc = sorted(set(need_acc))
        print(f"{m}: fetching {len(need_acc)} candidate sequences", file=sys.stderr)
        seqs = ez.fetch_fasta(need_acc)
        with tempfile.TemporaryDirectory() as tmpdir:
            placed = anchor_all(m, seqs, tmpl[m], a.vsearch, tmpdir, a.anchor_id)
        gone = sorted(set(need_acc) - set(seqs))
        for x in gone:
            log.write(f"fetch\t{m}\t{x}\tmissing\tnot returned by efetch\n")
        for tgt, cands, cap in want_seq:
            kept = 0
            for c in cands:
                if kept >= cap:
                    break
                acc = c["accession"]
                if (m, acc) in used or acc not in seqs:
                    continue
                got, why = row(m, acc, seqs[acc], c["taxid"], td, "genbank", tgt,
                               placed.get(acc))
                if got is None:
                    log.write(f"refused\t{m}\t{acc}\t{tgt['name']}\t{why}\n")
                    continue
                rows.append(got)
                used.add((m, acc))
                kept += 1
            if kept == 0:
                log.write(f"empty\t{m}\t{tgt['rank']}\t{tgt['name']}\tno usable sequence\n")

    rows.sort(key=lambda r: (r["marker"], r["accession"]))
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write("# Pinned reference library. One row per sequence; accession carries its version,\n")
        fh.write("# which is the pin, and md5 is of the TRIMMED uppercase sequence, so a silent\n")
        fh.write("# revision of a record fails the build. build_reference.py rebuilds from this file\n")
        fh.write("# and searches nothing.\n")
        fh.write(f"# resolved: {time.strftime('%Y-%m-%d')}   wall clock: {time.time() - t0:.0f} s\n")
        fh.write(f"# taxdump: {os.path.basename(a.taxdump)} md5 {td_md5}\n")
        fh.write(f"# barcode_uk ranges: {a.ranges}\n")
        for m in markers:
            lo, hi = SLEN[m]
            fh.write(f"# search {m}: txid<genus>[Organism:exp] AND ({MARKER_TERMS[m]}) AND {lo}:{hi}[SLEN],"
                     f" whole hit set paged, sorted locally\n")
            fh.write(f"# trim {m}: vsearch --usearch_global against one Barcode UK template per family, best hit; floor {FLOOR[m]} bp\n")
        fh.write(f"# caps: {a.max_species} sequences per species target, "
                 f"{a.max_genus_species} species per genus target, per marker\n")
        fh.write(f"# targets: {len(sp_t)} species-rank, {len(gn_t)} genus-rank, {len(rejects)} rejected\n")
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in COLS) + "\n")
    if a.vendor:
        # the manifest is the recipe; the sequences themselves are vendored so that a record
        # withdrawn from GenBank cannot make the repository unbuildable. NCBI places no
        # restrictions on redistributing GenBank data
        os.makedirs(a.vendor, exist_ok=True)
        for m in markers:
            with gzip.open(os.path.join(a.vendor, f"{m}.fasta.gz"), "wt") as fh:
                for r in rows:
                    if r["marker"] == m:
                        fh.write(f">{r['accession']}\n{r['_seq']}\n")
    with open(a.rejects, "w", encoding="utf-8") as fh:
        fh.write("# target names that produced no entry, with the reason. This file is the\n")
        fh.write("# denominator of every coverage figure quoted for the library.\n")
        fh.write("name\tsource\treason\n")
        for n, s, why in sorted(rejects):
            fh.write(f"{n}\t{s}\t{why}\n")
    log.write(f"# entrez: {ez.n_req} requests, {ez.n_hit} served from cache, {time.time() - t0:.0f} s\n")
    log.close()
    for m in markers:
        n = sum(1 for r in rows if r["marker"] == m)
        b = sum(1 for r in rows if r["marker"] == m and r["source"] == "barcode_uk")
        print(f"{m}: {n} sequences ({b} barcode_uk, {n - b} genbank)", file=sys.stderr)
    print(f"entrez: {ez.n_req} requests, {ez.n_hit} cached, {time.time() - t0:.0f} s", file=sys.stderr)


if __name__ == "__main__":
    main()
