#!/usr/bin/env python3
#
# Classifies one marker's ASVs (arm 2, phase 3) from the blastn hits against the rebuilt
# reference library. Lowest common ancestor over the hits within a substitution of the best
# raw score; every rule may only make a call coarser.
# Floors, band and the rest: docs/reanalysis.md#phase-3-in-detail-classification
#
import argparse, collections, csv, os, sys

# From the published discrimination of the two markers, NOT a calibration run here; see
# docs/reanalysis.md#where-the-identity-floors-come-from.
PARAMS = {
    "rbcL": dict(species=100.0, genus=99.5, family=99.0, order=98.0,
                 min_coverage=80.0, max_evalue=1e-20, min_asv_len=400),
    "ITS2": dict(species=100.0, genus=95.0, family=92.0, order=90.0,
                 min_coverage=80.0, max_evalue=1e-20, min_asv_len=200),
}
BAND_SUBSTITUTIONS, MATCH_REWARD, MISMATCH_PENALTY = 1, 1, 2
MAX_GENERA, MAX_TARGET_SEQS = 3, 500

# What BLAST_ASV asks blastn for. qlen, qcovhsp and score are not optional: coverage and the tie
# band are computed from them, and a run without them would skip two gates in silence.
FIELDS = "qseqid sseqid pident length qlen qcovhsp mismatch gapopen evalue bitscore score"

RANKS = ("species", "genus", "family", "order")
STRICT = {"species": "species", "genus": "genus", "genus_group": "family", "family": "family",
          "order": "order", "order_spanning": "none", "unassigned": "none", "asv_too_short": "none"}
BUCKET = {"no_hit": "__unassigned_no_hit", "below_identity": "__unassigned_below_identity",
          "below_coverage": "__unassigned_below_coverage",
          "reference_without_taxonomy": "__unassigned_no_taxonomy",
          "order_spanning": "__order_spanning", "asv_too_short": "__asv_too_short"}
COLUMNS = ("asv taxon rank top_bitscore top_identity n_top_hits note "           # arm 1's header,
           "rank_strict top_coverage genus_singleton_in_library no_exact_match").split()

Hit = collections.namedtuple("Hit", "acc key bitscore pident coverage")


#
# A read count from the ASV table as an int, refusing anything that is not a whole number
def count(v):
    # R writes a round count of 1e5 or more as '1e+05', which plain int() cannot read.
    n = float(v)
    if n != int(n):
        sys.exit(f"ASV table holds a non-integer count: {v!r}")
    return int(n)


#
# Reads the blastn table and keeps, per query, the best line for each subject. Returns those
# hits, the query lengths, the number of short lines seen, and how many hits each query lost
# to the e-value gate.
def read_blast(path, max_evalue):
    idx = {name: i for i, name in enumerate(FIELDS.split())}
    hits, qlen, short = collections.defaultdict(dict), {}, 0
    dropped = collections.Counter()
    for line in open(path, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) < len(idx):
            # A half-written last line means a truncated search; counted, and the caller
            # refuses to run rather than classify it as a complete one.
            if line.strip():
                short += 1
            continue
        # Counted per query as well as dropped, so an ASV whose only hits fail the gate is not
        # then reported as no_hit.
        if float(f[idx["evalue"]]) > max_evalue:
            dropped[f[idx["qseqid"]]] += 1
            continue
        q, acc = f[idx["qseqid"]], f[idx["sseqid"]]
        bits, pident = float(f[idx["bitscore"]]), float(f[idx["pident"]])
        qlen[q] = int(f[idx["qlen"]])
        cov = float(f[idx["qcovhsp"]])
        h = Hit(acc, float(f[idx["score"]]), bits, pident, cov)
        old = hits[q].get(acc)
        # Coverage is part of the tie-break: two HSPs of one subject can share score and
        # identity yet cover different fractions, and "first line wins" would make the call
        # depend on the order blastn printed them in.
        if old is None or (h.key, h.pident, h.coverage) > (old.key, old.pident, old.coverage):
            hits[q][acc] = h
    return hits, qlen, short, dropped


#
# Sequence length per name in a fasta, so an ASV with no hit still has a length
def fasta_lengths(path):
    lengths, name = {}, None
    for line in open(path, encoding="utf-8"):
        if line.startswith(">"):
            name = line[1:].split()[0]
            lengths[name] = 0
        elif name:
            lengths[name] += len(line.strip())
    return lengths


#
# The label and the rank for a lowest common ancestor, which can be:
# species
# genus
# genus_group, i.e. "A/B/C spp." for a few genera of one family
# family
# order
# order_spanning, when the candidates do not share one order
def label_for(lca, rows, genera, max_genera):
    if lca == "species":
        return rows[0]["species"], "species", ""
    if lca == "genus":
        return f"{genera[0]} spp.", "genus", ""
    if lca == "family":
        if 0 < len(genera) <= max_genera:               # "A/B/C spp.", as in their published tables
            return "/".join(genera) + " spp.", "genus_group", ""
        # Past that many genera the group label stops being a name and becomes a list, so the call
        # is the family. Say which genera were dropped: this is where most of the resolution is
        # lost against the published tables, and without the note it is lost silently too. A wider
        # reference library ties more genera here, so the same rule bites harder than it did on
        # Barcode UK alone.
        if len(genera) > max_genera:
            return rows[0]["family"], "family", (f"{len(genera)} genera tie above --max-genera "
                                                 f"{max_genera}: " + "/".join(genera))
        return rows[0]["family"], "family", "no genus in the reference taxonomy for these hits"
    if lca == "order":
        return rows[0]["order"], "order", "several families in one order; rule not stated in the paper"
    return "NA", "order_spanning", "candidates span several orders"


#
# One step down the rank ladder, and the label that goes with it. A genus call demotes to the
# single FAMILY, not to a one-genus group, which would claim the same genus on family-level
# evidence; genus_group and family share a floor, so both of those drop to order. An empty
# field is never agreement, so a row with no genus falls through to unassigned.
def coarser(rank, rows, genera):
    families = sorted({r["family"] for r in rows if r["family"]})
    orders = sorted({r["order"] for r in rows if r["order"]})
    if rank == "species" and genera:
        return "genus", f"{genera[0]} spp."
    if rank in ("species", "genus") and len(families) == 1:
        return "family", families[0]
    if rank in ("species", "genus", "genus_group", "family") and len(orders) == 1:
        return "order", orders[0]
    return "unassigned", "NA"


#
# Classifies one ASV. Returns its output row and the bucket its reads are accounted for in,
# or None when the ASV was assigned. It can fail before the hits are weighed (too short, no
# hit) or after (below the identity or the coverage floor, reference without taxonomy).
def classify(asv, subjects, length, tax, genus_size, a, floor, n_over_evalue=0):
    row = dict(asv=asv, taxon="NA", rank="unassigned", top_bitscore="0.0", top_identity="0.0",
               n_top_hits=0, note="", rank_strict="none", top_coverage="NA",
               genus_singleton_in_library="NA", no_exact_match="NA",
               )

    # An ASV much shorter than the marker's insert is a phase-2 artefact, not a plant; reported
    # rather than classified, and its reads stay in the accounting.
    if length is not None and length < a.min_asv_len:
        row.update(rank="asv_too_short", note=f"asv length {length} < {a.min_asv_len}")
        return row, BUCKET["asv_too_short"]
    if not subjects:
        row["note"] = ("no_hit" if not n_over_evalue else
                       f"no_hit; {n_over_evalue} hit(s) removed by --max-evalue {a.max_evalue:g}")
        return row, BUCKET["no_hit"]

    ordered = sorted(subjects.values(), key=lambda h: (-h.key, h.acc))   # deterministic order
    best_raw = ordered[0]
    row.update(top_bitscore=f"{best_raw.bitscore:.1f}", top_identity=f"{best_raw.pident:.3f}",
               top_coverage=f"{best_raw.coverage:.1f}")

    # Admissible: above the order floor, covered, and in the taxonomy.
    adm = [h for h in ordered if h.pident >= floor["order"]
           and h.coverage >= a.min_coverage and h.acc in tax]
    if not adm:
        if best_raw.pident < floor["order"]:
            reason = "below_identity"
        elif best_raw.coverage is not None and best_raw.coverage < a.min_coverage:
            reason = "below_coverage"
        else:
            reason = "reference_without_taxonomy"
        row["note"] = reason
        return row, BUCKET[reason]

    # Candidates within the tie band, then their lowest common ancestor.
    cand = [h for h in adm if h.key >= adm[0].key - a.band_points]
    rows = [tax[h.acc] for h in cand]
    genera = sorted({r["genus"] for r in rows if r["genus"]})
    lca = next((r for r in RANKS
                if len({t[r] for t in rows}) == 1 and rows[0][r]), None)
    ident = max(h.pident for h in cand)                 # the best any candidate reaches
    notes = []

    taxon, rank, note = label_for(lca, rows, genera, a.max_genera)
    fail = BUCKET["below_identity"]          # where an unassigned row's reads are accounted for
    if rank == "order_spanning" and len({r["order"] for r in rows if r["order"]}) < 2:
        # The candidates do not DISAGREE about the order: the reference records none for them,
        # and phases 4-5 read __order_spanning as a tie across REAL orders.
        rank, taxon, note = "unassigned", "NA", "no order in the reference taxonomy"
        fail = BUCKET["reference_without_taxonomy"]
    if note:
        notes.append(note)
    # Demotion only: a floor says how far the LIBRARY can be trusted at that rank and can never
    # introduce a name, so every step is a loss of resolution and is recorded.
    while rank in ("species", "genus", "genus_group", "family", "order"):
        limit = floor["family"] if rank == "genus_group" else floor[rank]
        if ident >= limit:
            break
        was = rank
        rank, taxon = coarser(rank, rows, genera)
        notes.append(f"demoted {was}->{rank} (best identity {ident:.2f} < {limit})")
    if rank == "unassigned" and fail == BUCKET["below_identity"]:
        notes.append("below_identity")

    # "no" must mean an exact, FULL-LENGTH match to ONE reference: identity covers the aligned
    # part only, so maxima taken separately would invent a match neither reference is.
    exact = any(h.pident >= 100.0 and h.coverage >= 100.0 for h in cand)
    # The three top_ columns describe ONE reference, the best candidate: maxima taken separately
    # over the candidate set can come from different subjects. A row that never gets a candidate
    # set keeps the best HIT's three instead, or zeros when there was no hit at all.
    row.update(taxon=taxon, rank=rank, rank_strict=STRICT[rank],
               top_bitscore=f"{cand[0].bitscore:.1f}", top_identity=f"{cand[0].pident:.3f}",
               top_coverage=f"{cand[0].coverage:.1f}",
               n_top_hits=len(cand), no_exact_match="no" if exact else "yes",
               note="; ".join(notes))
    if rank in ("species", "genus"):
        g = genera[0] if genera else None
        row["genus_singleton_in_library"] = "yes" if genus_size.get(g, 0) == 1 else "no"
    if rank == "order_spanning":
        return row, BUCKET["order_spanning"]
    if rank == "unassigned":
        return row, fail
    return row, None


#
# Reads the taxonomy, the ASV table and the BLAST hits, classifies every ASV, and writes the
# per-ASV table, the two taxon tables, the rank summary and the thresholds that were in force
def main():
    ap = argparse.ArgumentParser(
            description="""Classifies a marker's ASVs from blastn hits against the reference library""",
            formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--blast", required=True, help=f"blastn -outfmt '6 {FIELDS}'")
    ap.add_argument("--taxonomy", required=True, help="<marker>.taxonomy.tsv (accession species genus family order lineage)")
    ap.add_argument("--asv-table", required=True, help="phase-2 <marker>.asv_table.tsv, samples in rows")
    ap.add_argument("--marker", required=True, choices=sorted(PARAMS))
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--asv-fasta", help="<marker>.asvs.fasta, so ASVs without a hit also get a length")
    # Thresholds default to the marker's calibrated value, so one flag is one sensitivity run.
    for rank in RANKS:
        ap.add_argument(f"--{rank}-floor", type=float, dest=f"floor_{rank}",
                        help=f"minimum identity for a {rank}-level call")
    ap.add_argument("--min-coverage", type=float, help="minimum qcovhsp of an admissible hit")
    ap.add_argument("--max-evalue", type=float)
    ap.add_argument("--min-asv-len", type=int)
    ap.add_argument("--band-substitutions", type=int, default=BAND_SUBSTITUTIONS,
                    help="how far below the best raw score a hit still votes, in substitutions")
    ap.add_argument("--match-reward", type=int, default=MATCH_REWARD)
    ap.add_argument("--mismatch-penalty", type=int, default=MISMATCH_PENALTY)
    ap.add_argument("--max-genera", type=int, default=MAX_GENERA,
                    help="genera still written as 'A/B/C spp.'; every published group label has <= 3")
    ap.add_argument("--max-target-seqs", type=int, default=MAX_TARGET_SEQS, help="as given to blastn")
    a = ap.parse_args()

    for key, value in PARAMS[a.marker].items():
        name = f"floor_{key}" if key in RANKS else key
        if getattr(a, name) is None:
            setattr(a, name, value)
    floor = {r: getattr(a, f"floor_{r}") for r in RANKS}
    # blastn writes its own penalty negative (-penalty -2); taken literally the band would be
    # NEGATIVE and every candidate set empty
    if a.band_substitutions < 0 or a.match_reward <= 0 or a.mismatch_penalty <= 0:
        sys.exit("--band-substitutions must be >= 0 and --match-reward/--mismatch-penalty > 0; the "
                 "mismatch penalty is a magnitude, so blastn's '-penalty -2' is '--mismatch-penalty 2'")
    a.band_points = a.band_substitutions * (a.match_reward + a.mismatch_penalty)
    tax_rows = csv.DictReader(open(a.taxonomy, encoding="utf-8"), delimiter="\t")
    missing = [c for c in ("accession",) + RANKS if c not in (tax_rows.fieldnames or [])]
    if missing:
        sys.exit(f"{a.taxonomy} lacks column(s): {', '.join(missing)}")
    tax = {r["accession"]: r for r in tax_rows}
    genus_size = collections.Counter()
    for genus, species in {(r["genus"], r["species"]) for r in tax.values() if r["genus"]}:
        genus_size[genus] += 1

    table = [r for r in csv.reader(open(a.asv_table, encoding="utf-8"), delimiter="\t") if r]
    names = table[0][1:]
    # A name used twice stops the table being a matrix, and the read accounting cannot see it.
    for what, repeated in (("sample(s)", collections.Counter(r[0] for r in table[1:])),
                           ("ASV column(s)", collections.Counter(names))):
        twice = sorted(n for n, k in repeated.items() if k > 1)
        if twice:
            sys.exit(f"{a.asv_table}: {what} listed more than once: {', '.join(twice)}")
    samples = sorted({r[0] for r in table[1:]})
    reads = collections.defaultdict(collections.Counter)                 # asv -> sample -> reads
    for r in table[1:]:
        for name, v in zip(names, r[1:]):
            if v and v != "0":
                reads[name][r[0]] += count(v)

    hits, qlen, short, over_evalue = read_blast(a.blast, a.max_evalue)
    if short:
        sys.exit(f"{a.blast}: {short} line(s) with fewer than {len(fields.split())} fields; the "
                 f"search looks truncated, and classifying it would look like a complete run")
    lengths = fasta_lengths(a.asv_fasta) if a.asv_fasta else qlen
    # the cap the search was given is a separate parameter in nextflow.config and is now passed in,
    # so this only catches a hand-run search that used a larger one
    observed = max((len(s) for s in hits.values()), default=0)
    if observed > a.max_target_seqs:
        sys.exit(f"{a.blast}: a query has {observed} subjects but --max-target-seqs is "
                 f"{a.max_target_seqs}; pass the value the search actually used.")

    os.makedirs(a.out, exist_ok=True)
    taxa = collections.defaultdict(collections.Counter)
    buckets = collections.defaultdict(collections.Counter)
    # An ASV in the fasta that reached neither the table header nor the search is a phase-2
    # bug, and must show up as a row rather than as nothing at all.
    asvs = sorted(set(names) | set(hits) | set(lengths))                 # deterministic order
    n_asv, n_read = collections.Counter(), collections.Counter()         # (level, value) -> count
    with open(os.path.join(a.out, f"{a.marker}.asv_taxonomy.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for asv in asvs:
            row, bucket = classify(asv, hits.get(asv, {}), lengths.get(asv), tax, genus_size,
                                   a, floor, over_evalue.get(asv, 0))
            fh.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")
            (buckets[bucket] if bucket else taxa[row["taxon"]]).update(reads[asv])
            # An ASV-level percentage hides abundance and a read-level one hides diversity, so
            # every line of the summary carries both
            keys = [("rank", row["rank"]), ("rank_strict", row["rank_strict"])]
            for flag in ("no_exact_match", "genus_singleton_in_library"):
                if row[flag] == "yes":
                    keys.append(("flag", flag))
            if row["top_coverage"] != "NA" and float(row["top_coverage"]) < 90.0 and bucket is None:
                keys.append(("flag", "top_coverage_below_90"))
            for k in keys:
                n_asv[k] += 1
                n_read[k] += sum(reads[asv].values())

    #
    # Writes a taxon x sample matrix
    def write(path, rows):
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("taxon\t" + "\t".join(samples) + "\n")
            for k in sorted(rows):
                fh.write(k + "\t" + "\t".join(str(rows[k][s]) for s in samples) + "\n")

    write(os.path.join(a.out, f"{a.marker}.taxon_table.tsv"), taxa)
    write(os.path.join(a.out, f"{a.marker}.taxon_table.full.tsv"), {**taxa, **buckets})

    # The read accounting must close, or a table would quietly carry fewer reads than phase 2.
    for s in samples:
        got = sum(t[s] for t in taxa.values()) + sum(b[s] for b in buckets.values())
        want = sum(reads[n][s] for n in names)
        if got != want:
            sys.exit(f"{a.marker} {s}: {got} reads accounted for, {want} in the ASV table")

    total_reads = sum(sum(c.values()) for c in reads.values())
    assigned = sum(sum(t.values()) for t in taxa.values())
    with open(os.path.join(a.out, f"{a.marker}.rank_summary.tsv"), "w", encoding="utf-8") as fh:
        fh.write("marker\tlevel\tvalue\tasvs\treads\tpct_reads\n")
        for level, value in sorted(n_asv):
            r = n_read[(level, value)]
            fh.write(f"{a.marker}\t{level}\t{value}\t{n_asv[(level, value)]}\t{r}\t"
                     f"{100.0 * r / total_reads if total_reads else 0.0:.2f}\n")

    # The thresholds that were in force, so a table cannot be read under cut-offs it was not
    # made with. A gate the input could not support is written 'not applied', never as its
    # nominal value.
    params = [("marker", a.marker)]
    params += [(f"{r}_floor", floor[r]) for r in RANKS]
    params += [
        ("identity_gate_for_admissibility", floor["order"]),
        ("min_coverage", a.min_coverage),
        ("max_evalue", f"{a.max_evalue:g}"),
        ("evalue_applied_by", "this script"),
        ("hits_removed_by_max_evalue", sum(over_evalue.values())),
        ("min_asv_len", a.min_asv_len),
        ("band_substitutions", a.band_substitutions if a.band_points or a.band_substitutions == 0
         else "not applied"),
        ("band_raw_points", a.band_points),
        ("band_basis", "raw score"),
        ("match_reward", a.match_reward), ("mismatch_penalty", a.mismatch_penalty),
        ("max_genera", a.max_genera),
        ("max_target_seqs_declared", a.max_target_seqs),
        ("max_subjects_per_query_observed", observed),
        ("blast_fields", FIELDS),
        ("n_samples", len(samples)), ("n_asvs", len(asvs)), ("n_taxa", len(taxa)),
        ("reads_total", total_reads), ("reads_assigned", assigned),
        ("reads_unassigned", total_reads - assigned),
    ]
    with open(os.path.join(a.out, f"{a.marker}.classify_params.tsv"), "w", encoding="utf-8") as fh:
        fh.write("parameter\tvalue\n")
        for k, v in params:
            fh.write(f"{k}\t{v}\n")

    print(f"{a.marker}: {len(asvs)} ASVs, {len(taxa)} taxa, {assigned} reads assigned, "
          f"{total_reads - assigned} unassigned", file=sys.stderr)


if __name__ == "__main__":
    main()
