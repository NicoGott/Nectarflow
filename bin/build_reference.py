#!/usr/bin/env python3
#
# Builds the reference library from assets/reference_manifest.tsv. It searches
# nothing: every accession, every trim coordinate and every taxonomic label is
# already fixed in the manifest, so the library does not depend on what NCBI
# returns on the day. Writes <marker>.fasta and <marker>.taxonomy.tsv, in
# manifest order, for makeblastdb to index.
#
# By default the sequences come from the vendored copy in assets/reference/.
# GenBank places no restrictions on redistributing its data, and a manifest is
# a recipe rather than an artefact: a record withdrawn in five years would
# otherwise leave the repository unbuildable. --refetch pulls the records from
# GenBank again and reports every difference, which is the audit path.
#
import argparse, csv, gzip, hashlib, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
TAX = ("accession", "species", "genus", "family", "order", "lineage")


#
# Reads the manifest, skipping its comment header
def read_manifest(path):
    return list(csv.DictReader((l for l in open(path, encoding="utf-8") if not l.startswith("#")), delimiter="\t"))


#
# Reads a FASTA, plain or gzipped, into {accession: sequence}
def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    out, name, seq = {}, None, []
    with op(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = "".join(seq)
                name, seq = line[1:].split()[0], []
            elif line.strip():
                seq.append(line.strip().upper())
    if name:
        out[name] = "".join(seq)
    return out


#
# Pulls the records named in the manifest from GenBank and cuts each one to the
# coordinates the manifest records. epost truncates long id lists without
# saying so, so the chunk is small and anything that did not come back is asked
# for again on its own rather than by re-posting the whole chunk
def refetch(rows, email, cache):
    if not email:
        raise SystemExit("build_reference: --refetch needs --email")
    from entrez import Entrez
    ez = Entrez(email=email, cache=cache)
    accs = sorted({r["accession"] for r in rows})
    raw = ez.fetch_fasta(accs)
    missing = [a for a in accs if a not in raw]
    for a in missing:
        raw.update(ez.fetch_fasta([a], chunk=1))
    still = [a for a in accs if a not in raw]
    if still:
        raise SystemExit(f"build_reference: {len(still)} accessions were not returned, first {still[:5]}")
    return raw


#
# Cuts one record the way the manifest says
def cut(row, seq):
    s = seq.translate(CMP)[::-1] if row["strand"] == "-" else seq
    return s[int(row["start"]) - 1:int(row["end"])]


#
# Checks a finished sequence against the manifest. A mismatch is a silently
# revised GenBank record, or a vendored copy that is not what was resolved, and
# it fails the build by name
def check(row, s):
    md5 = hashlib.md5(s.encode()).hexdigest()
    if md5 != row["md5"]:
        raise SystemExit(f"build_reference: {row['accession']} ({row['marker']}) does not match the manifest: "
                         f"md5 {md5} for {len(s)} bp, manifest says {row['md5']} for {row['length']} bp")
    return s


#
# Writes one FASTA and one taxonomy table per marker
def main():
    ap = argparse.ArgumentParser(description="Builds the reference library from the pinned manifest")
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--vendored", default="", help="directory of trimmed <marker>.fasta.gz; default assets/reference")
    ap.add_argument("--refetch", action="store_true", help="pull the records from GenBank again and diff")
    ap.add_argument("--email", default="", help="required with --refetch")
    ap.add_argument("--cache", default=".entrez_cache")
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    rows = read_manifest(a.manifest)
    markers = sorted({r["marker"] for r in rows})
    vend = a.vendored or os.path.join(os.path.dirname(os.path.abspath(a.manifest)), "reference")

    raw = refetch(rows, a.email, a.cache) if a.refetch else {}
    for m in markers:
        mine = [r for r in rows if r["marker"] == m]
        entries, diffs = [], 0
        have = {} if a.refetch else read_fasta(os.path.join(vend, f"{m}.fasta.gz"))
        for r in mine:
            acc = r["accession"]
            if a.refetch:
                if acc not in raw:
                    raise SystemExit(f"build_reference: {acc} ({m}) was not returned by GenBank")
                s = check(r, cut(r, raw[acc]))
            else:
                if acc not in have:
                    raise SystemExit(f"build_reference: {acc} ({m}) is not in the vendored FASTA")
                s = check(r, have[acc])
            entries.append((acc, s, r))
        if a.refetch and os.path.exists(os.path.join(vend, f"{m}.fasta.gz")):
            old = read_fasta(os.path.join(vend, f"{m}.fasta.gz"))
            diffs = sum(1 for acc, s, _ in entries if old.get(acc) != s)
            print(f"{m}: refetch audit, {diffs} of {len(entries)} differ from the vendored copy", file=sys.stderr)
        n_in = len(entries)
        with open(os.path.join(a.outdir, f"{m}.fasta"), "w") as fa, \
             open(os.path.join(a.outdir, f"{m}.taxonomy.tsv"), "w") as tx:
            tx.write("\t".join(TAX) + "\n")
            for acc, s, r in entries:
                fa.write(f">{acc}\n{s}\n")
                tx.write(acc + "\t" + "\t".join(r[k] for k in TAX[1:]) + "\n")
        bu = sum(1 for r in mine if r["source"] == "barcode_uk")
        print(f"{m}: {len(entries)} sequences written ({n_in} in the manifest, {bu} barcode_uk)", file=sys.stderr)


if __name__ == "__main__":
    main()
