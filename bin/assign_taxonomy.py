#!/usr/bin/env python3
#
# Turns blastn hits into one taxonomic assignment per sequence cluster, the way
# Lowe et al. 2022 (J Appl Ecol, Methods 2.6) describe it, in place of the
# group's manual spreadsheet curation.
# Classification floors and the discard rules: docs/replica.md#7-blast-and-assign-blast_assign
#
import argparse, csv, collections


#
# Reads the BLAST table and the reference taxonomy and writes one row per
# cluster. The rank is one of:
# species
# genus
# genus_group
# family
# order
# discarded
def main():
    ap = argparse.ArgumentParser(description="Assigns a taxon to each cluster from its top BLAST hits",
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--blast", required=True, help="blastn -outfmt '6 std' (qseqid sseqid ... bitscore last)")
    ap.add_argument("--taxonomy", required=True, help="accession/species/genus/family/order TSV")
    ap.add_argument("--out", required=True)
    ap.add_argument("--percentile", type=float, default=1.0)
    ap.add_argument("--max-genera", type=int, default=3)
    a = ap.parse_args()

    tax = {r["accession"]: r for r in csv.DictReader(open(a.taxonomy), delimiter="\t")}
    hits = collections.defaultdict(list)
    seen = collections.defaultdict(set)
    for line in open(a.blast):
        f = line.rstrip("\n").split("\t")
        # one entry per reference sequence (a reference can return several local alignments)
        if f[1] in seen[f[0]] or len(hits[f[0]]) >= 20:
            continue
        seen[f[0]].add(f[1])
        hits[f[0]].append((f[1], float(f[11]), float(f[2])))

    tops = sorted(max(b for _, b, _ in h) for h in hits.values())
    cutoff = tops[int(len(tops) * a.percentile / 100.0)] if tops else 0.0

    with open(a.out, "w") as out:
        out.write("cluster\ttaxon\trank\ttop_bitscore\ttop_identity\tn_top_hits\tnote\n")
        for q, h in hits.items():
            best = max(b for _, b, _ in h)
            top = [(s, i) for s, b, i in h if b == best and s in tax]
            ident = max((i for _, i in top), default=0.0)
            if best < cutoff:
                out.write(f"{q}\tNA\tdiscarded\t{best}\t{ident}\t{len(top)}\tbelow bit-score percentile\n"); continue
            if not top:
                out.write(f"{q}\tNA\tdiscarded\t{best}\t{ident}\t0\treference without taxonomy\n"); continue
            rows = [tax[s] for s, _ in top]
            sp, ge = {r["species"] for r in rows}, sorted({r["genus"] for r in rows} - {""})
            fa = {r["family"] for r in rows}
            note = ""
            # only the hits sharing the exact top bit score vote; one species -> species, one genus
            # -> genus, one family -> family or "A/B/C spp." while at most --max-genera genera (as in
            # their published tables), several families -> discarded
            # a genus-rank reference row carries no species, and an empty name is not a call
            if len(sp) == 1 and "" not in sp:
                taxon, rank = sp.pop(), "species"
            elif len(ge) == 1:
                taxon, rank = f"{ge[0]} spp.", "genus"
            elif len(fa) == 1:
                if len(ge) <= a.max_genera:
                    taxon, rank = "/".join(ge) + " spp.", "genus_group"
                else:
                    taxon, rank = fa.pop(), "family"
            else:
                taxon, rank, note = "NA", "discarded", "top hits span several families"
            out.write(f"{q}\t{taxon}\t{rank}\t{best}\t{ident}\t{len(top)}\t{note}\n")


if __name__ == "__main__":
    main()
