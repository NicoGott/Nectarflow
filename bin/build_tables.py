#!/usr/bin/env python3
#
# Builds the read-count tables from the vsearch cluster map (.uc) and the
# taxonomic assignments: cluster x sample, taxon x sample, and a per-sample
# tracking line. Singleton clusters are dropped and counted separately.
#
import argparse, collections, csv, re

SIZE = re.compile(r";size=(\d+)")


#
# Splits a read label into its sample and its abundance. Labels look like
# SAMPLE|123;size=17 and are written by the per-sample dereplication
def label(lbl):
    sample = lbl.split("|", 1)[0]
    m = SIZE.search(lbl)
    return sample, int(m.group(1)) if m else 1


#
# Reads the cluster map, sums the reads per centroid and per sample, then
# writes the three output tables for this marker
def main():
    ap = argparse.ArgumentParser(description="Builds the cluster and taxon read-count tables for one marker")
    ap.add_argument("--uc", required=True)
    ap.add_argument("--assignments", required=True)
    ap.add_argument("--marker", required=True)
    a = ap.parse_args()

    counts = collections.defaultdict(lambda: collections.Counter())   # centroid -> sample -> reads
    for line in open(a.uc):
        f = line.rstrip("\n").split("\t")
        # S is the centroid of a new cluster, H a member of an existing one
        if f[0] == "S":
            centroid, member = f[8], f[8]
        elif f[0] == "H":
            centroid, member = f[9], f[8]
        else:
            continue
        sample, size = label(member)
        counts[re.sub(r";size=\d+;?", "", centroid)][sample] += size

    samples = sorted({s for c in counts.values() for s in c})
    total = {s: sum(c[s] for c in counts.values()) for s in samples}
    singles = {k for k, c in counts.items() if sum(c.values()) < 2}
    dropped = {s: sum(counts[k][s] for k in singles) for s in samples}
    for k in singles:
        del counts[k]

    assign = {re.sub(r";size=\d+;?", "", r["cluster"]): r for r in csv.DictReader(open(a.assignments), delimiter="\t")}
    taxa = collections.defaultdict(collections.Counter)
    for k, c in counts.items():
        r = assign.get(k)
        if r and r["rank"] != "discarded":
            for s, n in c.items():
                taxa[r["taxon"]][s] += n

    #
    # Writes one table: row key in the first column, one column per sample
    def write(path, table, first):
        with open(path, "w") as fh:
            fh.write(first + "\t" + "\t".join(samples) + "\n")
            for k in sorted(table):
                fh.write(k + "\t" + "\t".join(str(table[k][s]) for s in samples) + "\n")

    write(f"{a.marker}.cluster_table.tsv", counts, "cluster")
    write(f"{a.marker}.taxon_table.tsv", taxa, "taxon")
    with open(f"{a.marker}.table_stats.tsv", "w") as fh:
        fh.write(f"sample\t{a.marker}_clustered\t{a.marker}_singletons_dropped\t{a.marker}_assigned\n")
        for s in samples:
            fh.write(f"{s}\t{total[s]}\t{dropped[s]}\t{sum(t[s] for t in taxa.values())}\n")


if __name__ == "__main__":
    main()
