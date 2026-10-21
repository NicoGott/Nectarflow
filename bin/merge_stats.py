#!/usr/bin/env python3
#
# Joins the per-step read-count tables named on the command line into one
# tracking table on stdout: one row per sample, one column per stage, in the
# order the columns were first seen, plus a TOTAL row.
#
import csv, sys, collections

if "--help" in sys.argv or "-h" in sys.argv or len(sys.argv) == 1:
    sys.exit("usage: merge_stats.py STATS.tsv [STATS.tsv ...] > read_tracking.tsv\n\n"
             "Joins per-step read-count tables into one tracking table on stdout: one row per\n"
             "sample, one column per stage in the order first seen, plus a TOTAL row.")

table, cols = collections.defaultdict(dict), []
for path in sys.argv[1:]:
    for row in csv.DictReader(open(path), delimiter="\t"):
        s = row.pop("sample")
        for k, v in row.items():
            if k not in cols:
                cols.append(k)
            table[s][k] = v
w = csv.writer(sys.stdout, delimiter="\t")
w.writerow(["sample"] + cols)
for s in sorted(table):
    w.writerow([s] + [table[s].get(c, "") for c in cols])
# a stage a sample never reached leaves an empty cell, which counts as zero here
tot = ["TOTAL"] + [sum(int(float(table[s].get(c) or 0)) for s in table) for c in cols]
w.writerow(tot)
