#!/usr/bin/env python3
#
# Collects cutadapt's per-sample JSON reports into one table: pairs in, pairs
# kept per marker, and the pairs no marker claimed. The sample and the marker
# come from the report file name, which is <sample>.<marker>.*.json
#
import argparse, collections, json, pathlib


#
# Reads every JSON report and writes the one summary table
def main():
    ap = argparse.ArgumentParser(description="Collects the cutadapt JSON reports into one table")
    ap.add_argument("jsons", nargs="+")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    total, kept = {}, collections.defaultdict(dict)
    markers = []
    for f in a.jsons:
        sample, marker = pathlib.Path(f).name.split(".")[:2]
        d = json.load(open(f))
        # every marker's report for a sample saw the same input, so this just overwrites itself
        total[sample] = d["read_counts"]["input"]
        kept[sample][marker] = d["read_counts"]["output"]
        if marker not in markers:
            markers.append(marker)
    markers.sort()

    with open(a.out, "w") as fh:
        fh.write("sample\traw_pairs\t" + "\t".join(f"{m}_pairs" for m in markers) + "\tunassigned_pairs\n")
        for s in sorted(total):
            row = [kept[s].get(m, 0) for m in markers]
            fh.write(f"{s}\t{total[s]}\t" + "\t".join(str(v) for v in row) + f"\t{total[s] - sum(row)}\n")


if __name__ == "__main__":
    main()
