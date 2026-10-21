#!/usr/bin/env python3
#
# Writes one sheet of the published workbook out as a plain taxon x sample TSV,
# so phase 6 can put the authors' table through the same analysis script as
# ours. Read counts only: nothing is renamed, dropped or rescaled.
#
import argparse, csv, re


#
# Reads the sheet, drops its leading metadata columns and transposes it: the
# sheet has one row per sample, the output one row per taxon
def main():
    ap = argparse.ArgumentParser(description="Writes a sheet of the published workbook out as a taxon x sample TSV")
    ap.add_argument("--published", required=True)
    ap.add_argument("--sheet", default="rbcL_ITS2_consensus_read_counts")
    ap.add_argument("--meta-cols", type=int, default=8, help="leading metadata columns of the sheet")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    import openpyxl
    # row 0 is a banner, row 1 the taxon names, the samples start at row 2
    rows = list(openpyxl.load_workbook(a.published, read_only=True)[a.sheet].iter_rows(values_only=True))
    taxa = [str(t).strip() for t in rows[1][a.meta_cols:]]
    samples, counts = [], []
    for r in rows[2:]:
        if not r[0]:
            continue
        samples.append(re.sub(r"_0*(\d+)$", r"_\1", str(r[0]).strip()))     # as compare_published.norm
        counts.append([int(float(v or 0)) for v in r[a.meta_cols:]])

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["taxon"] + samples)
        for j, taxon in enumerate(taxa):
            w.writerow([taxon] + [row[j] for row in counts])
    print(f"{len(taxa)} taxa x {len(samples)} samples from sheet {a.sheet}")


if __name__ == "__main__":
    main()
