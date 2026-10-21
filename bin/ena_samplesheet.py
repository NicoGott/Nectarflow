#!/usr/bin/env python3
#
# Builds the pipeline samplesheet (sample, run, seq_run, fastq_1, fastq_2,
# md5_1, md5_2, read_count) from an ENA or SRA project accession, by querying
# the ENA file report API. Single-end runs are skipped.
#
import argparse, csv, re, sys, urllib.request

URL = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}&result=read_run"
       "&fields=run_accession,sample_alias,library_layout,read_count,fastq_ftp,fastq_md5,library_name&format=tsv")


#
# Downloads the run table for the accession and writes one CSV row per paired
# run, sorted by sample alias
def main():
    ap = argparse.ArgumentParser(description="Builds a samplesheet from an ENA/SRA project accession")
    ap.add_argument("accession")
    ap.add_argument("-o", "--out", default="samplesheet.csv")
    ap.add_argument("--include", default="", help="comma-separated sample aliases to keep (default: all)")
    a = ap.parse_args()
    with urllib.request.urlopen(URL.format(acc=a.accession), timeout=120) as r:
        rows = list(csv.DictReader(r.read().decode().splitlines(), delimiter="\t"))
    keep = {s.strip() for s in a.include.split(",") if s.strip()}
    n, seen = 0, set()
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "run", "seq_run", "fastq_1", "fastq_2", "md5_1", "md5_2", "read_count"])
        for r in sorted(rows, key=lambda x: x["sample_alias"]):
            files, md5s = r["fastq_ftp"].split(";"), r["fastq_md5"].split(";")
            if r["library_layout"] != "PAIRED" or len(files) < 2:
                print(f"skip {r['run_accession']}: not paired", file=sys.stderr)
                continue
            if keep and r["sample_alias"] not in keep:
                continue
            # sample aliases become file names and vsearch labels: keep them shell-safe
            sample = r["sample_alias"].strip().replace(" ", "_").replace("|", "_")
            if sample in seen:      # same sample sequenced twice: keep both runs, make the name unique
                sample = f"{sample}_{r['run_accession']}"
            seen.add(sample)
            # sequencing run, out of library_name (e.g. "..._run_2"): the second arm learns DADA2
            # error models per run, the faithful arm ignores the column
            m = re.search(r"(run[_-]?\d+)$", (r.get("library_name") or "").strip())
            # the last two files are the R1/R2 pair; an index read, when present, comes first
            w.writerow([sample, r["run_accession"], m.group(1) if m else "run_1",
                        "https://" + files[-2], "https://" + files[-1],
                        md5s[-2], md5s[-1], r["read_count"]])
            n += 1
    print(f"{n} paired runs written to {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
