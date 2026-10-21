#!/usr/bin/env python3
"""Pull the per-insect metadata of Lowe et al. 2022 out of the Dryad workbook, as a plain TSV.

Written for the one column a sample code cannot carry: which insect the pollen came off.
See docs/reanalysis.md, "Phase 6 in detail: the ecology script"."""
import argparse, csv, datetime, re

SHEET = "Pollen_sample_metadata"                 # the only sheet this script may open
HEADER_ROWS = 2                                  # title row, then the column names
# our name -> column name in their sheet. The columns ecology.R also derives from the sample name
# are taken too, so the two derivations can be cross-checked against each other.
COLUMNS = {
    "sample_published": "Sample Code", "season": "Season", "month": "Month", "year": "Year",
    "date": "Date", "habitat": "Habitat", "site": "Site", "transect": "Transect",
    "insect_group": "Pollinator Group", "insect_binomial": "Pollinator Binomial",
    "insect_order": "Order", "insect_family": "Family", "sex": "Sex",
    "observed_on": "Observed_on",
}
NAME = re.compile(r"^([A-Z]{2,4})([0-9]+)_([0-9]{1,2})([0-9]{2})_([0-9]+)$")


def norm(code):
    # identical rule to compare_published.norm; duplicated rather than imported because Nextflow
    # stages each bin/ script on its own and a cross-import would not resolve in the work directory
    return re.sub(r"_0*(\d+)$", r"_\1", str(code).strip())


def read_seq_runs(path):
    """sample -> seq_run from the samplesheet. A sample split over several rows is one run; if the
    rows disagree the sample gets no run rather than an arbitrary one of them."""
    seen = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            if "seq_run" not in row:
                return {}
            s, run = norm(row.get("sample", "")), (row.get("seq_run") or "").strip()
            if not s or not run:
                continue
            if s in seen and seen[s] != run:
                seen[s] = ""
            else:
                seen.setdefault(s, run)
    return {s: r for s, r in seen.items() if r}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--published", required=True, help="Supporting_Data.xlsx (Dryad)")
    ap.add_argument("--samplesheet", help="the run samplesheet, for the seq_run column")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    import openpyxl

    wb = openpyxl.load_workbook(a.published, read_only=True)
    rows = list(wb[SHEET].iter_rows(values_only=True))
    head = [str(c).strip() if c is not None else "" for c in rows[HEADER_ROWS - 1]]
    missing = [v for v in COLUMNS.values() if v not in head]
    if missing:
        raise SystemExit(f"{a.published}: sheet {SHEET} has no column(s) {missing}")
    idx = {k: head.index(v) for k, v in COLUMNS.items()}

    out = []
    for r in rows[HEADER_ROWS:]:
        if not r or not r[0]:
            continue
        rec = {k: r[i] for k, i in idx.items()}
        rec["sample"] = norm(rec["sample_published"])
        if isinstance(rec["date"], datetime.datetime):
            rec["date"] = rec["date"].date().isoformat()
        for k, v in rec.items():
            rec[k] = "" if v is None else str(v).strip()
        m = NAME.match(rec["sample"])
        # a name that does not parse would silently lose its derived columns in ecology.R, so say so
        rec["name_parses"] = "yes" if m else "no"
        out.append(rec)
    out.sort(key=lambda d: d["sample"])

    # Which sequencing run each insect was on. It is a technical fact, not an ecological one, but
    # it belongs here: a study variable that turns out to be a function of the run cannot be told
    # apart from it, and ecology.R can only say so about a column it can see.
    batch = read_seq_runs(a.samplesheet) if a.samplesheet else {}
    for rec in out:
        rec["seq_run"] = batch.get(rec["sample"], "")
    n_batch = sum(1 for rec in out if rec["seq_run"])

    fields = ["sample", "sample_published", "name_parses"] + [k for k in COLUMNS if k != "sample_published"]
    if a.samplesheet:
        fields.append("seq_run")
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(out)

    print(f"{len(out)} insects written to {a.out}; "
          f"{sum(d['name_parses'] == 'no' for d in out)} sample codes do not parse"
          + (f"; sequencing run known for {n_batch}" if a.samplesheet else ""))


if __name__ == "__main__":
    main()
