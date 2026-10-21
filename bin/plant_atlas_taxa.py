#!/usr/bin/env python3
#
# Writes assets/plant_atlas_2020_gb.tsv: every vascular plant taxon the BSBI
# records in Great Britain, with its native status. Source is the Plant Atlas
# 2020 native-status dataset, Stroh, Walker, Humphrey, Pescott & Burkmar 2024,
# doi:10.5281/zenodo.11097371, CC BY 4.0 - the Botanical Society's own successor
# to Preston, Pearman & Dines (2002), with nomenclature from the BSBI
# Distribution Database, which tracks Stace (2019).
#
# Every filter below is a judgement call, so each is one commented step and the
# counts are printed. Names are written exactly as the source has them; nothing
# is normalised here. The dataset covers Britain, Ireland and the Channel
# Islands, so reading statusGB rather than overallStatus is load-bearing.
#
import argparse, csv, hashlib, io, sys, urllib.request, zipfile

ZENODO = "https://zenodo.org/records/11097371/files/PA2020_NativeStatusData.zip?download=1"
MD5 = "a5f60929e6118b77094110d6ca77bb48"
ROWS = 3497
CITE = ("Stroh, P.A., Walker, K.J., Humphrey, T.A., Pescott, O.L. & Burkmar, R.J. (2024) "
        "Plant Atlas 2020: native status data. Zenodo. https://doi.org/10.5281/zenodo.11097371")


#
# Returns the zip as bytes, from disk or from Zenodo, and refuses to go on if
# the md5 is not the one Zenodo publishes for this version
def load(path):
    if path:
        blob = open(path, "rb").read()
    else:
        with urllib.request.urlopen(ZENODO, timeout=300) as r:
            blob = r.read()
    md5 = hashlib.md5(blob).hexdigest()
    if md5 != MD5:
        raise SystemExit(f"plant_atlas_taxa: md5 {md5} is not the published {MD5}")
    return blob, md5


#
# Joins status.csv to nameLookup.csv on the record id and applies the filters.
# Returns the kept rows and a count per statusGB value
def derive(blob):
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        root = next(n.rsplit("/", 1)[0] for n in z.namelist() if n.endswith("status.csv"))
        status = list(csv.DictReader(io.TextIOWrapper(z.open(f"{root}/status.csv"), "utf-8")))
        names = list(csv.DictReader(io.TextIOWrapper(z.open(f"{root}/nameLookup.csv"), "utf-8")))
    if len(status) != ROWS or len(names) != ROWS:
        raise SystemExit(f"plant_atlas_taxa: expected {ROWS} rows, got {len(status)} and {len(names)}")
    # 1. join: status.ddbid is nameLookup.ddbidOrig (ddbidLookup is the same id with . replaced by _)
    by_id = {r["ddbidOrig"]: r for r in names}
    counts, out = {}, []
    for s in status:
        n = by_id.get(s["ddbid"])
        if n is None:
            raise SystemExit(f"plant_atlas_taxa: no name for {s['ddbid']}")
        gb = s["statusGB"]
        counts[gb] = counts.get(gb, 0) + 1
        # 2. drop statusGB = absent: 27 taxa that occur in Ireland or the Channel
        #    Islands but not in Great Britain, so they cannot be at the site
        if gb == "absent":
            continue
        # 3. everything else is kept, hybrids and 'native or alien' included, each
        #    carrying its own status so the coverage table can be cut by it. The
        #    hybrid binomials mostly have no NCBI node; they are rejected at
        #    resolution, where the reason is recorded, not dropped silently here
        out.append({"taxon": n["canonical"], "rank": n["rank"], "status_gb": gb,
                    "status_overall": s["overallStatus"], "ddbid": s["ddbid"]})
    return out, counts


#
# Reads the zip, derives the table and writes it with the attribution and the
# modification note CC BY 4.0 asks for
def main():
    ap = argparse.ArgumentParser(description="Writes the BSBI Plant Atlas 2020 Great Britain taxon list")
    ap.add_argument("--zip", help="PA2020_NativeStatusData.zip; downloaded from Zenodo if not given")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    blob, md5 = load(a.zip)
    rows, counts = derive(blob)
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write(f"# {CITE}\n")
        fh.write("# Licence: CC BY 4.0, https://creativecommons.org/licenses/by/4.0/\n")
        fh.write(f"# Source file: PA2020_NativeStatusData.zip, md5 {md5}, {ROWS} rows\n")
        fh.write("# Modified: status.csv joined to nameLookup.csv on ddbid; statusGB = 'absent' dropped;\n")
        fh.write("#           columns reduced to canonical name, rank, statusGB, overallStatus, ddbid.\n")
        fh.write("# Names are as the source writes them, unnormalised. See bin/plant_atlas_taxa.py\n")
        fh.write("taxon\trank\tstatus_gb\tstatus_overall\tddbid\n")
        for r in sorted(rows, key=lambda r: (r["taxon"], r["ddbid"])):
            fh.write("\t".join(r[k] for k in ("taxon", "rank", "status_gb", "status_overall", "ddbid")) + "\n")
    print(f"{len(rows)} taxa written to {a.out} (of {ROWS} in the source)", file=sys.stderr)
    for k, v in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"  {v:5d}  {k}{'   [dropped]' if k == 'absent' else ''}", file=sys.stderr)


if __name__ == "__main__":
    main()
