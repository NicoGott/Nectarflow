#!/usr/bin/env python3
#
# Measures what the reference library can and cannot name. This is a report,
# never an input: it reads the authors' published results to ask how much of
# what they found the library would now miss, and it runs after the manifest is
# written and hashed, so nothing it says can have steered the build. The
# manifest md5 it was run against is printed at the top for that reason.
#
# Reports, per marker: how many target taxa carry a sequence, the library size,
# coverage split by Plant Atlas native status, and - as a measurement after the
# fact - how many of the taxa Lowe et al. 2022 report would be absent.
#
import argparse, collections, csv, hashlib, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from refnames import binomial, genus, normalise

MARKERS = ("rbcL", "ITS2")


#
# Reads a TSV, skipping the comment header
def read_tsv(path):
    return list(csv.DictReader((l for l in open(path, encoding="utf-8") if not l.startswith("#")), delimiter="\t"))


#
# The names a library can answer with: every species, genus, family and order
# label it holds, per marker
def labels(rows):
    out = {}
    for m in sorted({r["marker"] for r in rows}):
        mine = [r for r in rows if r["marker"] == m]
        rank = set()
        for r in mine:
            rank |= {x for x in r["lineage"].split(";") if x.endswith(("aceae", "eae", "ales"))}
        out[m] = {"species": {r["species"] for r in mine if r["species"]} |
                             {binomial(r["target_name"]) for r in mine if binomial(r["target_name"])},
                  "genus": {r["genus"] for r in mine if r["genus"]} |
                           {genus(r["target_name"]) for r in mine if genus(r["target_name"])},
                  "family": {r["family"] for r in mine if r["family"]},
                  "order": {r["order"] for r in mine if r["order"]},
                  "rank_name": rank}
    return out


#
# Can the library reach one published label? Their vocabulary is 'Genus
# species', 'Genus spp.', 'A/B spp.' and bare family or tribe names, so each is
# asked of the rank it is written at
def reach(name, lab):
    n = " ".join(str(name).replace("\xa0", " ").split())
    head = n.split()[0] if n.split() else ""
    # family, tribe and order labels are answered from the frozen lineages, which carry all three
    if head.endswith(("aceae", "eae", "ales")):
        return head in lab["rank_name"] or head in lab["family"] or head in lab["order"]
    if " spp." in n or "/" in n:
        return any(g.strip() in lab["genus"] for g in n.replace(" spp.", "").split("/") if g.strip())
    b = binomial(n)
    if b and b in lab["species"]:
        return True
    return genus(n) in lab["genus"]


#
# One row of the coverage table: what Barcode UK alone reaches and what the
# whole library reaches, at species rank and at genus rank or better
def line(key, want, bu_sp, bu_ge, lab):
    n = len(want)
    f = lambda k: f"{k:5d} ({100 * k / n:3.0f}%)"
    return (f"{key:16s} {n:5d} targets | barcode_uk species {f(sum(1 for b in want if b in bu_sp))} "
            f"genus+ {f(sum(1 for b in want if b in bu_sp or b.split()[0] in bu_ge))} "
            f"| library species {f(sum(1 for b in want if b in lab['species']))} "
            f"genus+ {f(sum(1 for b in want if b in lab['species'] or b.split()[0] in lab['genus']))}")


#
# Prints every figure the honesty statement rests on
def main():
    ap = argparse.ArgumentParser(description="Measures what the reference library can and cannot name")
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--plant-atlas", required=True)
    ap.add_argument("--rejects")
    ap.add_argument("--published", help="Supporting_Data.xlsx; read only to measure the gap")
    ap.add_argument("--shipped", help="the old refdb directory, for the size comparison")
    a = ap.parse_args()
    rows = read_tsv(a.manifest)
    md5 = hashlib.md5(open(a.manifest, "rb").read()).hexdigest()
    print(f"manifest {a.manifest}  md5 {md5}  {len(rows)} rows")
    if "PARTIAL RESOLVE" in open(a.manifest, encoding="utf-8").readline():
        print("WARNING: this is a partial resolve. The GenBank fill columns below are a sample,\n"
              "         not the library. Barcode UK is complete.")
    lab = labels(rows)

    print("\n== library ==")
    for m in sorted(lab):
        mine = [r for r in rows if r["marker"] == m]
        bu = sum(1 for r in mine if r["source"] == "barcode_uk")
        L = sorted(int(r["length"]) for r in mine)
        print(f"{m}: {len(mine)} sequences ({bu} barcode_uk, {len(mine) - bu} genbank); "
              f"{len(lab[m]['species'])} species, {len(lab[m]['genus'])} genera, {len(lab[m]['family'])} families; "
              f"length min {L[0]} med {L[len(L) // 2]} max {L[-1]}")
        if a.shipped:
            p = os.path.join(a.shipped, f"{m}.fasta")
            if os.path.exists(p):
                n = sum(1 for l in open(p) if l.startswith(">"))
                print(f"     shipped library: {n} sequences  ({len(mine) - n:+d})")

    print("\n== targets reached ==  (species = the target's own binomial is in the library;"
          " genus+ = it or a congener is)")
    atlas = read_tsv(a.plant_atlas)
    status = {}
    for r in atlas:
        b = binomial(r["taxon"])
        if b:
            status[b] = r["status_gb"]
    for m in sorted(lab):
        bu_sp = {r["species"] for r in rows if r["marker"] == m and r["source"] == "barcode_uk" and r["species"]}
        bu_ge = {r["genus"] for r in rows if r["marker"] == m and r["source"] == "barcode_uk" and r["genus"]}
        print(f"-- {m}")
        for key in ("native", "archaeophyte", "neophyte", "native or alien"):
            want = sorted(b for b, s in status.items() if s == key)
            if not want:
                continue
            print("   " + line(key, want, bu_sp, bu_ge, lab[m]))
    if a.rejects:
        rj = read_tsv(a.rejects)
        print(f"\n   rejected target names: {len(rj)}  " +
              "  ".join(f"{k}={v}" for k, v in collections.Counter(r["reason"].split(":")[0] for r in rj).most_common()))

    if a.published:
        import openpyxl
        wb = openpyxl.load_workbook(a.published, read_only=True)
        print("\n== measured after the fact, against the authors' results ==")
        cons = {" ".join(str(r[0]).replace("\xa0", " ").split())
                for r in list(wb["consensus_identifications"].iter_rows(values_only=True))[2:] if r[0]}
        flor = {" ".join(str(r[0]).replace("\xa0", " ").split())
                for r in list(wb["floral_surveys"].iter_rows(values_only=True))[2:] if r[0]}
        union = {k: set().union(*(lab[m][k] for m in lab)) for k in
                 ("species", "genus", "family", "order", "rank_name")}
        for what, names in (("published consensus names", sorted(cons)),
                            ("floral survey taxa (held out)", sorted(flor))):
            miss = [n for n in names if not reach(n, union)]
            print(f"   {what:32s} {len(names) - len(miss):4d} / {len(names):4d} reachable; "
                  f"missing: {', '.join(miss[:12]) if miss else 'none'}")


if __name__ == "__main__":
    main()
