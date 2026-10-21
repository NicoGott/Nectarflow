#!/usr/bin/env python3
#
# One parse of every source, cached as CSV/JSON in explore/data/, so figure.py reads tidy
# files instead of re-opening the workbook. Run once: python3 explore/prep.py
#
import csv, json, collections, os, statistics, sys, openpyxl

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUN  = os.environ.get("REPLICA_RUN") or sys.exit("set REPLICA_RUN to the results directory")
F, A = f"{RUN}/faithful", f"{RUN}/reanalysis"
OUT  = f"{REPO}/explore/data"
os.makedirs(OUT, exist_ok=True)

SPEC = {"species": 4, "genus": 3, "genus_group": 2, "family": 1, "above": 0}
NAME = {4: "species", 3: "genus", 2: "genus group", 1: "family", 0: "above family"}

# genus or family name -> family, from the reference taxonomy both arms searched
fam = {}
for m in ("rbcL", "ITS2"):
    for r in csv.DictReader(open(f"{A}/06_reference/refdb/{m}.taxonomy.tsv"), delimiter="\t"):
        if r["family"]:
            if r["genus"]:
                fam.setdefault(r["genus"], r["family"])
            fam.setdefault(r["family"], r["family"])


# Tribes sit between genus and family. Their tables use four of them, and without
# this map they parse as binomials and land in the species bucket, which overstates
# how finely the published analysis resolved. A tribe is the same statement as their
# "A/B/C spp." -- declining to choose among the genera of a clade -- so it is counted
# as a genus group, under the family the tribe belongs to.
TRIBE = {"Anthemideae": "Asteraceae", "Cichorieae": "Asteraceae",
         "Mentheae": "Lamiaceae", "Genisteae": "Fabaceae"}


# a label ("Rubus spp.", "Heracleum sphondylium", "Asteraceae") -> (rank, family or None)
def parse(label):
    n = label.strip()
    if n.endswith("aceae"):
        return "family", n
    if n in TRIBE:
        return "genus_group", TRIBE[n]
    if n.endswith(("ales", "idae", "phyta", "opsida")):
        return "above", None
    if n.endswith("spp."):
        b = n[:-4].strip()
        gs = [g.strip() for g in b.split("/")] if "/" in b else [b]
        return ("genus_group" if len(gs) > 1 else "genus"), next((fam[g] for g in gs if g in fam), None)
    return "species", fam.get(n.split()[0])


def wide(path):                       # taxon x sample TSV -> {sample: [(label, reads)]}
    rows = list(csv.reader(open(path), delimiter="\t"))
    head, out = rows[0][1:], collections.defaultdict(list)
    for r in rows[1:]:
        for s, v in zip(head, r[1:]):
            if v and v != "0":
                out[s].append((r[0], int(v)))
    return out


wb = openpyxl.load_workbook(f"{REPO}/assets/published/Supporting_Data.xlsx", read_only=True)


def sheet(name, skip=8):              # sample x taxon workbook sheet -> {sample: [(label, reads)]}
    rows = [x for x in wb[name].iter_rows(values_only=True)]
    labels, out = [str(c) for c in rows[1][skip:]], collections.defaultdict(list)
    for r in rows[2:]:
        if not r or not r[0]:
            continue
        for nm, v in zip(labels, r[skip:]):
            if isinstance(v, (int, float)) and v > 0:
                out[str(r[0])].append((nm, int(v)))
    return out


ours  = wide(f"{A}/04_consensus/arm2/consensus_taxon_table.tsv")
arm1  = wide(f"{F}/08_tables/rbcL.taxon_table.tsv")
arm1b = wide(f"{F}/08_tables/ITS2.taxon_table.tsv")
theirs = sheet("rbcL_ITS2_consensus_read_counts")

# ---- 1. per (insect, family): the finest rank each side reached, and the reads behind it
def finest(dets):
    b = {}
    for label, reads in dets:
        rk, f = parse(label)
        if not f:
            continue
        cur = b.get(f, (-1, 0, ""))
        b[f] = (max(cur[0], SPEC[rk]), cur[1] + reads, label if SPEC[rk] > cur[0] else cur[2])
    return b


both = sorted(set(ours) & set(theirs))
with open(f"{OUT}/cells.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["insect", "family", "their_rank", "our_rank", "their_reads", "our_reads",
                "their_label", "our_label", "seen_by"])
    for s in both:
        a, t = finest(ours[s]), finest(theirs[s])
        for f in sorted(set(a) | set(t)):
            ta = t.get(f); oa = a.get(f)
            w.writerow([s, f,
                        ta[0] if ta else "", oa[0] if oa else "",
                        ta[1] if ta else 0, oa[1] if oa else 0,
                        ta[2] if ta else "", oa[2] if oa else "",
                        "both" if ta and oa else ("theirs" if ta else "ours")])

# ---- 2. read flow through the stages, all three
tot = [r for r in csv.DictReader(open(f"{F}/read_tracking.tsv"), delimiter="\t")
       if r["sample"] == "TOTAL"][0]
d2 = {}
for m in ("rbcL", "ITS2"):
    rows = list(csv.DictReader(open(f"{A}/02_asv/{m}.dada2_track.tsv"), delimiter="\t"))
    for c in ("cutadapt_pairs", "filtered_pairs", "denoised", "merged", "nonchimeric"):
        d2[f"{m}_{c}"] = sum(int(float(r[c])) for r in rows if r[c])

rank2 = collections.defaultdict(dict)
for m in ("rbcL", "ITS2"):
    for r in csv.DictReader(open(f"{A}/03_assignment/{m}.rank_summary.tsv"), delimiter="\t"):
        if r["level"] == "rank":
            rank2[m][r["value"]] = int(r["reads"])

flow = dict(
    published=dict(raw=40800709, named=22510682, rbcL_named=11305697, ITS2_named=11204985),
    arm1={k: int(tot[k]) for k in tot if k != "sample"},
    arm2=dict(d2, **{f"{m}_rank_{k}": v for m, d in rank2.items() for k, v in d.items()}),
)
json.dump(flow, open(f"{OUT}/flow.json", "w"), indent=1)

# ---- 3. rank marginals per marker and side, by reads and by label
def marg(dets_by_sample):
    rr, ll = collections.Counter(), collections.Counter()
    seen = collections.defaultdict(set)
    for s, dets in dets_by_sample.items():
        for label, reads in dets:
            rk, _ = parse(label)
            rr[rk] += reads; seen[rk].add(label)
    for k, v in seen.items():
        ll[k] = len(v)
    return rr, ll


rows = []
for side, data, marker in (("published", sheet("rbcL_consensus_ID_read_counts"), "rbcL"),
                           ("published", sheet("ITS2_consensus_ID_read_counts"), "ITS2"),
                           ("arm1", arm1, "rbcL"), ("arm1", arm1b, "ITS2"),
                           ("arm2", wide(f"{A}/03_assignment/rbcL.taxon_table.tsv"), "rbcL"),
                           ("arm2", wide(f"{A}/03_assignment/ITS2.taxon_table.tsv"), "ITS2")):
    rr, ll = marg(data)
    for rk in ("species", "genus", "genus_group", "family", "above"):
        rows.append(dict(side=side, marker=marker, rank=rk, reads=rr[rk], labels=ll[rk]))
with open(f"{OUT}/marginals.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, ["side", "marker", "rank", "reads", "labels"]); w.writeheader(); w.writerows(rows)

print(f"cells.csv     {sum(1 for _ in open(f'{OUT}/cells.csv')) - 1} rows")
print(f"marginals.csv {len(rows)} rows")
print(f"flow.json     arm1 named {int(tot['rbcL_assigned']) + int(tot['ITS2_assigned']):,}")
