#!/usr/bin/env python3
#
# Draws the three data figures in docs/figures from a finished reanalysis run.
#   make_figures.py --results <outdir> --published Supporting_Data.xlsx --out docs/figures
#
import argparse, collections, csv, os, re
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

MONTHS = {4: "Apr", 5: "May", 6: "Jun", 7: "Jul", 8: "Aug", 9: "Sep", 10: "Oct"}
COL = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860", "#DA8BC3", "#6ACC64", "#D9D9D9"]
BEE, FLY = "#C9A227", "#3B7EA1"
SEASONS = {"Spring": "#5B8C3E", "Summer": "#D9A520", "Autumn": "#9C4A2F"}
# wind-pollinated: they should come out as avoided, which is the check that the index works
WIND = {"Poaceae", "Urticaceae", "Fagaceae", "Betulaceae", "Cyperaceae", "Juncaceae"}
norm = lambda s: re.sub(r"_0*(\d+)$", r"_\1", s.strip())

#
# genus or any rank name -> family, from the reference taxonomy of both markers
def family_map(results):
    fam = {}
    for m in ("rbcL", "ITS2"):
        for r in csv.DictReader(open(f"{results}/06_reference/refdb/{m}.taxonomy.tsv"), delimiter="\t"):
            ranks = (r.get("lineage") or "").split(";")
            f = next((x for x in ranks if x.endswith("aceae")), None)
            if f:
                fam.setdefault(r["genus"], f)
                for x in ranks[ranks.index(f):]:   # only the family and what sits under it: an
                    fam.setdefault(x, f)           # order holds many families and belongs to none
                                                   # of them, so charging its reads to one of its
                                                   # families would invent composition. Same rule
                                                   # as compare_published.py.
    return fam

#
# a taxon label ("Ranunculus spp.", "Heracleum sphondylium", "Apiaceae") -> family, or None
def to_family(label, fam):
    if label.endswith("aceae"):
        return label
    base = label.replace(" spp.", "").strip()
    for g in (base.split("/") if "/" in base else [base.split()[0]]):
        if g.strip() in fam:
            return fam[g.strip()]
    return fam.get(base.split()[0])

#
# per insect: the share of its load on each family, plus its date, season and class
def loads(results, meta_path, fam):
    meta = {norm(r["sample"]): r for r in csv.DictReader(open(meta_path), delimiter="\t")}
    rows = list(csv.reader(open(f"{results}/04_consensus/arm2/consensus_taxon_table.tsv"), delimiter="\t"))
    samples = [norm(s) for s in rows[0][1:]]
    counts = collections.defaultdict(dict)
    for r in rows[1:]:
        for s, v in zip(samples, r[1:]):
            if v and v != "0":
                counts[s][r[0]] = float(v)
    out = []
    for s, c in counts.items():
        m = meta.get(s)
        if not m or not c or m["insect_group"] == "NA":
            continue
        total = sum(c.values())
        share = collections.Counter()
        for k, v in c.items():
            f = to_family(k, fam)
            if f:
                share[f] += v / total
        # their "Month" is the sampling round, not the month of capture: half the rows disagree
        out.append(dict(sample=s, month=int(m["date"][5:7]), season=m["season"],
                        cls="Hoverfly" if m["insect_group"] == "Hoverfly" else "Bee", share=share))
    return out

#
# figure 1: stacked composition month by month, bees in one panel and hoverflies in the other
def composition(ins, top, out):
    panels = ["Bee", "Hoverfly"]
    cell, n = collections.defaultdict(collections.Counter), collections.Counter()
    for i in ins:
        n[(i["month"], i["cls"])] += 1
        for f, v in i["share"].items():
            cell[(i["month"], i["cls"])][f] += v
    nc = len(panels)
    fig, axes = plt.subplots(1, nc, figsize=(5.0 * nc, 5.0), sharey=True)
    axes = np.atleast_1d(axes).ravel()
    for ax, cl in zip(axes, panels):
        months = [b for b in sorted(MONTHS) if n[(b, cl)] >= 5]   # never draw a bar from a handful
        bottom = [0.0] * len(months)
        for fi, f in enumerate(top + ["other"]):
            vals = []
            for m in months:
                c, k = cell[(m, cl)], n[(m, cl)]
                vals.append((c[f] if f in top else sum(c.values()) - sum(c[t] for t in top)) / k)
            ax.bar(range(len(months)), vals, .78, bottom=bottom, color=COL[fi],
                   label=f if ax is axes[0] else None)
            bottom = [b + v for b, v in zip(bottom, vals)]
        ax.set_xticks(range(len(months)))
        ax.set_xticklabels(["%s\nn=%d" % (MONTHS[m], n[(m, cl)]) for m in months], fontsize=13.0)
        ax.set_title("%s (n=%d)" % (cl, sum(n[(m, cl)] for m in months)), fontsize=15.9)
        ax.set_ylim(0, 1)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_ylabel("share of an insect's pollen load")
    fig.legend(loc="lower center", ncol=5, frameon=False, fontsize=13.0)
    fig.suptitle("Which plant families the pollen came from, month by month", fontsize=17.4)
    fig.tight_layout(rect=[0, .13, 1, .93])
    save(fig, out, "composition_by_month.png")

#
# Jacobs' selectivity index: +1 is used and not available, -1 the reverse, 0 is use in proportion
def jacobs(used, available):
    d = used + available - 2 * used * available
    return (used - available) / d if d else 0.0

#
# figure 2: used against available, from the authors' own floral surveys (ten most abundant plants
# per transect and season, as a proportion of ground area in flower)
def preference(ins, fam, workbook, out):
    import openpyxl
    ws = openpyxl.load_workbook(workbook, read_only=True)["floral_surveys"]
    cover, tot = collections.defaultdict(float), collections.defaultdict(float)
    for r in [x for x in ws.iter_rows(values_only=True)][2:]:
        f = to_family(r[0], fam) if r[0] else None
        if f:
            cover[(r[2], r[3], f)] += float(r[4])
            tot[(r[2], r[3])] += float(r[4])
    # the sheet records area in flower, not shares: a transect's rows sum to anything from 0.03 to
    # 3.4, so take each transect as a share of its own cover and average the transects of a season
    avail = collections.defaultdict(float)
    nt = collections.Counter(se for se, tr in tot)
    for (se, tr, f), v in cover.items():
        avail[(se, f)] += v / tot[(se, tr)] / nt[se]
    use, n = collections.defaultdict(float), collections.Counter()
    for i in ins:
        n[(i["season"], i["cls"])] += 1
        for f, v in i["share"].items():
            use[(i["season"], i["cls"], f)] += v
    for key in list(n):
        for k in [k for k in use if (k[0], k[1]) == key]:
            use[k] /= n[key]

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 5.8), sharex=True)
    for ax, se in zip(axes, ("Spring", "Summer")):
        fams = sorted({f for (s2, f) in avail if s2 == se and avail[(s2, f)] >= .01})
        rows = [(f, jacobs(use.get((se, "Bee", f), 0.0), avail[(se, f)]),
                 jacobs(use.get((se, "Hoverfly", f), 0.0), avail[(se, f)]), avail[(se, f)]) for f in fams]
        rows.sort(key=lambda r: (r[1] + r[2]) / 2)
        y = np.arange(len(rows))
        ax.barh(y + .19, [r[1] for r in rows], .38, color=BEE, label="bees" if se == "Spring" else None)
        ax.barh(y - .19, [r[2] for r in rows], .38, color=FLY, label="hoverflies" if se == "Spring" else None)
        ax.set_yticks(y)
        ax.set_yticklabels(["%s%s  (%.0f%%)" % ("~ " if f in WIND else "", f, 100 * p)
                            for f, _, _, p in rows], fontsize=12.3)
        ax.axvline(0, color="#333", lw=.9)
        ax.set_xlim(-1.05, 1.05)
        ax.set_xlabel("avoided        Jacobs' index        preferred")
        ax.set_title("%s\n%d bees, %d hoverflies" % (se, n[(se, "Bee")], n[(se, "Hoverfly")]), fontsize=15.0)
        ax.grid(axis="x", color="#EEEEEE"); ax.set_axisbelow(True)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.tick_params(left=False)
    fig.legend(loc="lower center", ncol=2, frameon=False, fontsize=13.0)
    fig.suptitle("Do they use the flowers that are there? Pollen load against flowering area", fontsize=18.1)
    fig.tight_layout(rect=[0, .04, 1, .95])
    save(fig, out, "preference.png")

#
# figure 3: the same insects ordinated once per version of the table, as the paper's own Figure 5.
# The coordinates come from bin/ordination.R, not from phase 6, because phase 6 keeps every sample
# and two of them cannot be placed; see docs/figures.md. The panels answer only whether the three
# tables arrange the insects the same way. How closely they do is the Procrustes r in the title,
# and the agreement that does not depend on a two-dimensional fit at all is in docs/results.md.
def ordination(ordir, meta_path, out):
    from matplotlib.patches import Ellipse
    arms = ["published", "replica", "reanalysis"]
    meta = {norm(r["sample"]): r for r in csv.DictReader(open(meta_path), delimiter="\t")}
    fit = collections.defaultdict(dict)
    for r in csv.DictReader(open(f"{ordir}/ordination_fit.tsv"), delimiter="\t"):
        fit[r["table"]][r["metric"]] = r["value"]

    def coords(k):
        rows = list(csv.DictReader(open(f"{ordir}/{k}.nmds.tsv"), delimiter="\t"))
        return {norm(r["sample"]): (float(r["NMDS1"]), float(r["NMDS2"])) for r in rows}

    # centre, scale to unit sum of squares, rotate Y onto X: the rotation vegan::procrustes finds.
    # NMDS axes carry no orientation, so without this the panels cannot be compared by eye.
    def rotate(X, Y):
        X, Y = X - X.mean(0), Y - Y.mean(0)
        X, Y = X / np.sqrt((X ** 2).sum()), Y / np.sqrt((Y ** 2).sum())
        u, _, vt = np.linalg.svd(X.T @ Y)
        return Y @ (u @ vt).T

    tab = {k: coords(k) for k in arms}
    keep = [s for s in sorted(set.intersection(*(set(t) for t in tab.values()))) if s in meta]
    pts = {k: np.array([tab[k][s] for s in keep]) for k in arms}
    ref = rotate(pts["published"], pts["published"])
    drawn = {k: (ref if k == "published" else rotate(ref, pts[k])) for k in arms}
    lim = max(np.percentile(np.abs(P), 99.5) for P in drawn.values()) * 1.1

    fig, axes = plt.subplots(1, 3, figsize=(10.8, 5.6))
    for ax, k in zip(axes, arms):
        P = drawn[k]
        for cls, mk in (("Bee", "o"), ("Hoverfly", "^")):
            for se in SEASONS:
                i = [j for j, s in enumerate(keep)
                     if meta[s]["season"] == se
                     and ("Hoverfly" if meta[s]["insect_group"] == "Hoverfly" else "Bee") == cls]
                if i:
                    ax.scatter(P[i, 0], P[i, 1], s=16, marker=mk, c=SEASONS[se], alpha=.75,
                               linewidths=0, label=f"{se}, {cls.lower()}" if ax is axes[0] else None)
        for se in SEASONS:
            i = [j for j, s in enumerate(keep) if meta[s]["season"] == se]
            if len(i) >= 20:            # autumn is 11 insects: an ellipse on that is noise
                w, v = np.linalg.eigh(np.cov(P[i].T))
                ax.add_patch(Ellipse(P[i].mean(0), 2 * np.sqrt(w[-1]), 2 * np.sqrt(w[0]),
                                     angle=np.degrees(np.arctan2(v[1, -1], v[0, -1])),
                                     fill=False, edgecolor=SEASONS[se], lw=1.5, alpha=.9))
        t = "%s\nstress %.3f" % (k, float(fit[k]["stress"]))
        r = fit.get(f"published vs {k}", {}).get("procrustes_r")
        ax.set_title(t + (f", r = {float(r):.3f}" if r else ""), fontsize=15.0)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        for s in ("top", "right", "bottom", "left"):
            ax.spines[s].set_visible(False)
    fig.legend(loc="lower center", ncol=3, frameon=False, fontsize=13.0)
    fig.suptitle("The same insects, ordinated from each version of the table", fontsize=18.0)
    fig.tight_layout(rect=[0, .13, 1, .93])
    save(fig, out, "ordination.png")

#
# figure 4: the paper's Figure 3 question -- do the pollinator groups carry different plants? --
# asked of all three tables at once. Composition is charged to plant families rather than to taxon
# labels, because the three tables do not share a vocabulary: 191 published labels against 784 and
# 422. Family is the coarsest thing all three agree on and so the only fair common currency.
def composition_by_table(tables, meta_path, fam, out):
    meta = {norm(r["sample"]): r for r in csv.DictReader(open(meta_path), delimiter="\t")}
    groups = ["Bumblebee", "Honeybee", "Solitary bee", "Hoverfly"]
    names = list(tables)

    #  table -> group -> family -> summed share of a load, and the insects behind it
    share, n = {}, {}
    for k, path in tables.items():
        rows = [l.rstrip("\n").split("\t") for l in open(path)]
        samples = [norm(s) for s in rows[0][1:]]
        cell = collections.defaultdict(collections.Counter)
        seen = collections.Counter()       # counted per table, not once: a table that drops an
        per = collections.defaultdict(collections.Counter)   # insect must divide by its own n
        whole = collections.Counter()      # the whole load, including labels with no family,
        for r in rows[1:]:                 # so this denominator matches the one loads() uses
            f = to_family(r[0], fam)
            for s, v in zip(samples, r[1:]):
                if v and v != "0":
                    whole[s] += float(v)
                    if f:
                        per[s][f] += float(v)
        for s, c in per.items():
            m = meta.get(s)
            if not m or m["insect_group"] not in groups or not c:
                continue
            for f, v in c.items():
                cell[m["insect_group"]][f] += v / whole[s]
                if m["insect_group"] != "Hoverfly":
                    cell["bees"][f] += v / whole[s]     # the paper's own split is bees against
            seen[m["insect_group"]] += 1                # hoverflies, and the README quotes it
            seen["bees"] += m["insect_group"] != "Hoverfly"
        share[k], n[k] = cell, seen

    # the families to name: the biggest across every table and group, so one legend serves all
    pooled = collections.Counter()
    for cell in share.values():
        for g in groups:
            for f, v in cell[g].items():
                pooled[f] += v
    top = [f for f, _ in pooled.most_common(8)]

    # the same numbers as a table, because the README quotes four of them and a figure is a poor
    # place to read a third decimal off
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "family_share_by_group.tsv"), "w") as fh:
        fh.write("family\tinsect_group\t" + "\t".join(names) + "\n")
        rows = groups + ["bees"]
        for g in rows:                     # the denominators, so a divergence is visible
            fh.write("(insects)\t" + g + "\t" + "\t".join(str(n[k][g]) for k in names) + "\n")
        for f in top:
            for g in rows:
                fh.write(f + "\t" + g + "\t" +
                         "\t".join(f"{share[k][g][f] / (n[k][g] or 1):.4f}" for k in names) + "\n")

    fig, axes = plt.subplots(1, 4, figsize=(11.4, 5.4), sharey=True)
    for ax, g in zip(axes, groups):
        for bi, k in enumerate(names):
            c = share[k][g]
            tot = sum(c.values()) or 1.0
            bottom = 0.0
            for fi, f in enumerate(top + ["other"]):
                v = (c[f] if f in top else tot - sum(c[t] for t in top)) / tot
                ax.bar(bi, v, .72, bottom=bottom, color=COL[fi],
                       label=f if (ax is axes[0] and bi == 0) else None)
                bottom += v
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names, fontsize=13.0, rotation=20, ha="right")
        ax.set_title("%s\nn = %d" % ("Non-corbiculate bee" if g == "Solitary bee" else g,
                                      n[names[0]][g]), fontsize=15.9)
        ax.set_ylim(0, 1)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_ylabel("mean share of a pollen load", fontsize=13.0)
    fig.legend(loc="lower center", ncol=5, frameon=False, fontsize=13.0)
    fig.suptitle("Which plant families each pollinator group carried, under each version of the table", fontsize=18.1)
    fig.tight_layout(rect=[0, .13, 1, .93])
    save(fig, out, "composition_by_table.png")

#
# matplotlib cannot always write straight to a Windows mount, so go via /tmp
def save(fig, out, name):
    import shutil
    tmp = f"/tmp/{name}"
    fig.savefig(tmp, dpi=200)
    os.makedirs(out, exist_ok=True)
    shutil.copy(tmp, os.path.join(out, name))
    plt.close(fig)
    print(os.path.join(out, name))


def main():
    ap = argparse.ArgumentParser(description="Draw the data figures from a reanalysis run")
    ap.add_argument("--results", required=True, help="the reanalysis outdir")
    ap.add_argument("--metadata", required=True, help="insect metadata from ecology_metadata.py")
    ap.add_argument("--published", required=True, help="Supporting_Data.xlsx, for the floral surveys")
    ap.add_argument("--ordination", help="output dir of bin/ordination.R; the ordination figure is skipped without it")
    ap.add_argument("--published-table", help="the published sheet as a taxon x sample TSV, from bin/published_table.py")
    ap.add_argument("--out", default="docs/figures")
    a = ap.parse_args()

    fam = family_map(a.results)
    ins = loads(a.results, a.metadata, fam)
    share = collections.Counter()
    for i in ins:
        for f, v in i["share"].items():
            share[f] += v / len(ins)
    top = [f for f, _ in share.most_common(8)]
    # the paper's own split, bees against hoverflies, and the finer one the metadata allows
    composition(ins, top, a.out)
    preference(ins, fam, a.published, a.out)
    if a.ordination:
        ordination(a.ordination, a.metadata, a.out)
    if a.published_table:
        composition_by_table({"published": a.published_table,
                              "replica": f"{a.results}/04_consensus/arm1/consensus_taxon_table.tsv",
                              "reanalysis": f"{a.results}/04_consensus/arm2/consensus_taxon_table.tsv"},
                             a.metadata, fam, a.out)


if __name__ == "__main__":
    main()
