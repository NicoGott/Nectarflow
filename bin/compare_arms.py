#!/usr/bin/env python3
#
# Compares the two arms' taxon tables with each other and with the published ones (phase 5).
# Arm 1 against arm 2 is computed here; the two published columns are read straight out of the
# per-sample files compare_published.py wrote.
#
import argparse, collections, csv, os, re, statistics, sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import compare_published as cp      # norm / genera / spearman; its openpyxl import is inside main()

# The accounting rows of arm 2's <marker>.taxon_table.full.tsv: reads, not taxa, so they never
# enter a composition.
BUCKET = "__"
PER_SAMPLE = ("sample in_published reads_arm1 reads_arm2 reads_published "
              "unclassified_arm1 unclassified_arm2 taxa_arm1 taxa_arm2 "
              "bc_arm1_vs_arm2 bc_arm1_vs_published bc_arm2_vs_published spearman_arm1_vs_arm2 "
              "top_arm1 top_arm2 top_published "
              "top_arm1_agrees_arm2 top_arm1_agrees_published top_arm2_agrees_published").split()
# n_scored is beside n_samples because the medians and counts are over the scored samples only,
# while n_samples counts every sample compared; one denominator for both would be wrong.
SUMMARY = ("comparison marker resolution n_samples n_scored median_bray_curtis_similarity mean_bray_curtis_similarity "
           "samples_ge_0.9 samples_lt_0.5 dominant_taxon_identical median_read_count_diff_pct").split()


#
# Reads a taxon x sample table. Returns the samples, the table, how many accounting rows were
# skipped, the reads on them per sample, and how many duplicate taxon rows were summed
def read_table(path):
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        sys.exit(f"{path} is empty")
    header = rows[0]
    samples = [cp.norm(s) for s in header[1:]]
    if len(set(samples)) != len(samples):
        # Two columns normalising to one name would be summed into one sample unseen; the
        # duplicate second run of GTF3_919_2 keeps a suffixed name precisely to avoid this.
        dup = sorted({s for s in samples if samples.count(s) > 1})
        sys.exit(f"{path}: sample columns collide after normalisation: {', '.join(dup)}")
    table, buckets, dropped, merged = {}, 0, collections.Counter(), 0
    for lineno, r in enumerate(rows[1:], start=2):
        if not r or not r[0]:
            continue
        if len(r) != len(header):
            sys.exit(f"{path} line {lineno} ({r[0]!r}) has {len(r)} fields, the header has "
                     f"{len(header)}; a short row would silently read as zero reads")
        if r[0].startswith(BUCKET):
            buckets += 1
            for s, v in zip(samples, r[1:]):
                dropped[s] += float(v or 0)
            continue
        if r[0] in table:                       # summed, as compare_published.py sums them
            merged += 1
        c = table.setdefault(r[0], collections.Counter())
        for s, v in zip(samples, r[1:]):
            c[s] += float(v or 0)
    return sorted(samples), table, buckets, dropped, merged


#
# A {name: family} map, from three sources in increasing authority:
# the first '-aceae' name in a lineage, but only when every lineage carrying that element
#   agrees, so a tribe resolves and an order correctly does not
# the genus and family columns
# an arm's own taxon_lineage.tsv, keyed by the label it wrote
def family_map(paths):
    lineage = collections.defaultdict(set)
    fam = {}
    for path in paths:
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("lineage"):
                    names = [x for x in r["lineage"].split(";") if x]
                    here = next((n for n in names if n.endswith("aceae")), None)
                    for n in names:
                        lineage[n].add(here)
                if r.get("family"):
                    if r.get("genus"):
                        fam[r["genus"]] = r["family"]
                    if r.get("taxon"):                  # a taxon_lineage.tsv: the arm's own statement
                        fam[r["taxon"]] = r["family"]
                    fam[r["family"]] = r["family"]
    for name, fs in lineage.items():
        if name not in fam and len(fs) == 1 and None not in fs:
            fam[name] = fs.pop()
    return fam


#
# The family of a label, or None when its genera disagree or none is known
def family_of(label, fam):
    if label in fam:
        return fam[label]
    fs = {fam.get(g) for g in cp.genera(label)}
    return fs.pop() if len(fs) == 1 and None not in fs else None


#
# Sums a table to family level. A label without a family keeps its own label and is listed
def roll_up(table, fam):
    out, unmapped = collections.defaultdict(collections.Counter), set()
    for taxon, counts in table.items():
        f = family_of(taxon, fam)
        if f is None:                               # keeps its own label, and is counted
            unmapped.add(taxon)
            f = taxon
        out[f].update(counts)
    return dict(out), sorted(unmapped)


#
# The rank a label claims, read off its shape. Valid only because arm 2 keeps arm 1's
# vocabulary: 'Genus epithet' / 'Genus spp.' / 'Genus sp.' / 'A/B spp.' / a family or a higher name
#
# Bray-Curtis SIMILARITY and Spearman rho between two samples. Similarity, not the dissimilarity
# vegan returns: sum(min(p_i, q_i)) over proportions, so 1 is the same composition and 0 shares
# nothing. The statistics in ecology.R use the dissimilarity, as adonis2 and metaMDS require.
# NA for both when either sample is empty: an empty blank is not a composition and must not score 0
def bray_curtis(p, q):
    tp, tq = sum(p.values()), sum(q.values())
    if tp <= 0 or tq <= 0:
        return None, None
    present = sorted({t for t in p if p[t] > 0} | {t for t in q if q[t] > 0})
    x = [p.get(t, 0) / tp for t in present]
    y = [q.get(t, 0) / tq for t in present]
    return sum(min(u, v) for u, v in zip(x, y)), cp.spearman(x, y)


#
# Which resolution a compare_published.py file was written at, taken from its name
def published_resolution(path):
    return "family" if path and ".family." in os.path.basename(path) else "published_label"


#
# Reads a compare_published.py per-sample file, keyed by normalised sample name
def read_vs_published(path):
    with open(path) as fh:
        lines = fh.readlines()
    rows = csv.DictReader([l for l in lines if not l.startswith("#")], delimiter="\t")
    return {cp.norm(r["sample"]): r for r in rows}


#
# Formats a number, with None and NaN both coming out as NA
def fmt(v, digits=3):
    return "NA" if v is None or v != v else f"{v:.{digits}f}"


#
# One row of the summary table. The medians and counts are over the scored samples only, which
# is why n_scored sits beside n_samples
def summary_row(comparison, marker, resolution, bcs, agree, n, diffs):
    ok = [b for b in bcs if b is not None]
    return {"comparison": comparison, "marker": marker, "resolution": resolution, "n_samples": n,
            "n_scored": len(ok),
            "median_bray_curtis_similarity": fmt(statistics.median(ok) if ok else None),
            "mean_bray_curtis_similarity": fmt(statistics.fmean(ok) if ok else None),
            "samples_ge_0.9": sum(1 for b in ok if b >= 0.9),
            "samples_lt_0.5": sum(1 for b in ok if b < 0.5),
            "dominant_taxon_identical": agree,
            "median_read_count_diff_pct": fmt(statistics.median(diffs) if diffs else None, 2)}


#
# Writes a TSV, with optional '#' comment lines before the header and after the last row
def write_tsv(path, columns, rows, head=None, tail=None):
    with open(path, "w") as fh:
        if head:
            fh.write("# " + head + "\n")
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in columns) + "\n")
        if tail:
            fh.write("# " + tail + "\n")


#
# Relabels a table through a name_translation.tsv, so an arm is read in the authors'
# vocabulary. Returns the new table and how many labels actually moved
def translate(table, path):
    if not path:
        return table, 0
    mapping = {r["replica_taxon"]: r["published_label"]
               for r in csv.DictReader(open(path), delimiter="\t") if r.get("published_label")}
    out, hit = collections.defaultdict(collections.Counter), 0
    for taxon, counts in table.items():
        label = mapping.get(taxon, taxon)
        hit += label != taxon
        for sample, v in counts.items():
            out[label][sample] += v
    return out, hit


#
# Does the whole comparison and writes the four tables: per sample, the summary, the medians
# per site group, and where the resolution went between the arms
def run(a):
    samples1, label1, b1, drop1, m1 = read_table(a.arm1)
    samples2, label2, b2, drop2, m2 = read_table(a.arm2)
    fam = family_map(a.lineage)
    # One stem for the three files written beside --out. `a.out.replace(".tsv", ...)` does
    # NOTHING when --out has no .tsv suffix: all four tables would land on --out, each
    # overwriting the last.
    stem = a.out[:-4] if a.out.endswith(".tsv") else a.out
    resolution = "family" if a.rank == "family" else "published_label"
    note = f"marker={a.marker} rank={a.rank} samples_arm1={len(samples1)} samples_arm2={len(samples2)}"
    if m1 or m2:
        note += f" duplicate_taxon_rows_summed={m1}+{m2}"
    t1, t2 = label1, label2
    if a.rank == "label" and (a.translate1 or a.translate2):
        # Both arms into the authors' vocabulary, so this row can be read beside the published ones.
        t1, n1 = translate(label1, a.translate1)
        t2, n2 = translate(label2, a.translate2)
        resolution = "published_label (translated)"
        note += f" translated_labels={n1}+{n2}"
    if a.rank == "family":
        if not fam:
            sys.exit("--rank family needs --lineage (a reference taxonomy or taxon_lineage TSV)")
        t1, u1 = roll_up(label1, fam)
        t2, u2 = roll_up(label2, fam)
        note += (f" labels_without_a_family={len(u1)}+{len(u2)}"
                 + (f" ({';'.join(sorted(set(u1) | set(u2))[:10])})" if u1 or u2 else ""))
    # NA here means "the table passed does not account for unclassified reads" (the plain
    # taxon_table.tsv, and arm 1's always), not "none were dropped".
    note += (f" accounting_rows={b1}+{b2}"
             f" unclassified_reads=arm1:{int(sum(drop1.values())) if b1 else 'NA'}"
             f" arm2:{int(sum(drop2.values())) if b2 else 'NA'}")

    pub1 = read_vs_published(a.vs_published_arm1) if a.vs_published_arm1 else {}
    pub2 = read_vs_published(a.vs_published_arm2) if a.vs_published_arm2 else {}
    # Each published column is named after the file it came from, never after this run's --rank:
    # a family-level median read as a label-level one overstates agreement by exactly what the
    # roll-up buys. compare_published.py does not record its --rank inside the file, so the tag
    # in the name is what says so.
    res1, res2 = (published_resolution(p) for p in (a.vs_published_arm1, a.vs_published_arm2))
    # The UNION, not the intersection: both arms drop a sample with no surviving reads, so
    # intersecting would delete exactly the samples one arm lost, which is a finding. Such a
    # sample scores NA here and shows as a zero read count.
    shared = sorted(set(samples1) | set(samples2))
    in_pub = set(pub1) | set(pub2)
    note += f" compared={len(shared)} in_both={len(set(samples1) & set(samples2))}"
    # reads_published and top_published are ONE column each, so they describe ONE of the two
    # files: the one written at this run's own resolution. Agreement with the published top
    # taxon is only meaningful there, since our roll-up changes our labels but not theirs.
    pub_cols, pub_from = (pub1, "arm1") if pub1 else (pub2, "arm2")
    if res2 == resolution and res1 != resolution:
        pub_cols, pub_from = pub2, "arm2"
    if pub1 or pub2:
        # The medians below are over all compared samples, controls included, while the
        # published comparisons cover the 369 insects only.
        note += (f" also_in_published={len(in_pub & set(shared))}"
                 f" published_columns=arm1:{res1 if pub1 else 'NA'} arm2:{res2 if pub2 else 'NA'}"
                 f" published_shared_columns=from_{pub_from}")
    # The published columns come from files written off a table this script never sees; if it is
    # not the same table, every published number belongs to another comparison. Tolerance of one
    # read: compare_published.py writes 1 for a sample with no reads at all, and the two totals
    # are float sums of the same integers in a different order.
    for arm, pub, table in (("arm1", pub1, t1), ("arm2", pub2, t2)):
        bad = [s for s, r in pub.items() if s in shared
               and abs(float(r["reads_replica"]) - sum(v[s] for v in table.values())) > 1]
        if bad:
            sys.exit(f"{arm}: {len(bad)} samples have different read totals in --{arm}'s table and in "
                     f"its --vs-published file ({', '.join(sorted(bad)[:3])}); the two were not "
                     f"computed from the same taxon table")

    rows, bcs, agree, diffs = [], [], 0, []
    for s in shared:
        p = collections.Counter({k: v[s] for k, v in t1.items() if v[s]})
        q = collections.Counter({k: v[s] for k, v in t2.items() if v[s]})
        bc, rho = bray_curtis(p, q)
        r1, r2, rp = pub1.get(s, {}), pub2.get(s, {}), pub_cols.get(s, {})
        top1 = max(sorted(p), key=p.get) if p else "NA"          # sorted first: ties break by name
        top2 = max(sorted(q), key=q.get) if q else "NA"
        reads1, reads2 = sum(p.values()), sum(q.values())
        # Exactly the samples the median Bray-Curtis beside it is over: a sample one arm lost
        # scores NA there, so counting it as -100% here would mix two denominators in one row.
        if bc is not None:
            diffs.append(100.0 * (reads2 - reads1) / reads1)
        bcs.append(bc)
        agree += top1 == top2 != "NA"
        rows.append(dict(
            sample=s, in_published=(s in in_pub) if (pub1 or pub2) else "NA",
            reads_arm1=int(reads1), reads_arm2=int(reads2),
            reads_published=rp.get("reads_published", "NA"),
            unclassified_arm1=int(drop1[s]) if b1 else "NA",
            unclassified_arm2=int(drop2[s]) if b2 else "NA",
            taxa_arm1=len(p), taxa_arm2=len(q),
            bc_arm1_vs_arm2=fmt(bc), spearman_arm1_vs_arm2=fmt(rho),
            bc_arm1_vs_published=r1.get("bray_curtis_similarity", "NA"),
            bc_arm2_vs_published=r2.get("bray_curtis_similarity", "NA"),
            top_arm1=top1, top_arm2=top2, top_published=rp.get("top_published", "NA"),
            top_arm1_agrees_arm2=top1 == top2 != "NA",
            top_arm1_agrees_published=r1.get("top_agrees", "NA"),
            top_arm2_agrees_published=r2.get("top_agrees", "NA")))

    ok = [b for b in bcs if b is not None]
    tail = (f"median_bc_arm1_vs_arm2={fmt(statistics.median(ok) if ok else None)} "
            f"dominant_taxon_agrees={agree}/{len(shared)} "
            f"samples_scored={len(ok)} (a sample empty in one arm scores NA)")
    write_tsv(a.out, PER_SAMPLE, rows, note, tail)

    summary = [summary_row("arm1_vs_arm2", a.marker, resolution, bcs, agree, len(shared), diffs)]
    for name, pub, res in (("arm1_vs_published", pub1, res1), ("arm2_vs_published", pub2, res2)):
        if not pub:
            continue
        v = [float(r["bray_curtis_similarity"]) for r in pub.values()]
        d = [100.0 * (float(r["reads_replica"]) - float(r["reads_published"])) / float(r["reads_published"])
             for r in pub.values() if float(r["reads_published"])]
        summary.append(summary_row(name, a.marker, res, v,
                                   sum(1 for r in pub.values() if r["top_agrees"] == "True"),
                                   len(pub), d))
    write_tsv(a.summary or stem + ".summary.tsv", SUMMARY, summary)

    # Per site group: the library gap is a garden-bed effect, so a study-wide median hides it.
    groups = collections.defaultdict(list)
    for r in rows:
        m = re.match(r"[A-Za-z]+", r["sample"])
        groups[m.group(0) if m else "NA"].append(r)
    gr = []
    for g, rs in sorted(groups.items()):
        #
        # Median of one column over this group, skipping the unscored samples
        def med(col):
            v = [float(r[col]) for r in rs if r[col] not in ("NA", "nan")]
            return fmt(statistics.median(v) if v else None)
        gr.append(dict(group=g, n_samples=len(rs), median_bc_arm1_vs_arm2=med("bc_arm1_vs_arm2"),
                       median_bc_arm1_vs_published=med("bc_arm1_vs_published"),
                       median_bc_arm2_vs_published=med("bc_arm2_vs_published")))
    write_tsv(stem + ".by_group.tsv",
              ["group", "n_samples", "median_bc_arm1_vs_arm2",
               "median_bc_arm1_vs_published", "median_bc_arm2_vs_published"], gr, note)

    print(f"{a.marker} ({a.rank}): {len(shared)} samples, median arm1-vs-arm2 "
          f"{fmt(statistics.median(ok) if ok else None)}, dominant taxon agrees {agree}/{len(shared)}",
          file=sys.stderr)
    return 0


#
# Command line: one taxon table per arm, plus whatever is needed to roll up or translate
def main():
    ap = argparse.ArgumentParser(
            description="""Compares the two arms' taxon tables with each other and with the published ones""",
            formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arm1", help="arm 1 <marker>.taxon_table.tsv")
    ap.add_argument("--arm2", help="arm 2 <marker>.taxon_table.tsv")
    ap.add_argument("--lineage", nargs="*", default=[],
                    help="reference <marker>.taxonomy.tsv and/or <marker>.taxon_lineage.tsv; "
                         "pooled, so both arms roll up through the same map")
    ap.add_argument("--translate1", help="arm 1's *.name_translation.tsv, to map its labels into "
                    "the published vocabulary before comparing (--rank label only)")
    ap.add_argument("--translate2", help="likewise for arm 2")
    ap.add_argument("--rank", choices=("label", "family"), default="label",
                    help="compare the labels as they are, or both arms rolled up to family")
    ap.add_argument("--marker", default="NA", help="name written into the summary rows")
    ap.add_argument("--vs-published-arm1", help="that arm's per-sample file from compare_published.py")
    ap.add_argument("--vs-published-arm2", help="likewise for arm 2")
    ap.add_argument("--out", help="per-sample TSV; the summary, by-group and rank-shift files sit beside it")
    ap.add_argument("--summary", help="summary TSV (default: --out with .summary.tsv)")
    a = ap.parse_args()

    for need in ("arm1", "arm2", "out"):
        if not getattr(a, need):
            ap.error(f"--{need} is required")
    sys.exit(run(a))


if __name__ == "__main__":
    main()
