#!/usr/bin/env python3
#
# Reports what is in the extraction and PCR blanks (phase 4), from the per-marker ASV tables
# and the assignments phase 3 wrote. Report only: nothing here subtracts, filters or marks
# anything for removal.
#
import argparse, collections, csv, os, re, statistics, sys, textwrap

SIZE = re.compile(r";size=\d+;?")
OUT_FILES = ("controls_overview.tsv", "controls_top_taxa.tsv", "controls_asvs.tsv",
             "controls_summary.txt", "controls_params.tsv")
# The order the conditions are tried in; first match wins. Written to controls_params.tsv too.
PATTERN_ORDER = ("control_only", "control_dominant", "bleed_like", "widespread", "minor")


#
# An ASV or cluster name with vsearch's ';size=N;' stripped off
def feature_id(name):
    return SIZE.sub("", name.strip())


#
# Reads a count table either way round and returns (samples, {feature: {sample: reads}}).
# The first column NAME decides the orientation; the shape cannot, since every study has more
# features than samples. Streamed line by line, non-zero cells only: arm 1's full ITS2 table
# costs 4.0 GB resident if the split rows are held (measured).
def read_matrix(path):
    counts = collections.defaultdict(collections.Counter)
    header = features = None
    samples, seen = [], set()
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            if header is None:
                header = row
                if header[0] == "sample":
                    features = [feature_id(f) for f in header[1:]]
                    dup = [f for f, n in collections.Counter(features).items() if n > 1]
                else:
                    samples = header[1:]
                    dup = [s for s, n in collections.Counter(samples).items() if n > 1]
                if dup:
                    sys.exit(f"{path}: header names {dup[0]} twice ({len(dup)} such names)")
                continue
            if len(row) != len(header):
                sys.exit(f"{path}: line {lineno} has {len(row)} fields, the header has {len(header)}")
            name = feature_id(row[0]) if features is None else row[0]
            if name in seen:
                sys.exit(f"{path}: line {lineno} repeats the name {name}")
            seen.add(name)
            if features is None:                              # rows are features, columns samples
                for s, v in zip(samples, row[1:]):
                    if v and v != "0":
                        counts[name][s] += int(v)
            else:                                             # rows are samples, columns features
                samples.append(name)
                for feat, v in zip(features, row[1:]):
                    if v and v != "0":
                        counts[feat][name] += int(v)
    if header is None:
        sys.exit(f"{path}: empty table")
    return samples, counts


#
# One {feature: (taxon, rank)} map PER MARKER. A file belongs to the marker its basename starts
# with; one matching no marker is shared by all of them
def read_assignments(paths, markers):
    # Arm 1 cluster names collide between the markers (109,517 in the full run), and a merged
    # dict would give one marker's colliding clusters the other marker's plant.
    shared, own = {}, {m: {} for m in markers}
    for path in paths:
        stem = os.path.basename(path).split(".")[0]
        into = own[stem] if stem in own else shared
        if into is shared and len(markers) > 1:
            print(f"controls_report: {os.path.basename(path)} names none of the markers "
                  f"({', '.join(markers)}); its rows are used for all of them", file=sys.stderr)
        with open(path) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            if not rd.fieldnames:
                continue
            key = rd.fieldnames[0]
            for r in rd:
                name = feature_id(r[key] or "")
                if not name:
                    sys.exit(f"{path}: a row has nothing in the `{key}` column")
                into[name] = ((r.get("taxon") or "NA").strip(), (r.get("rank") or "").strip())
    return {m: (dict(shared, **own[m]) if shared else own[m]) for m in markers}


#
# The taxon label for a feature, or a '__' placeholder naming why there is none. Unnamed
# features keep their place, so the per-control read counts still add up to the table
def label_of(assign, feat):
    if feat not in assign:
        return "__no_assignment"
    taxon, rank = assign[feat]
    return taxon if taxon and taxon != "NA" else "__" + (rank or "unassigned")


#
# A percentage to three decimals, or NA when there is nothing to divide by
def pct(part, whole):
    return f"{100.0 * part / whole:.3f}" if whole else "NA"


#
# A number for output, never in exponent form: `%g` prints 1,234,567 as `1.23457e+06`
def num(x):
    return str(int(x)) if float(x).is_integer() else f"{x:.1f}"


#
# The largest (value, name) pair, ties going to the FIRST name, which is the order every sorted
# list here uses; plain `max` on the tuple breaks a tie the other way
def largest(pairs, default=(0, "none")):
    best = min(((-v, k) for v, k in pairs), default=(-default[0], default[1]))
    return -best[0], best[1]


#
# Describes one row of counts against the blanks and the real samples, and names the first
# condition it meets:
# control_only     seen in no real sample
# control_dominant more abundant in the blanks together than in any one insect
# bleed_like       abundant in one sample and a trace of it in the blanks
# widespread       present in at least --widespread-pct of insects
# minor            none of these
def classify(row, controls, real, total, opt):
    ctrl_reads = sum(row[c] for c in controls)
    hits = [s for s in real if row[s] >= opt.presence_min_reads]
    best_reads, best_sample = largest(((row[s], s) for s in hits), (0, "NA"))
    max_rra, rra_sample = largest(((row[s] / total[s], s) for s in hits if total[s]), (0.0, "NA"))
    ratio = ctrl_reads / best_reads if best_reads else None
    prevalence = 100.0 * len(hits) / len(real) if real else 0.0
    # Bleed is a statement about ONE sample: both halves must hold in the SAME sample, or it is
    # not a bleed story. Multiplied out, so an empty sample cannot divide.
    bleed = any(total[s] and row[s] >= opt.bleed_rra * total[s]
                and ctrl_reads <= opt.bleed_ratio * row[s] for s in hits)
    if not hits:
        pattern = "control_only"
    elif ctrl_reads > best_reads:
        pattern = "control_dominant"
    elif bleed:
        pattern = "bleed_like"
    elif prevalence >= opt.widespread_pct:
        pattern = "widespread"
    else:
        pattern = "minor"
    return dict(ctrl_reads=ctrl_reads, hits=hits, max_rra=max_rra, rra_sample=rra_sample,
                best_reads=best_reads, best_sample=best_sample, ratio=ratio,
                prevalence=prevalence, pattern=pattern)


#
# Everything reported for one marker: the ASV rows that touch a blank, the per-control
# overview, and the top taxa per control. Called once per --asv-table
def analyse(marker, path, assign, control_names, opt):
    samples, counts = read_matrix(path)
    present = set(samples)
    in_table = [c for c in control_names if c in present]
    named = set(control_names)
    real = [s for s in samples if s not in named]

    total = collections.Counter()
    for c in counts.values():
        total.update(c)
    marker_reads = sum(total.values())

    # One pass for everything per feature; a second walk per (control, taxon) pair costs minutes
    # on arm 1's million clusters and gives the same numbers.
    label = {f: label_of(assign, f) for f in counts}
    tax = collections.defaultdict(collections.Counter)
    tax_rank, asvs_in, with_ctrl = {}, collections.defaultdict(collections.Counter), []
    for f, c in counts.items():
        lab = label[f]
        tax[lab].update(c)
        tax_rank.setdefault(lab, assign.get(f, ("", ""))[1])
        if any(c[x] for x in in_table):
            with_ctrl.append(f)
            for x in in_table:
                if c[x]:
                    asvs_in[x][lab] += 1

    asv_rows = []
    for f in sorted(with_ctrl, key=lambda f: (-sum(counts[f][c] for c in in_table), f)):
        d = classify(counts[f], in_table, real, total, opt)
        real_reads_f = sum(counts[f][s] for s in real)
        asv_rows.append([marker, f, label[f], assign.get(f, ("", ""))[1],
                         str(d["ctrl_reads"]),
                         ";".join(f"{c}:{counts[f][c]}" for c in in_table if counts[f][c]),
                         max((pct(counts[f][c], total[c]) for c in in_table), default="NA",
                             key=lambda x: float(x) if x != "NA" else -1.0),
                         str(len(d["hits"])), pct(len(d["hits"]), len(real)), str(d["best_reads"]),
                         d["best_sample"], f"{100.0 * d['max_rra']:.3f}", d["rra_sample"],
                         str(real_reads_f),
                         pct(d["ctrl_reads"], d["ctrl_reads"] + real_reads_f),
                         f"{d['ratio']:.5f}" if d["ratio"] is not None else "NA", d["pattern"]])

    # Every taxon the blanks hold, not only the ones listed per control, so a taxon spread over
    # two blanks is counted once.
    in_controls = {lab: classify(row, in_table, real, total, opt)
                   for lab, row in tax.items() if sum(row[c] for c in in_table)}

    overview, top_rows = [], []
    real_totals = sorted(total[s] for s in real)
    median_real = statistics.median(real_totals) if real_totals else 0
    for c in control_names:
        if c not in present:
            # Reported, never skipped -- but as NA, not as 0: a blank missing from the table was
            # either emptied by filtering or never matched a column, and only the read-tracking
            # table can tell those apart.
            overview.append([marker, c, "no", "NA", "NA", "NA", str(len(real)), num(median_real),
                             "NA", "NA", "NA", "NA", "NA"])
            continue
        per_tax = {lab: row[c] for lab, row in tax.items() if row[c]}
        reads = sum(per_tax.values())
        n_asvs = sum(asvs_in[c].values())
        top = sorted(per_tax.items(), key=lambda kv: (-kv[1], kv[0]))
        overview.append([marker, c, "yes", str(reads), str(n_asvs), pct(reads, marker_reads),
                         str(len(real)), num(median_real), pct(reads, median_real),
                         str(len(per_tax)), top[0][0] if top else "NA",
                         str(top[0][1]) if top else "0", pct(top[0][1], reads) if top else "NA"])
        for lab, n in top[:opt.top]:
            d = in_controls[lab]
            top_rows.append([marker, c, lab, tax_rank.get(lab, ""), str(asvs_in[c][lab]),
                             str(n), pct(n, reads), pct(len(d["hits"]), len(real)),
                             f"{100.0 * d['max_rra']:.3f}", d["rra_sample"],
                             str(d["best_reads"]), d["best_sample"], str(d["ctrl_reads"]),
                             f"{d['ratio']:.5f}" if d["ratio"] is not None else "NA", d["pattern"]])
        # The remainder keeps every column's meaning, so `asvs` still sums to this control's ASV
        # count in controls_overview.tsv and the reads are carried rather than dropped.
        rest = sum(n for _, n in top[opt.top:])
        if rest:
            top_rows.append([marker, c, "__other_taxa", "",
                             str(sum(asvs_in[c][lab] for lab, _ in top[opt.top:])), str(rest),
                             pct(rest, reads), "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"])

    ctrl_reads = sum(total[c] for c in in_table)
    return dict(rows=asv_rows, overview=overview, top=top_rows, marker=marker,
                controls=in_table, absent=[c for c in control_names if c not in present],
                marker_reads=marker_reads, ctrl_reads=ctrl_reads, n_real=len(real), median_real=median_real, taxa=in_controls,
                biggest=largest((total[c], c) for c in in_table),
                top_taxon=largest((d["ctrl_reads"], lab) for lab, d in in_controls.items()))


#
# The prose summary, one sentence per marker plus the standing caveats
def paragraph(results, control_names, opt):
    bits = []
    for r in results:
        if not r["controls"]:
            # Not "the cleanest result a blank can give": a misspelt --controls, an unsequenced
            # blank and four empty ones all look alike from the table alone.
            bits.append(f"For {r['marker']} not one of the {len(control_names)} named blanks "
                        f"({', '.join(control_names)}) appears in the table of {r['n_real']} "
                        f"samples, so this marker has no blank to report; that is what "
                        f"{len(control_names)} empty blanks look like and equally what a name "
                        f"matching no sample looks like, "
                        f"and the two are told apart in the read-tracking table, not here")
            continue
        kinds = collections.Counter(d["pattern"] for d in r["taxa"].values())
        top = r["taxa"].get(r["top_taxon"][1], {})
        bits.append(
            f"For {r['marker']} the blanks hold {r['ctrl_reads']:,} reads "
            f"({pct(r['ctrl_reads'], r['marker_reads'])}% of the marker's {r['marker_reads']:,}) "
            f"in {len(r['rows'])} ASVs carrying {len(r['taxa'])} taxa; the largest blank is "
            f"{r['biggest'][1]} with {r['biggest'][0]:,} reads, "
            f"{pct(r['biggest'][0], r['median_real'])}% of a median insect's "
            f"{num(r['median_real'])}; the most abundant thing in the blanks is {r['top_taxon'][1]} "
            f"({r['top_taxon'][0]:,} reads), which reaches "
            f"{100.0 * top.get('max_rra', 0.0):.1f}% of one real sample "
            f"({top.get('rra_sample', 'NA')}) and is present in "
            f"{top.get('prevalence', 0.0):.1f}% of insects; sorting those taxa by the first "
            f"condition each one meets, in the order they are listed in controls_params.tsv, "
            f"{kinds.get('control_only', 0)} are seen in no real sample, "
            f"{kinds.get('control_dominant', 0)} are more abundant in the blanks together than in "
            f"any single insect, {kinds.get('bleed_like', 0)} more look like index bleed from a "
            f"sample where they are abundant, {kinds.get('widespread', 0)} more are spread over at "
            f"least {opt.widespread_pct:g}% of insects, and {kinds.get('minor', 0)} are none of "
            f"these (a taxon meeting two conditions is counted only under the first, so these are "
            f"not counts of taxa meeting each condition)")
        if r["absent"]:
            bits[-1] += (f" ({', '.join(r['absent'])} do not appear in the table at all: either "
                         f"nothing survived for them or the name matches no sample, which this "
                         f"script cannot distinguish)")
    return textwrap.fill(
        ". ".join(bits) + ". Nothing was subtracted, filtered or removed: these tables are input to "
        "a decision, not the decision. `pattern` reports the first arithmetic condition a row met "
        f"(bleed = at least {100 * opt.bleed_rra:g}% of some real sample and at most "
        f"{100 * opt.bleed_ratio:g}% of that sample's count for it; widespread = present in at "
        f"least {opt.widespread_pct:g}% of real samples, counting {opt.presence_min_reads} read as "
        "presence), and a blank that is small in total can still matter if its reads sit on one "
        "taxon in one run, so read controls_asvs.tsv before deciding anything.", 96) + "\n"


#
# Analyses every marker and writes the five output files
def run(tables, markers, assign_paths, control_names, out, opt):
    assign = read_assignments(assign_paths, markers)
    results = [analyse(m, p, assign[m], control_names, opt) for m, p in zip(markers, tables)]
    os.makedirs(out, exist_ok=True)

    #
    # Writes one TSV of the given header and rows
    def write(name, header, rows):
        with open(os.path.join(out, name), "w") as fh:
            fh.write("\t".join(header) + "\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")

    write(OUT_FILES[0],
          ("marker control in_table reads asvs pct_of_marker_reads real_samples "
           "median_real_sample_reads pct_of_median_real_sample taxa top_taxon top_taxon_reads "
           "top_taxon_pct").split(),
          [r for res in results for r in res["overview"]])
    # `asvs`, `reads` and `pct_of_control_reads` are THIS control's; the next five describe the
    # taxon across the study; the last two are over all blanks summed. The highest share of a
    # sample and the highest read count need not be the same sample, so each names its own.
    write(OUT_FILES[1],
          ("marker control taxon rank asvs reads pct_of_control_reads real_prevalence_pct "
           "max_real_rra_pct max_real_rra_sample max_real_reads max_real_sample "
           "all_controls_reads all_controls_to_max_real_ratio pattern_all_controls").split(),
          [r for res in results for r in res["top"]])
    write(OUT_FILES[2],
          ("marker asv taxon rank control_reads per_control max_control_rra_pct "
           "real_samples_present real_prevalence_pct max_real_reads max_real_sample "
           "max_real_rra_pct max_real_rra_sample total_real_reads control_share_pct "
           "control_to_max_real_ratio pattern").split(),
          [r for res in results for r in res["rows"]])
    with open(os.path.join(out, OUT_FILES[3]), "w") as fh:
        fh.write(paragraph(results, control_names, opt))
    # Every judgement call, beside the tables it was made in: the `pattern` columns are
    # unreadable without these numbers and there is no other way to recover them.
    with open(os.path.join(out, OUT_FILES[4]), "w") as fh:
        fh.write("parameter\tvalue\n")
        for k, v in (("presence_min_reads", opt.presence_min_reads),
                     ("bleed_rra", opt.bleed_rra), ("bleed_ratio", opt.bleed_ratio),
                     ("widespread_pct", opt.widespread_pct), ("top_taxa_listed", opt.top),
                     ("pattern_precedence", ">".join(PATTERN_ORDER)),
                     ("pattern_scope", "named controls summed, not per control"),
                     ("controls", ",".join(control_names)), ("markers", ",".join(markers)),
                     ("subtraction_applied", "none")):
            fh.write(f"{k}\t{v}\n")
        for m, p in zip(markers, tables):
            fh.write(f"table_{m}\t{os.path.basename(p)}\n")
        for p in assign_paths:
            fh.write(f"assignments\t{os.path.basename(p)}\n")
    for res in results:
        print(f"{res['marker']}: {res['ctrl_reads']} reads in {len(res['rows'])} ASVs over "
              f"{len(res['controls'])} blanks present ({len(res['absent'])} absent), "
              f"{res['n_real']} real samples", file=sys.stderr)
        if not res["controls"]:
            # Said out loud: the quiet reading ("clean blanks") is usually wrong, and a
            # misspelt blank name is counted among the real samples rather than dropped.
            print(f"controls_report: WARNING {res['marker']} matched none of the named controls "
                  f"({', '.join(control_names)}); a misspelt name looks exactly like this and "
                  f"leaves the blank counted as one of the {res['n_real']} real samples",
                  file=sys.stderr)
    return results


#
# Command line: one ASV table per marker, the assignments, and the names of the blanks
def main():
    ap = argparse.ArgumentParser(
            description="""Reports what is in the blanks; nothing is subtracted or removed""",
            formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--asv-table", nargs="+", help="one table per marker, either orientation")
    ap.add_argument("--assignments", nargs="*", default=[], help="assignments.tsv / asv_taxonomy.tsv")
    ap.add_argument("--controls", help="comma-separated sample names, e.g. Extr_ctrl_1,H2O_1")
    ap.add_argument("--markers", help="comma-separated names for the tables (default: basename up to the first dot)")
    ap.add_argument("--out", help="output directory")
    ap.add_argument("--top", type=int, default=10, help="taxa listed per control (the rest are summed)")
    ap.add_argument("--presence-min-reads", type=int, default=1, help="reads counted as present in a real sample")
    ap.add_argument("--bleed-rra", type=float, default=0.05, help="share of a real sample above which bleed is plausible")
    ap.add_argument("--bleed-ratio", type=float, default=0.01, help="control reads / that sample's reads, below which")
    ap.add_argument("--widespread-pct", type=float, default=50.0, help="prevalence above which an ASV is study-wide")
    a = ap.parse_args()

    for need in ("asv_table", "controls", "out"):
        if not getattr(a, need):
            ap.error(f"--{need.replace('_', '-')} is required")
    markers = (a.markers.split(",") if a.markers
               else [os.path.basename(p).split(".")[0] for p in a.asv_table])
    if len(markers) != len(a.asv_table):
        ap.error("--markers must name one marker per --asv-table")
    controls = [c.strip() for c in a.controls.split(",") if c.strip()]
    run(a.asv_table, markers, a.assignments, controls, a.out, a)


if __name__ == "__main__":
    main()
