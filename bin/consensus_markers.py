#!/usr/bin/env python3
#
# Combines the two markers' taxon tables (phase 3 of arm 2) into one composition per insect.
# A label is demoted only when the OTHER marker offers exactly one finest label that holds it;
# both tables are then relabelled and summed.
# The rules in full: docs/reanalysis.md#phase-3-in-detail-two-marker-consensus
#
import argparse, collections, csv, functools, os, sys

MARKERS = ("rbcL", "ITS2")
# Tab-separated, not CSV: nothing that writes them quotes, so a field beginning with a double quote
# (a cultivar epithet) must be read back verbatim rather than unquoted into a different name.
TSV = dict(delimiter="\t", quoting=csv.QUOTE_NONE)


#
# The genus names a label mentions:
# 'Cirsium/Centaurea spp.' -> {'Cirsium','Centaurea'}; 'Asteraceae' -> {'Asteraceae'}
@functools.cache                       # called O(labels^2) times; the frozenset is safe to share
def genera(label):
    label = label.replace(" spp.", "").strip()
    if "/" in label:
        return frozenset(g.strip() for g in label.split("/") if g.strip())
    return frozenset(label.split()[:1])


#
# Is this label a binomial, i.e. one named species? 'Genus sp.' is not: it names a set, being a
# reference record never identified below genus rather than a call this pipeline made
def is_species(label):
    # Read as a binomial it would inflate every rank tally and hide real disagreements from
    # consensus_conflicts.tsv as 'two named species'.
    return ("/" not in label and " spp." not in label and not label.endswith(" sp.")
            and len(label.split()) > 1)


#
# The lineages of both markers' reference taxonomies, pooled: which names sit above which,
# and which family each name was ever seen under
class Lineages:
    #
    # Reads the taxonomy TSVs and indexes their lineage, genus, family and order columns
    def __init__(self, paths):
        self.genera, self.orders = set(), set()
        self.above = collections.defaultdict(set)      # name -> the names above it in a lineage
        self.families = collections.defaultdict(set)   # name -> every family it was ever seen under
        for path in paths:
            with open(path, encoding="utf-8") as fh:
                rows = list(csv.DictReader(fh, **TSV))
            for r in rows:
                lineage = [x for x in (r.get("lineage") or "").split(";") if x]
                if r.get("genus"):
                    self.genera.add(r["genus"])
                if r.get("order"):
                    self.orders.add(r["order"])
                for i, name in enumerate(lineage):
                    self.above[name].update(lineage[:i])
                    if r.get("family"):
                        self.families[name].add(r["family"])

    #
    # The family a single name sits in, or None. Only when every lineage carrying the name
    # agrees, which keeps families and tribes and rejects 'rosids'
    def family_of_name(self, name):
        # Orders are excluded by name: an order holding a single family in this library must
        # not inherit it.
        fs = self.families.get(name, ())
        return next(iter(fs)) if len(fs) == 1 and name not in self.orders else None

    #
    # The family a whole label sits in, or None when its genera disagree
    def family_of(self, label):
        fs = {self.family_of_name(g) for g in genera(label)}
        return fs.pop() if len(fs) == 1 and fs != {None} else None

    #
    # Does the single name x appear above every genus of y in some lineage?
    def is_ancestor(self, x, y):
        names = genera(x)
        if len(names) != 1:
            return False
        n = next(iter(names))
        return bool(genera(y)) and all(n in self.above.get(g, ()) for g in genera(y))

    #
    # Does x hold every plant y could be? A binomial holds only itself; a genus or genus group
    # holds the genera it names; a tribe holds what sits under it. A family-level label holds
    # nothing unless families_too, and an order-level label never holds anything
    def contains(self, x, y, families_too=False):
        # A family label is an unresolved bin: letting it absorb pulls every genus of its
        # family into one row.
        if x == y:
            return True
        if is_species(x):
            return False
        if genera(x) <= self.genera:                      # x names genera
            return bool(genera(y)) and genera(y) <= genera(x)
        name = next(iter(genera(x)), None)                # x names a rank above genus
        family = self.family_of_name(name) if name else None
        if family is None or (family == name and not families_too):
            return False
        return self.is_ancestor(x, y)


#
# The label-merging rule, and the one place a label may change. Each label is replaced by the
# FINEST label the OTHER marker offers that strictly holds it, then that replacement may itself
# be replaced, and so on until nothing moves. Returns the final mapping and, per label, which
# of unchanged / ambiguous / genus / genus_group / supra_generic happened to it.
def demote(labels, marker_of, lin, families_too=False):
    step, ambiguous = {}, set()

    #
    # Strict containment: two labels CAN hold each other ('Helianthus sp.' and 'Helianthus
    # spp.'), and swapping one for the other is a rename, not a consensus. Strictness is also
    # what makes the chain below finite, a strict partial order having no cycles.
    def holds(x, y):
        return lin.contains(x, y, families_too) and not lin.contains(y, x, families_too)

    for l in labels:
        # Only the OTHER marker may force a demotion, and only to the FINEST label it offers
        # that still holds l: demoting further would move these reads into a bin the other
        # marker's own matching row is not in. Two incomparable holders is no honest choice.
        cand = [c for c in labels if c != l and marker_of[c] - marker_of[l] and holds(c, l)]
        fine = [c for c in cand if not any(holds(c, d) for d in cand)]
        step[l] = fine[0] if len(fine) == 1 else l
        if len(fine) > 1:
            ambiguous.add(l)   # kept, but not for the reason a label nothing could hold is kept
    # The bound is the longest chain this many labels can form; anything smaller would report a
    # long but honest chain as a cycle.
    out = dict(step)
    for _ in range(len(labels) + 1):
        nxt = {l: step[out[l]] for l in labels}
        if nxt == out:
            break
        out = nxt
    else:
        raise SystemExit("consensus_markers: the demotion chain did not settle (cycle in contains?)")
    rule = {}                                  # named after the rank of the label that replaced l
    for l in labels:
        c, gs = out[l], genera(out[l])
        rule[l] = ("ambiguous" if c == l and l in ambiguous else
                   "unchanged" if c == l else
                   "supra_generic" if not gs <= lin.genera else
                   "genus" if len(gs) == 1 else "genus_group")
    return out, rule


#
# One cell of a taxon table as a whole read count, or exit
def reads(v, where):
    # '10.0' is accepted (a writer may have round-tripped through float); '10.9' is not, since
    # a table of shares floored to integers would read as a composition of almost nothing
    # rather than as the wrong file.
    try:
        x = float(v) if v.strip() else 0.0
    except ValueError:
        raise SystemExit(f"consensus_markers: {where}: '{v}' is not a read count")
    if x < 0 or x != int(x):
        raise SystemExit(f"consensus_markers: {where}: {v} is not a whole read count")
    return int(x)


#
# Reads a taxon x sample table into (sample names, {taxon: {sample: reads}})
def read_table(path):
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.reader(fh, **TSV))
    samples = [s.strip() for s in rows[0][1:]]
    table = {}
    for n, r in enumerate(rows[1:], 2):
        if not r or not any(f.strip() for f in r):
            continue
        where = f"{path} line {n}"
        # A ragged row padded with zeros would arrive as a plausible-looking composition that
        # nothing downstream could still see.
        if len(r) != len(samples) + 1:
            raise SystemExit(f"consensus_markers: {where} has {len(r)} fields, not {len(samples) + 1}")
        table[r[0].strip()] = {s: reads(v, where) for s, v in zip(samples, r[1:])}
    return samples, table


#
# Writes a taxon x sample matrix, rows sorted by name
def write_matrix(path, first, table, samples, fmt):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(first + "\t" + "\t".join(samples) + "\n")
        for k in sorted(table):
            fh.write(k + "\t" + "\t".join(fmt(table[k].get(s, 0)) for s in samples) + "\n")


#
# Formats a share for output
def share(v):
    return "0" if not v else f"{v:.10g}"   # 10 digits: a sample's shares still sum to 1 on a reread


#
# Reads both markers' tables, applies the demotion rule, and writes the consensus table, the
# two relative-abundance matrices, the label map, the detections, the nested labels, the
# conflicts and the parameters
def main():
    ap = argparse.ArgumentParser(
            description="""Combines the rbcL and ITS2 taxon tables into one consensus table""",
            formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rbcl", required=True, help="rbcL taxon x sample table (04_tables/rbcL.taxon_table.tsv)")
    ap.add_argument("--its2", required=True, help="ITS2 taxon x sample table")
    ap.add_argument("--taxonomy", nargs="+", required=True,
                    help="reference taxonomy TSVs (accession/species/genus/family/order/lineage), both markers")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--conflict-min", type=float, default=0.01,
                    help="share of its marker's reads in a sample a label must exceed to be reported as a conflict")
    ap.add_argument("--merge-family-labels", action="store_true",
                    help="let a family-level label absorb everything of that family the other marker named, "
                         "as the authors did for Apiaceae and Poaceae")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    lin = Lineages(a.taxonomy)

    samples, marker_tables = set(), {}
    for marker, path in zip(MARKERS, (a.rbcl, a.its2)):
        cols, table = read_table(path)
        marker_tables[marker] = table
        samples |= set(cols)
    samples = sorted(samples)
    labels = sorted({t for table in marker_tables.values() for t in table})
    marker_of = {l: frozenset(m for m in MARKERS if l in marker_tables[m]) for l in labels}

    consensus, rule_of = demote(labels, marker_of, lin, a.merge_family_labels)

    # relabel, then sum, exactly as they did
    summed = collections.defaultdict(collections.Counter)
    per_marker = {m: collections.defaultdict(collections.Counter) for m in MARKERS}
    for marker, table in marker_tables.items():
        for taxon, counts in table.items():
            for s, n in counts.items():
                summed[consensus[taxon]][s] += n
                per_marker[marker][consensus[taxon]][s] += n

    # Relabelling and summing cannot lose a read; a mapping missing a label would.
    given = sum(n for t in marker_tables.values() for c in t.values() for n in c.values())
    kept = sum(n for c in summed.values() for n in c.values())
    if given != kept:
        raise SystemExit(f"consensus_markers: relabelling changed the read total ({given} -> {kept})")

    # .get from here on: per_marker is a defaultdict, and asking it for a taxon the marker never
    # saw would insert an empty row that a later pass would count as a detection
    nothing = collections.Counter()
    got = lambda m, t: per_marker[m].get(t, nothing)

    totals = {m: {s: sum(c[s] for c in per_marker[m].values()) for s in samples} for m in MARKERS}
    grand = {s: sum(totals[m][s] for m in MARKERS) for s in samples}
    write_matrix(f"{a.out}/consensus_taxon_table.tsv", "taxon", summed, samples, str)

    with open(f"{a.out}/consensus_map.tsv", "w", encoding="utf-8") as fh:
        fh.write("marker\tmarker_taxon\tconsensus_label\trule\tfamily\n")
        for marker in MARKERS:
            for taxon in sorted(marker_tables[marker]):
                fh.write(f"{marker}\t{taxon}\t{consensus[taxon]}\t{rule_of[taxon]}\t{lin.family_of(taxon) or 'NA'}\n")

    rules_of = collections.defaultdict(set)
    for l in labels:
        rules_of[consensus[l]].add(rule_of[l])
    with open(f"{a.out}/consensus_detection.tsv", "w", encoding="utf-8") as fh:
        fh.write("taxon\tsample\tdetected_by\trbcL_reads\tITS2_reads\trule\n")
        for t in sorted(summed):
            rules = ";".join(sorted(rules_of[t]))
            for s in samples:
                n = {m: got(m, t)[s] for m in MARKERS}
                if not sum(n.values()):
                    continue
                by = "both" if all(n.values()) else next(m for m in MARKERS if n[m])
                fh.write(f"{t}\t{s}\t{by}\t{n['rbcL']}\t{n['ITS2']}\t{rules}\n")

    # A consensus row can sit INSIDE another one: only the other marker may coarsen a label, so
    # if both markers report 'Borago officinalis' and both also report 'Borago spp.', neither
    # moves and one plant keeps two rows. The row count is an upper bound on the plants.
    reads_of = {t: sum(c.values()) for t, c in summed.items()}
    nested = [(t, h) for t in sorted(summed) for h in sorted(summed)
              if t != h and lin.contains(h, t, a.merge_family_labels)]
    nested_rows = {t for t, _ in nested}
    nested_reads = sum(reads_of[t] for t in nested_rows)
    with open(f"{a.out}/consensus_nested_labels.tsv", "w", encoding="utf-8") as fh:
        fh.write("taxon\tholder\trelation\treads_taxon\treads_holder\n")
        for t, h in nested:
            # 'equivalent' when the two hold each other ('Helianthus sp.' and 'Helianthus
            # spp.'): one set of plants under two names, which no rule here may rename away.
            rel = "equivalent" if lin.contains(t, h, a.merge_family_labels) else "inside"
            fh.write(f"{t}\t{h}\t{rel}\t{reads_of[t]}\t{reads_of[h]}\n")

    # What the demotion rule did not reconcile: two abundant labels of one family, one per
    # marker, in different rows - either two plants, or one plant bracketed differently and
    # counted twice. Two binomials, and labels of different families, are simply two plants.
    n_conflicts = n_unplaced = n_dominant = 0
    with open(f"{a.out}/consensus_conflicts.tsv", "w", encoding="utf-8") as fh:
        fh.write("sample\trbcL_label\tITS2_label\trbcL_reads\tITS2_reads\n")
        for s in samples:
            hot = {m: sorted(t for t, c in marker_tables[m].items()
                             if totals[m][s] and c.get(s, 0) / totals[m][s] > a.conflict_min) for m in MARKERS}
            for x in hot["rbcL"]:
                f = lin.family_of(x)
                for y in hot["ITS2"]:
                    if consensus[x] == consensus[y] or (is_species(x) and is_species(y)):
                        continue
                    g = lin.family_of(y)
                    if f is None or g is None:
                        # The lineages cannot say whether these two are even in one family;
                        # counted, so the conflict count never reads as the disagreement count.
                        n_unplaced += 1
                    elif f == g:
                        fh.write(f"{s}\t{x}\t{y}\t{marker_tables['rbcL'][x].get(s, 0)}"
                                 f"\t{marker_tables['ITS2'][y].get(s, 0)}\n")
                        n_conflicts += 1
            # The coarsest form of the same question: do the markers agree on what this insect
            # mostly carried? Sorting first keeps ties deterministic.
            top = {m: (max(sorted(marker_tables[m]), key=lambda t: marker_tables[m][t].get(s, 0))
                       if totals[m][s] else None) for m in MARKERS}
            if all(top.values()) and consensus[top["rbcL"]] != consensus[top["ITS2"]]:
                n_dominant += 1

    counts = collections.Counter(rule_of.values())
    empty = [s for s in samples if not grand[s]]
    unplaced_labels = sum(1 for l in labels if lin.family_of(l) is None)
    # The settings that qualify the tables, since no matrix carries a trace of how it was built.
    # No date and no commit: two runs of the same inputs must give byte-identical files.
    with open(f"{a.out}/consensus_params.tsv", "w", encoding="utf-8") as fh:
        fh.write("parameter\tvalue\n")
        for k, v in (("rbcl_table", a.rbcl), ("its2_table", a.its2),
                     ("taxonomy", ";".join(a.taxonomy)),
                     ("merge_family_labels", str(a.merge_family_labels).lower()),
                     ("conflict_min", f"{a.conflict_min:g}"),
                     ("samples", len(samples)), ("samples_without_reads", len(empty)),
                     ("marker_labels_rbcL", len(marker_tables["rbcL"])),
                     ("marker_labels_ITS2", len(marker_tables["ITS2"])),
                     ("marker_labels_union", len(labels)), ("consensus_rows", len(summed)),
                     ("reads_in", given), ("reads_out", kept),
                     ("labels_unchanged", counts["unchanged"]),
                     ("labels_ambiguous", counts["ambiguous"]),
                     ("labels_to_genus", counts["genus"]),
                     ("labels_to_genus_group", counts["genus_group"]),
                     ("labels_to_supra_generic", counts["supra_generic"]),
                     ("labels_without_family", unplaced_labels),
                     ("rows_inside_another_row", len(nested_rows)),
                     ("reads_in_rows_inside_another_row", nested_reads),
                     ("pct_reads_in_rows_inside_another_row",
                      f"{100 * nested_reads / kept:.2f}" if kept else "NA"),
                     ("conflict_pairs_reported", n_conflicts),
                     ("conflict_pairs_family_unknown", n_unplaced),
                     ("samples_dominant_label_disagrees", n_dominant)):
            fh.write(f"{k}\t{v}\n")

    print(f"{len(labels)} marker labels ({len(marker_tables['rbcL'])} rbcL, "
          f"{len(marker_tables['ITS2'])} ITS2) -> {len(summed)} consensus rows over {len(samples)} "
          f"samples; demoted " + ", ".join(f"{counts[r]} to {r}" for r in
                                           ("genus", "genus_group", "supra_generic"))
          + f", {counts['ambiguous']} ambiguous, {unplaced_labels} the lineages cannot place in a "
          f"family; {len(nested_rows)} rows ("
          + (f"{100 * nested_reads / kept:.1f}%" if kept else "NA")
          + " of reads) sit inside another row; {} sample-level conflicts ({} pairs unplaceable), "
            "dominant label disagrees in {} of {} samples".format(n_conflicts, n_unplaced,
                                                                  n_dominant, len(samples)),
          file=sys.stderr)


if __name__ == "__main__":
    main()
