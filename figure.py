#!/usr/bin/env python3
"""Draws the two Sankey-style figures the README carries, into docs/figures/.

Two figures. The first is the read budget: replica and reanalysis facing each other
across a shared column of fates, so every raw read on both sides lands in the
same five rows and the two can be read off against each other. Arm 1 stands in
for the published analysis, which it closely reproduces.

The fates are defined by PURPOSE, not by tool, because the two arms discard at
different steps with different software. Grouping by purpose is what makes a
shared column honest; grouping by tool name would not line up.

The second figure leaves reads behind entirely and counts detections, one per
insect and plant family, to show how the two analyses named the same plants.
"""
import collections, csv, html, json, os

REPO = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(REPO, "explore")

flow = json.load(open(f"{ROOT}/data/flow.json"))
a1, a2 = flow["arm1"], flow["arm2"]
S = lambda d, *k: sum(d[x] for x in k)
fmt = lambda v: f"{v:,}"
RAW = a1["raw_pairs"]

# the read-level join, one row per (replica fate, reanalysis fate)
JOIN = {(r["arm1"], r["arm2"]): int(r["reads"])
        for r in csv.DictReader(open(f"{ROOT}/data/crosstab.csv"))}

marg = collections.defaultdict(int)
for r in csv.DictReader(open(f"{ROOT}/data/marginals.csv")):
    marg[(r["side"], r["rank"])] += int(r["reads"])

RANK_LABEL = {"species": "species", "genus": "genus", "genus_group": "genus group",
              "family": "family", "above": "above family"}
TINT = {"species": 1.00, "genus": 0.76, "genus_group": 0.54, "family": 0.36, "above": 0.22}
INK1, INK2 = "#2176C7", "#D96A12"

# ------------------------------------------------- the shared column of fates
A1 = dict(named=S(a1, "rbcL_assigned", "ITS2_assigned"),
          quality=RAW - a1["trimmed_pairs"],
          primer=a1["merged"] - S(a1, "rbcL_length", "ITS2_length"),
          sequence=(a1["trimmed_pairs"] - a1["merged"])
                   + S(a1, "rbcL_singletons_dropped", "ITS2_singletons_dropped"))
# whatever is left reached the naming rules and got no name out of them
A1["naming"] = RAW - sum(A1.values())
_nc2 = S(a2, "rbcL_nonchimeric", "ITS2_nonchimeric")
A2 = dict(named=sum(marg[("arm2", k)] for k in RANK_LABEL),
          quality=S(a2, "rbcL_cutadapt_pairs", "ITS2_cutadapt_pairs")
                  - S(a2, "rbcL_filtered_pairs", "ITS2_filtered_pairs"),
          primer=RAW - S(a2, "rbcL_cutadapt_pairs", "ITS2_cutadapt_pairs"),
          sequence=(S(a2, "rbcL_filtered_pairs", "ITS2_filtered_pairs")
                    - S(a2, "rbcL_merged", "ITS2_merged"))
                   + (S(a2, "rbcL_merged", "ITS2_merged") - _nc2),
          naming=_nc2 - sum(marg[("arm2", k)] for k in RANK_LABEL))
assert sum(A1.values()) == RAW and sum(A2.values()) == RAW, (sum(A1.values()), sum(A2.values()))
# the join is built per read by explore/trace, from the run's own intermediates, so its
# margins are an independent derivation of the five totals above and must agree exactly
for _k in A1:
    assert sum(v for (a, _b), v in JOIN.items() if a == _k) == A1[_k], _k
    assert sum(v for (_a, b), v in JOIN.items() if b == _k) == A2[_k], _k
assert sum(JOIN.values()) == RAW

LAB = {"named": "used to name a plant", "quality": "read quality",
       "primer": "no marker, or too short for one", "sequence": "merging and real sequences",
       "naming": "the naming rules declined"}

FATES = [
    ("named", "used to name a plant", True,
     "BLAST, then the published rules: among the hits sharing the top bit score, one species gives a "
     "species, one genus a genus, up to three genera an “A/B/C spp.” group.",
     "BLAST, then the same shape of rule over a tie band one substitution wide, with an identity "
     "floor per marker and rank that can only ever make a name coarser."),
    ("quality", "read quality", False,
     "Trimmomatic: adapter clipping, a sliding quality window, and a 200 bp minimum.",
     "DADA2 filterAndTrim: expected errors over maxEE 2, any ambiguous base, truncation at the first "
     "Q2 base, a read under the truncation length, or a PhiX match."),
    ("primer", "no marker, or too short for one", False,
     "After merging, the read matched neither marker's primer regex, or fell under that marker's "
     "minimum length. The filter is a floor only.",
     "cutadapt anchored neither primer, or it anchored one and the trimmed pair fell under "
     "--minimum-length 50."),
    ("sequence", "merging, and what counts as a real sequence", False,
     f"FLASH found no acceptable overlap for {fmt(a1['trimmed_pairs'] - a1['merged'])} pairs; "
     f"vsearch --minsize 2 then dropped "
     f"{fmt(S(a1, 'rbcL_singletons_dropped', 'ITS2_singletons_dropped'))} reads sitting in clusters "
     "of one, which without an error model cannot be told from sequencing error.",
     f"DADA2 denoising and pair merging lost "
     f"{fmt(S(a2, 'rbcL_filtered_pairs', 'ITS2_filtered_pairs') - S(a2, 'rbcL_merged', 'ITS2_merged'))} "
     "together (the track table does not separate them), and the bimera model removed "
     f"{fmt(S(a2, 'rbcL_merged', 'ITS2_merged') - _nc2)} chimeras."),
    ("naming", "the naming rules declined", False,
     "Reads whose top hits spanned several families, so the rule gives no name, and reads below the "
     "1st-percentile bit-score cut-off.",
     "The result spanned more than one order, matched nothing, or the variant was too short to search."),
]

# ---------------------------------------------------------------- the crossing
RK = {4: "species", 3: "genus", 2: "genus_group", 1: "family"}
CR = ["species", "genus", "genus_group", "family"]
_rows = list(csv.DictReader(open(f"{ROOT}/data/cells.csv")))
XC, ONLY_P, ONLY_A = collections.Counter(), collections.Counter(), collections.Counter()
for r in _rows:
    if r["seen_by"] == "both":
        XC[(RK[int(r["their_rank"])], RK[int(r["our_rank"])])] += 1
    elif r["seen_by"] == "theirs":
        ONLY_P[RK[int(r["their_rank"])]] += 1
    else:
        ONLY_A[RK[int(r["our_rank"])]] += 1
BP = {r: sum(v for k, v in XC.items() if k[0] == r) + ONLY_P[r] for r in CR}
BA = {r: sum(v for k, v in XC.items() if k[1] == r) + ONLY_A[r] for r in CR}
assert sum(XC.values()) + sum(ONLY_P.values()) == sum(BP.values())
assert sum(XC.values()) + sum(ONLY_A.values()) == sum(BA.values())


def mix(hexc, t):
    r, g, b = (int(hexc[i:i + 2], 16) for i in (1, 3, 5))
    f = 0.14 + 0.86 * t
    return "#%02x%02x%02x" % tuple(round(255 - f * (255 - c)) for c in (r, g, b))


def band(x0, x1, t0, b0, t1, b1):
    m = (x0 + x1) / 2
    return (f"M{x0:.1f},{t0:.1f} C{m:.1f},{t0:.1f} {m:.1f},{t1:.1f} {x1:.1f},{t1:.1f} "
            f"L{x1:.1f},{b1:.1f} C{m:.1f},{b1:.1f} {m:.1f},{b0:.1f} {x0:.1f},{b0:.1f} Z")


# ================================================================== figure one
FW, FSCALE = 1620, 460 / RAW
TOP, ROWGAP, NW = 78, 34, 13
XR1, XF1, XF2, XR2 = 74, 548, 1072, 1546     # four bars: raw, fate, fate, raw
GUT = (XF1 + XF2) / 2
# The middle is a cross-tabulation, not a sequence: a read does not pass from one
# pipeline's fate into the other's, it simply has one fate in each. So the ribbons
# carry neither pipeline's colour, only whether the two agree about a read.
MATCH, DIFFER = "#DEDAD4", "#8A837B"


def fates():
    """Four bars. The replica's raw reads fan into its fates, the reanalysis's into
    its own, and the middle is the cross-tabulation of the two.

    Every ribbon is a real set of reads, the middle ones from a read by read join, so
    the band between "read quality" and "used to name a plant" is the reads
    Trimmomatic removed and the reanalysis went on to name. Nothing flows across the
    middle, because a read does not pass from one pipeline into the other; it has one
    fate in each. A split only ever happens at a node, never along a ribbon.
    """
    g = []
    rows, y = {}, TOP
    for key, *_ in FATES:
        h1, h2 = A1[key] * FSCALE, A2[key] * FSCALE
        rows[key] = (y, h1, h2)
        y += max(h1, h2, 58) + ROWGAP
    col_h = y - ROWGAP - TOP
    raw_top = TOP + (col_h - RAW * FSCALE) / 2

    for x, ink in ((XR1, INK1), (XR2 - NW, INK2)):
        g.append(f'<rect class="node" x="{x}" y="{raw_top:.1f}" width="{NW}" '
                 f'height="{RAW * FSCALE:.1f}" fill="{ink}"/>')
    for x, anc, ink, name in ((XR1, "start", INK1, "replica &#183; their pipeline rebuilt"),
                              (XR2, "end", INK2, "reanalysis &#183; current methods")):
        g.append(f'<text class="arm" x="{x}" y="24" text-anchor="{anc}" fill="{ink}">{name}</text>')
        g.append(f'<text class="armn" x="{x}" y="43" text-anchor="{anc}">'
                 f'{fmt(RAW)} raw read pairs</text>')
    g.append(f'<text class="armn" x="{GUT}" y="43" text-anchor="middle">'
             f'the same reads, as each pipeline classified them</text>')

    # the middle, drawn first so the fate bars sit on top of it. Each replica row pays
    # out in the row order of its counterparts and each reanalysis row takes in the
    # same order, which is what keeps the ribbons from doubling back on themselves.
    ORDER = [k for k, *_ in FATES]
    out = {k: rows[k][0] for k in ORDER}
    inn = {k: rows[k][0] for k in ORDER}
    for a in ORDER:
        for b in ORDER:
            n = JOIN[(a, b)]
            v = n * FSCALE
            col = MATCH if a == b else DIFFER
            g.append(f'<path class="link" d="{band(XF1 + NW, XF2 - NW, out[a], out[a] + v, inn[b], inn[b] + v)}" '
                     f'fill="{col}" opacity=".72"><title>replica: {html.escape(LAB[a])} &#183; '
                     f'reanalysis: {html.escape(LAB[b])} &#183; {fmt(n)} reads</title></path>')
            out[a] += v
            inn[b] += v

    c1, c2, deepest = raw_top, raw_top, 0.0
    for key, label, kept, why1, why2 in FATES:
        y0, h1, h2 = rows[key]
        t = 1.0 if kept else 0.40
        col1, col2 = mix(INK1, t), mix(INK2, t)
        g.append(f'<path class="flow" d="{band(XR1 + NW, XF1, c1, c1 + h1, y0, y0 + h1)}" '
                 f'fill="{col1}" opacity=".78"><title>replica &#183; {html.escape(why1)}</title></path>')
        g.append(f'<path class="flow" d="{band(XR2 - NW, XF2, c2, c2 + h2, y0, y0 + h2)}" '
                 f'fill="{col2}" opacity=".78"><title>reanalysis &#183; {html.escape(why2)}</title></path>')
        g.append(f'<rect class="node" x="{XF1}" y="{y0:.1f}" width="{NW}" '
                 f'height="{max(h1, 2):.1f}" fill="{col1}"/>')
        g.append(f'<rect class="node" x="{XF2 - NW}" y="{y0:.1f}" width="{NW}" '
                 f'height="{max(h2, 2):.1f}" fill="{col2}"/>')
        c1 += h1
        c2 += h2

        d = A2[key] - A1[key]
        more = (f"the reanalysis {'names' if kept else 'discards'} {fmt(abs(d))} "
                f"{'more' if d > 0 else 'fewer'}")
        g.append(f'<text class="fate-t" x="{GUT}" y="{y0 - 20:.1f}" text-anchor="middle">'
                 f'{html.escape(label)}</text>')
        g.append(f'<text class="fate-v" x="{XF1 + NW + 10}" y="{y0 + max(h1, 12) / 2 + 4:.1f}" '
                 f'fill="{INK1}">{fmt(A1[key])}</text>')
        g.append(f'<text class="fate-v" x="{XF2 - NW - 10}" y="{y0 + max(h2, 12) / 2 + 4:.1f}" '
                 f'text-anchor="end" fill="{INK2}">{fmt(A2[key])}</text>')
        g.append(f'<text class="fate-d" x="{GUT}" y="{y0 - 5:.1f}" text-anchor="middle">{more}</text>')
        deepest = max(deepest, y0 + max(h1, h2))
    return "\n".join(g), max(TOP + col_h, deepest) + 24


def crossing():
    CSCALE, CX0, CX1, CNW, CGAP = 430 / sum(BP.values()), 470, 1030, 14, 4
    ORD = {r: i for i, r in enumerate(CR)}
    g, top, deep = [], 56, 0.0
    nodes = {}
    for side, bar, x, ink in (("p", BP, CX0, "#7A3E9D"), ("a", BA, CX1, INK2)):
        y = top
        nodes[side] = {}
        for r in CR:
            h = bar[r] * CSCALE
            g.append(f'<rect class="node" x="{x - (CNW if side == "a" else 0)}" y="{y:.1f}" '
                     f'width="{CNW}" height="{h:.1f}" fill="{mix(ink, TINT[r])}"/>')
            nodes[side][r] = [y, y + h]
            y += h + CGAP
    cur_p = {r: nodes["p"][r][0] for r in CR}
    cur_a = {r: nodes["a"][r][0] for r in CR}
    for a in CR:
        for b in CR:
            v = XC[(a, b)] * CSCALE
            if v <= 0:
                continue
            col = SAME if a == b else (FINER if ORD[b] < ORD[a] else COARSER)
            g.append(f'<path class="flow" d="{band(CX0 + CNW, CX1 - CNW, cur_p[a], cur_p[a] + v, cur_a[b], cur_a[b] + v)}" '
                     f'fill="{col}" opacity=".66"><title>{RANK_LABEL[a]} &#8594; {RANK_LABEL[b]} &#183; '
                     f'{XC[(a, b)]:,} detections</title></path>')
            cur_p[a] += v
            cur_a[b] += v
    for cur, only, x, sgn, lab in ((cur_p, ONLY_P, CX0 + CNW, 1, "only the published analysis"),
                                   (cur_a, ONLY_A, CX1 - CNW, -1, "only the reanalysis")):
        base = max(max(v[1] for v in nodes["p"].values()),
                   max(v[1] for v in nodes["a"].values())) + (34 if sgn > 0 else 86)
        xe, y = x + sgn * 108, base
        for r in CR:
            v = only[r] * CSCALE
            if v <= 0:
                continue
            g.append(f'<path class="flow" d="{band(x, xe, cur[r], cur[r] + v, y, y + v)}" '
                     f'fill="{SOLO}" opacity=".6"/>')
            cur[r] += v
            y += v
        g.append(f'<rect class="node" x="{xe}" y="{base:.1f}" width="{sgn * 10}" '
                 f'height="{y - base:.1f}" fill="{SOLO}"/>')
        g.append(f'<text class="band-sub" x="{xe + sgn * 16}" y="{(base + y) / 2 + 4:.1f}" '
                 f'text-anchor="{"start" if sgn > 0 else "end"}">detected by {lab} &#183; '
                 f'{sum(only.values()):,}</text>')
        deep = max(deep, y)
    last = -1e9
    for side, bar, x, sgn in (("p", BP, CX0, -1), ("a", BA, CX1, 1)):
        tot, last = sum(bar.values()), -1e9
        for r in CR:
            y0, y1 = nodes[side][r]
            ly = max((y0 + y1) / 2, last + 34)
            last = ly
            if abs(ly - (y0 + y1) / 2) > 2:
                g.append(f'<path class="leader" d="M{x + sgn * (CNW if sgn > 0 else 0):.1f},'
                         f'{(y0 + y1) / 2:.1f} L{x + sgn * (CNW + 14):.1f},{ly:.1f}"/>')
            anc = "end" if sgn < 0 else "start"
            g.append(f'<text class="band" x="{x + sgn * (CNW + 18)}" y="{ly - 4:.1f}" '
                     f'text-anchor="{anc}">{RANK_LABEL[r]}</text>')
            g.append(f'<text class="band-sub" x="{x + sgn * (CNW + 18)}" y="{ly + 13:.1f}" '
                     f'text-anchor="{anc}">{bar[r]:,} detections &#183; {100 * bar[r] / tot:.0f}%</text>')
    g.append(f'<text class="arm" x="{CX0 - CNW}" y="{top - 16}" text-anchor="end" fill="#7A3E9D">'
             f'the published analysis named it</text>')
    g.append(f'<text class="arm" x="{CX1 + CNW}" y="{top - 16}" fill="{INK2}">the reanalysis named it</text>')
    return "\n".join(g), max(last, deep, max(v[1] for v in nodes["a"].values())) + 26


FINER, SAME, COARSER, SOLO = "#E8710A", "#C4C0BB", "#B07CC6", "#9E9E9E"
fbody, fbot = fates()
cbody, cbot = crossing()
# Each figure is written on its own, with the stylesheet inlined, so the README can embed it.
CSS = """<style>
  text { font-family: ui-sans-serif, -apple-system, "Segoe UI", Inter, system-ui, sans-serif; }
  .arm { font-size:13px; font-weight:600; } .armn { font-size:12.5px; fill:#9c958e; }
  .fate-t { font-size:13.5px; font-weight:600; fill:#2a2724; }
  .fate-v { font-size:13px; font-weight:600; } .fate-v.dim { fill:#9c958e; font-weight:400; }
  .fate-d { font-size:11.5px; fill:#9c958e; }
  .band { font-size:13px; font-weight:600; fill:#2a2724; } .band-sub { font-size:11.5px; fill:#6d6761; }
  .leader { fill:none; stroke:#9c958e; stroke-width:1; opacity:.5; }
  .node { stroke:#ffffff; stroke-width:1.2; }
  .fate-t, .fate-v, .fate-d { paint-order:stroke; stroke:#ffffff;
    stroke-width:3.4px; stroke-linejoin:round; }
</style>"""
for name, body, w, h in (("read_fates", fbody, FW, fbot), ("naming_crossing", cbody, 1500, cbot)):
    with open(f"docs/figures/{name}.svg", "w", encoding="utf-8") as fh:
        fh.write(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h:.0f}" '
                 f'width="{w}" height="{h:.0f}">{CSS}'
                 f'<rect width="100%" height="100%" fill="#ffffff"/>{body}</svg>')
    print("wrote docs/figures/%s.svg" % name)

print(f"  arm1 {sum(A1.values()):,}   arm2 {sum(A2.values()):,}   raw {RAW:,}   both conserve")
