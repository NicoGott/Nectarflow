# Figures

Seven figures ship in `docs/figures/`. Four are drawn from a finished reanalysis run by
`bin/make_figures.py`, two by `figure.py` from both pipelines, and the last is the pipeline map,
rendered from a source file rather than drawn by hand. All of them are shown and read in the
[README](../README.md), and several again in [`results.md`](results.md).

## Regenerating them

```bash
# the ordination first: it is recomputed rather than taken from phase 6, see below. Its output
# goes next to the run, not to /tmp, because ordination_fit.tsv is the only record of the
# Procrustes fits and of which samples could not be placed.
bin/published_table.py --published assets/published/Supporting_Data.xlsx \
    --out results/reanalysis/07_ecology/published_taxon_table.tsv
bin/ordination.R --out=results/ordination \
    --tables=published=results/reanalysis/07_ecology/published_taxon_table.tsv,\
replica=results/reanalysis/04_consensus/arm1/consensus_taxon_table.tsv,\
reanalysis=results/reanalysis/04_consensus/arm2/consensus_taxon_table.tsv

bin/make_figures.py \
    --results    results/reanalysis \
    --metadata   results/reanalysis/07_ecology/insect_metadata.tsv \
    --published  assets/published/Supporting_Data.xlsx \
    --ordination results/ordination \
    --published-table results/reanalysis/07_ecology/published_taxon_table.tsv \
    --out        docs/figures
```

`bin/ordination.R` needs the `--conda_ecology` environment, for `vegan`. Without `--ordination`
the other figures are drawn and the ordination is skipped.

Everything comes from the reanalysis's consensus table (`04_consensus/arm2/`), the reference taxonomy
(`06_reference/refdb/`) and the insect metadata. `preference.png` additionally reads the floral
survey sheet out of the published workbook, which is the only figure that touches it. The PNGs in this repository are drawn from the run the README reports, and `ordination.png` from
that run's tables by way of `bin/ordination.R`.

The ordination is the one figure that recomputes a published result. Phase 6 ordinates every sample
it analyses, which is right for the results table and wrong for a picture: two insects, `GTF2_818_5`
and `HBC1_318_1`, carry a single taxon each and therefore sit at Bray-Curtis distance 1 from almost
every other insect. No ordination can place such a point, it is thrown to the edge of the
configuration, and it lowers the stress of whichever table carries it most. `bin/ordination.R`
drops any sample at distance 1 from 95% of the others, raises the try count from 100 to 400, writes
what it dropped, and reports the Procrustes agreement between each pair of configurations. Both
stresses are in `ordination_fit.tsv` beside the coordinates.

The last two are separate again, and are not drawn from a results directory:

```bash
export REPLICA_RUN=/path/to/results     # the directory holding faithful/ and reanalysis/
python3 explore/prep.py                  # -> explore/data/*.csv, flow.json
python3 explore/trace/arm1_read_fate.py  # -> per-read fate in the replica
Rscript  explore/trace/arm2_read_fate.R  --marker=rbcL --run=run_1 --wdir=... # and the other five
python3 explore/trace/join_arms.py       # -> explore/data/crosstab.csv
python3 figure.py                        # -> docs/figures/read_fates.svg, naming_crossing.svg
```

`figure.py` reads the tidy files `explore/prep.py` caches from a finished run of both pipelines, plus the
read-level join under `explore/trace/`. All of them take the run directory from `REPLICA_RUN`.
Those four files are committed under `explore/data/`, so `python3 figure.py` on its own redraws
both SVGs without a run; the chain above is only needed to regenerate the files themselves. Neither it nor anything under `explore/` is part of a
pipeline; they run by hand to make the two images the README carries.

The join is the slow part and the only step that recomputes anything. The replica's per-read fate is
recoverable from the intermediates the run already wrote, because read ids survive every step and
`vsearch --derep_fulllength` is an exact-sequence collapse. The reanalysis's is not: DADA2 discards the map
from read to variant, so `arm2_read_fate.R` re-runs the denoising of one marker-by-run group under
the run's own seed and parameters, keeps the map this time, and asserts that its merged total
matches the published track table before writing anything. All six groups reproduce exactly.
`join_arms.py` then puts the two pipelines side by side over the raw reads, which it can address as dense
SRA spot numbers, and its row and column totals are an independent derivation of the five fate
totals `figure.py` computes from the track tables; the generator asserts they agree.

The pipeline map is separate:

```bash
nf-metro validate assets/metro_map.mmd
# the map folds onto two rows, one subgraph each, which halves its width and so doubles the size
# of the labels at the width a README gives it. --y-spacing 100 keeps the published line from
# crossing the two arms where it would open a gap in them and look like a broken line
nf-metro render   assets/metro_map.mmd -o docs/figures/pipelines.svg --y-spacing 100
```

One SVG serves both themes: nf-metro writes its colours as CSS `light-dark()` pairs, so the map
follows whatever theme the page is read in. Nextflow will also draw the real graph, which is useful
for checking the map against the code rather than for reading:

```bash
nextflow run replica.nf    -preview -with-dag replica.mmd    --published assets/published/Supporting_Data.xlsx
nextflow run reanalysis.nf -preview -with-dag reanalysis.mmd --published assets/published/Supporting_Data.xlsx \
    --replica_results results
```

`-preview` resolves the channels without running anything. Both need `--published`, and the
reanalysis needs the replica's results as well, or the comparison and statistics stages never
appear.

## What each one shows

| file | what it shows |
|---|---|
| `composition_by_month.png` | the mean share of an insect's pollen load on each plant family, month by month, in two panels, bees and hoverflies |
| `composition_by_table.png` | which plant families each pollinator group carried, one bar per version of the table, so the paper's group contrasts and their sensitivity to the pipeline are in one picture |
| `ordination.png` | the same insects ordinated from each of the three tables, so the three arrangements can be compared. Coordinates from `bin/ordination.R` |
| `preference.png` | pollen load against flowering area, as a Jacobs selectivity index, per season and insect group. The wind-pollinated families are the built-in control: they should come out avoided |
| `pipelines.svg` | the two pipelines as a metro map, one line per arm, rendered from `assets/metro_map.mmd` |
| `read_fates.svg` | every raw read of both pipelines, meeting at a shared column of fates: named, or the reason it was discarded |
| `naming_crossing.svg` | the rank each analysis reached for the plants both detected, counted in detections rather than reads |

## Two things to know before reading them

**A family is only charged reads it can be charged.** A taxon label is rolled up to its family only
when the family or something below it is known; an order-rank label belongs to no single family in
the order and is excluded rather than assigned to one of them. The stacked bars therefore need not
reach 1.0, and the gap is the read mass that has no family. `compare_published.py` uses the same
rule, so the figures and the comparison tables agree. `composition_by_month` is drawn this way;
`composition_by_table` instead rescales each bar to the family-assignable mass, so its bars do reach
1.0 and the three columns are on one denominator.

**Read share is not pollen quantity.** Every figure here is proportions of sequences, which carry
the biases of PCR and of marker copy number. This applies as much to `preference.png`, where a
selectivity index built from read share is compared against flowering area measured in the field,
two quantities on different footings.
