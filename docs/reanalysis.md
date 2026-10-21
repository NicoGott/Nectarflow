# reanalysis: the same study, analysed the way it would be done today

`replica.nf` stays the attempt at the authors' own method. `reanalysis.nf` is the rebuilt pipeline,
offered as one alternative rather than as a better answer. Both read the same samplesheet and share
the download cache, so any difference between them comes from the analysis and not the data.

| Phase | What changes | State |
|---|---|---|
| 1 | cutadapt: cut the 6N spacer and both primers, split the markers | done |
| 2 | DADA2: denoise paired reads into ASVs, remove chimeras | done |
| 3 | Classify: BLAST against the same library, identity floors per rank, explicit rbcL/ITS2 consensus | done |
| 4 | Controls: what is in the extraction and PCR blanks | done |
| 5 | Compare: against the published tables and against the replica | done |
| 6 | Ecology: mvabund and vegan in R, on both pipelines and on the published table | done |

## Phase 1, cutadapt (`CUTADAPT`)

### What the reads actually look like

Measured on the raw FASTQ of two runs, not assumed:

| | R1 | R2 |
|---|---|---|
| rbcL | `NNNNNN` + `ATGTCACCACAAACAGAGACTAAAGC` (rbcLaf, 26 bp) | `AGGGGACGACCATACTTGTTCA` (rbcLr506, 22 bp), no spacer |
| ITS2 | `NNNNNN` + `ATGCGATACTTGGTGTGAAT` (ITS2F, 20 bp) | `CC` + `CGHYTGAYYTGRGGTCDC` (UniPlantR), no spacer |

The 6 bases in front of the forward primer vary between reads (each of A, C, G, T at every one of the
six positions, G-biased as usual for these runs); the 12 bases that follow are the rest of the primer,
identical in 99.5% of reads. The random spacer is on R1 only. In R2 the first base is already primer.

The replica follows their code, `HEADCROP:3`, which removes three of the six random bases and never removes
the primers. Three random bases in front of an otherwise identical sequence split it into up to 64
variants, which is why 16.9% of rbcL and 29.6% of ITS2 reads were dropped as study-wide singletons there.

### The change

One cutadapt run per marker on the raw pairs, with a linked adapter on each read:

- R1: `^NNNNNN<forward primer>...<reverse complement of the reverse primer>`
- R2: `^<reverse primer>...<reverse complement of the forward primer>`

The 5' part is anchored and required, the 3' part optional. The 3' part matters for ITS2, whose
amplicon (~345 bp) is shorter than 2x300, so each read runs into the far primer.
`--discard-untrimmed --pair-filter=any` keeps a pair only when both 5' primers were found, so marker
assignment and primer removal are the same decision, made on the raw pairs instead of after merging.
Two mismatches are allowed, as in the replica, which is what rescues the 2019 runs where one base inside the
rbcL forward primer is miscalled. No quality trimming happens here: DADA2 does its own filtering, and
it needs the quality profiles intact.

### Result on the 40-sample development set (4,258,563 raw pairs)

| | replica (reads with primers found, after Trimmomatic and FLASH) | reanalysis (pairs kept) |
|---|---|---|
| rbcL | 1,325,173 (31.1%) | 1,664,532 (39.1%), +25.6% |
| ITS2 | 1,925,810 (45.2%) | 2,450,361 (57.5%), +27.2% |
| unassigned | 7,197 after merging | 143,670 (3.4%) |

The reanalysis keeps about a quarter more reads per marker, because the primer check now happens before
quality filtering and merging rather than after. The two numbers are not yet comparable end to end:
DADA2 will remove some of these in phase 2, and only then can read retention be compared honestly.
The 3.4% unassigned is the bin to examine for inter-marker chimeras; nothing downstream depends
on the answer, so it is not examined here.

Checks on the trimmed output: no read still contains any of the four primers; rbcL R1 is at most
268 bp (300 - 6 - 26) and ITS2 R1 at most 274 bp (300 - 6 - 20), both exactly as expected; the shorter
ITS2 reads are the ones whose 3' primer was found and cut.


---

The numbers every phase produced are in [`results.md`](results.md).

---

## Phase 2 in detail: post-clustering curation (`CURATE`)

### Why a marker-gene table needs curating

The rDNA is not a single sequence in a genome. It sits in hundreds of copies and those copies are not
identical: [Song et al. 2012](https://doi.org/10.1371/journal.pone.0043971) pyrosequenced ITS2 from
178 plant species and found a mean of 35 variants per species. A denoiser that returns exact sequences
therefore returns several of them for one plant, correctly, and the table carries that plant as several
rows. The same point has been made for fungi, where ASVs are
[argued to be unusable as units](https://doi.org/10.1111/1755-0998.13329) in ITS-based community work.
Curating the table afterwards is the accepted answer, and it is a step in its own right rather than a
patch: denoising decides which sequences are real, curation decides which of them are separate
organisms.

### The rule

[LULU](https://doi.org/10.1038/s41467-017-01312-x) (Frøslev et al. 2017, *Nat Commun* 8:1188),
developed and validated on plant ITS2 metabarcoding from 130 Danish sites against field surveys of the
same sites. An ASV is taken to be another copy of a second ASV when all three hold:

- it is at least `min_match` per cent similar to it (default 84),
- it is found almost only where that one is found (`min_cooccurrence`, default 0.95 of its samples),
- it is never more abundant than it in a sample they share (`min_ratio`, default 1).

The more widespread ASV then takes its reads. The decision is made from co-occurrence and abundance
across the 374 libraries, never from the reference library, so curation cannot pull the tables toward
the published ones. An ASV with a distribution of its own is kept whatever it looks like, which is why
the *Centaurea* variant found in 109 insects, *Angelica* in 122 and *Rubus* in 101 all survive: they
are either a second rDNA copy common enough to stand alone or the local Welsh haplotype, and nothing
in this dataset separates those two.

The curation is run by [mumu](https://github.com/frederic-mahe/mumu) 1.1.4, Mahé's maintained C++
implementation of the method, at its defaults, which are the published ones (`--minimum_match 84
--minimum_relative_cooccurrence 0.95 --minimum_ratio 1 --minimum_ratio_type min`). It reads the
self-match list vsearch has already produced and writes the curated table, a full comparison log and
`<marker>.lulu_merges.tsv`, so every merge can be read one by one.

mumu's manual records three deliberate differences from the 2017 R package, and they are worth
knowing rather than hiding. The one that matters is partial overlap: in the R code a parent absent
from one of the child's samples gives a minimum ratio of zero and is therefore always rejected,
which cancels the effect of the co-occurrence option it is supposed to work with (lulu issue 8).
mumu allows the partial overlap the option was written for, so it merges somewhat more. It also
orders candidate parents by similarity first rather than by spread, and rejects a candidate only
*below* the similarity threshold rather than at it. `mumu --legacy` reproduces the R behaviour if
a comparison is ever wanted.

An earlier version of this pipeline carried a from-scratch Python implementation of LULU. It was
removed: the method is not ours, and a referee should not have to audit our version of it.

### What it removes

| | rbcL | ITS2 |
|---|---|---|
| ASVs in | 3,085 | 2,764 |
| ASVs out | 908 | 1,146 |
| reads moved onto another ASV | 3.8% | 6.7% |

Two thirds of the ASVs go and a twentieth of the reads move, which is the expected shape: what curation
removes is rare copies of abundant plants. Nothing is discarded: 473,270 rbcL and 1,116,501 ITS2 reads
are moved onto a parent and the marker totals are unchanged. The consensus table holds 422 rows.

The similarity threshold hardly matters, a good sign that the co-occurrence rule is carrying the
decision rather than the sequence rule: at 84, 95, 97 and 99 per cent the ITS2 table keeps 1,153,
1,451, 1,619 and 1,934 ASVs and moves 6.8%, 6.6%, 6.5% and 5.6% of reads. That sweep was run on the
development set and has not been repeated on the full run, where the 84 per cent column is 1,146 ASVs
and 6.7%. The published default is kept.

### Why the species floor is not relaxed instead

A near-match could be allowed to count as a species, on the argument that one mismatch in 300 bp is
within the variation of many species. The literature advises against it.
[San Martin et al. 2024](https://doi.org/10.7717/peerj.16567) cross-validated pollen assignment on
reference sequences held out of the database and found species-level identifications wrong in roughly
a quarter to a third of cases for ITS2, and around half for rbcL, **even at 100% identity**, while
genus-level calls were about 94% and 90% correct. The floor stays where it is. Curation reaches the same
observation another way, by merging a plant with itself without claiming to know which species it
is.

---

## Phase 4 in detail: the blanks

`bin/controls_report.py` reports what is in the extraction and PCR blanks and does nothing else:
it never subtracts, filters or marks anything for removal. That is deliberate, because a blank
holds two different things that need two different remedies, and only a reader can say which one
is on the screen.

* **Contamination** — a plant that entered the extraction buffer, the PCR mix or the lab. It sits
  in the blank at a read count of the same order as a real sample's, and it is usually spread
  across the study, because every sample met the same reagent.
* **Index bleed** — index hopping on a patterned flow cell, plus a little demultiplexing
  cross-talk. Reads of a plant that is genuinely abundant in *some other* sample of the run land
  in the blank at a small fraction of that sample's count. Here the blank measures a rate, not a
  contaminant, and subtracting its counts from every sample would remove real signal from the
  source sample.

So for each ASV seen in a control the script writes the four numbers that separate the two: its
count in each control, how many real samples carry it, its largest share of any single real
sample, and the ratio of its control count to its largest real count. An ASV in every insect with
a large control count is the first kind; an ASV that is 40% of one insect and 0.2% of that count
in the blank is the second.

### The `pattern` column

`pattern` names which arithmetic condition a row met. It is a description of the numbers, not a
verdict, and the thresholds behind it are flags, written to `controls_params.tsv` beside the
tables so that nobody can read a `pattern` column without the cut-offs that produced it.

| pattern | condition |
|---|---|
| `control_only` | never seen in a real sample at all |
| `control_dominant` | the blanks together hold more reads of it than any single real sample does |
| `bleed_like` | ≥ `--bleed-rra` of ONE real sample, and the blanks together hold ≤ `--bleed-ratio` of that same sample's count — both halves in the same sample, or it is not a bleed story |
| `widespread` | present in ≥ `--widespread-pct` of real samples |
| `minor` | none of the above |

Two things the tallies would otherwise overstate. Every condition is evaluated on the named blanks
*together* (their counts summed), never on one blank, so a row's pattern is the same in every
control it appears in. And the FIRST condition met is the one printed: a row can be widespread and
bleed_like at once and will be called bleed_like, so "N widespread" counts rows whose first match
was widespread, not rows present in that many samples. The numeric columns, not the label, are
what a reader should argue from. The precedence order is written into `controls_params.tsv` as
`pattern_precedence`, so it never has to be inferred from the counts.

The bleed test is a statement about ONE sample: the plant is abundant *there* and the blank holds
a small fraction of *that* sample's count for it. Reading the share off whichever sample maximises
it and the fraction off whichever sample has the most reads would let the two halves be met in two
different samples, which is not the stated condition and is not evidence of anything. It is
written as a multiplication so that a sample with no reads cannot divide.

### Per ASV and per taxon

Both levels are reported because they can disagree. The replica splits one plant across hundreds of
clusters (the `HEADCROP:3` random bases), which dilutes every per-ASV ratio, and only the taxon row
puts those reads back together. The reanalysis's ASVs are not fragmented, so there the two levels agree.

In `controls_top_taxa.tsv`, `asvs`, `reads` and `pct_of_control_reads` are *this* control's; the
four columns after them describe the taxon across the whole study; and the last two are computed
on all the blanks summed. The column names say which, because a reader who divides that row's
`reads` by `max_real_reads` and gets a different number from the ratio column would be right to
distrust the file. Taxa past `--top` are summed into an `__other_taxa` row that keeps the same
meaning in every column it fills, so `asvs` is still ASVs there and the column still sums to the
control's ASV count in `controls_overview.tsv`; the reads are carried, not dropped.

Features with no assignment row keep their place in the accounting as `__no_assignment` instead of
vanishing: a blank made of sequences the library cannot name is a different story from a blank full
of *Trifolium*, and the per-control read counts must still add up to the table.

### Traps worth knowing about

* **A blank absent from the table reads `in_table = no` and NA, never 0.** A blank missing from the
  table is not a measurement of zero: DADA2 drops a sample once nothing survives filtering, and a
  name that matches no column of the table is absent in exactly the same way. The script cannot
  tell those apart, and writing 0 would let "our blanks were clean" be read off a sample that may
  never have reached the pipeline. Whether it is a clean blank is settled against the read-tracking
  table, not here. For the same reason, a marker matching none of the named controls produces a
  loud warning: a misspelt `--controls` name looks exactly like four genuinely empty blanks, and
  the misspelt blank is left counted among the real samples.
* **One assignment map per marker, not one map for every file.** The replica names a cluster after the
  read it was seeded from (`Extr_ctrl_1|6`), and that name is a *different* cluster in each marker:
  109,517 of them collide between rbcL and ITS2 in the full run, the blanks' own clusters among
  them. Merged into one dict, the file read last wins and every colliding cluster of the other
  marker gets the wrong plant — `Extr_ctrl_1|6` is *Ranunculus* spp. in rbcL and *Veronica
  serpyllifolia* in ITS2. A file belongs to the marker its basename starts with
  (`rbcL.assignments.tsv` → rbcL); a file matching no marker goes to a pool shared by all of them,
  which is what a single combined file should do.
* **The table is read one line at a time.** The replica's full ITS2 cluster table is 762 MB of
  998,860 × 374 fields, and holding it as a list of split rows costs 4.0 GB of resident memory
  (measured) before a single count is added up — more than a Nextflow task is usually given. Only
  the non-zero cells are kept, well under 1% of them.
* **Malformed input stops the run.** A row whose field count differs from the header, or a name
  used twice, is an error rather than something to zip away or silently sum: either one means the
  table is not the table the script was pointed at. The first column's *name* decides the
  orientation; guessing from the shape would break on a study with fewer samples than ASVs, which
  is every study.
* **Ties go to the first name** in the sorted order every listing here uses. Plain `max` on a
  `(value, name)` tuple breaks a tie the other way, so the same data could name one sample in a
  table and another in the paragraph beside it.
* **Read counts are never printed in exponent form.** `%g` prints a median of 1,234,567 reads as
  `1.23457e+06`, which is not a number a reader can put in a table.
* **Everything per feature is computed in one pass** — its label, the taxon × sample matrix, which
  features a blank holds at all, and how many ASVs of each taxon each blank holds. Asking the last
  two by walking the whole feature set again for every (control, taxon) pair takes minutes on arm
  1's million clusters and gives the same numbers.
* **The summary paragraph quotes only numbers the tables also carry**, and says what was not done.
  The tempting sentence for a marker with no blanks in the table — "the cleanest result a blank can
  give" — is not supportable, and is not written.

## Phase 5 in detail: comparing the pipelines

`bin/compare_arms.py` compares the two pipelines' taxon tables with each other and with the published
ones. The replica and the reanalysis answer the same question about the same insects with different machinery, so
the honest headline is not "which arm is closer to the published tables" but "where do they differ,
and does the difference survive a roll-up to family". The script produces both, per marker:

* **the replica vs the reanalysis, computed here**: Bray-Curtis similarity of relative read abundance per insect
  (1 = identical composition), Spearman over the taxa present, and the dominant taxon of each pipeline.
* **the replica vs published and the reanalysis vs published, NOT computed here.** Those numbers are read out of
  the per-sample files `compare_published.py` already wrote (`--vs-published-arm1/-arm2`), so both
  pipelines' agreement with the Dryad sheets comes from exactly one code path, and this file never has
  to open the workbook. `compare_published.py` is the only reader of the published read counts;
  `published_table.py` opens the workbook too, for the combined table.

Four files are written per run: the per-sample table at `--out`, plus `.summary.tsv`,
`.by_group.tsv` beside it.

### Scope, which every number is conditional on

The compared composition is the reads each pipeline put on a taxon. Reads an arm could not classify are
not a composition and are excluded, so a high similarity says the classified parts agree, not that
the pipelines classified the same amount. Pass each pipeline's `<marker>.taxon_table.full.tsv` instead of
`<marker>.taxon_table.tsv` and the accounting rows (`__unassigned_*`, `__order_spanning`,
`__asv_too_short`) are reported per sample in `unclassified_arm1/2` while still being kept out of
the composition; with the plain table those columns read NA, meaning "this table does not say", not
"nothing was dropped". The replica's table has no such rows at all: what it discarded is in its own
`table_stats.tsv`, and the reanalysis's in `<marker>.taxon_table.full.tsv` and `<marker>.rank_summary.tsv`.
A count of accounting rows says nothing about how much data they hold, so the header reports both
the row count and the reads on them.

`--rank family` rolls BOTH pipelines' labels up to family before comparing, because the reanalysis deliberately
reports coarser names than the replica (its identity floors demote a third of rbcL and over 40% of ITS2
species reads), and a label-level number alone would read that design decision as disagreement.
`docs/results.md` reports the same pair for the replica against the published tables: 0.988 at label
level, 0.997 at family level (rbcL; ITS2 is 0.992 and 0.994). The roll-up uses the reference lineages — the same library in both
pipelines, so one pooled map makes the two roll-ups identical rather than merely similar — and a label
with no family in the lineages keeps its own label and is counted in the header comment rather than
dropped.

Deliberately NOT done: rolling each pair down to the coarser of our label and theirs. That would
define the relabelling using the answer key and flatter us by construction.

### The published columns and their resolution

The roll-up to family in the published comparison is done by `compare_published.py`, not by this
script, so its resolution must be read off the file and never assumed from this script's `--rank`.
The pipeline runs the replica's published comparison at label level and the reanalysis's at both, so at
`--rank family` the two columns come from different resolutions, and a family-level median (replica:
0.997) reported as a label-level one (0.990) overstates agreement by exactly what the roll-up buys.
`compare_published.py` does not record its `--rank` inside the file, so `published_resolution()`
reads it off the `.family.`/`.label.` tag the pipeline puts in the file name, and the header names
each published column's resolution beside the numbers. If the two ever differ they are not a
like-for-like pair and must not be reported as one.

`reads_published` and `top_published` are ONE column each, so they can only describe ONE of the two
published files. Taking whichever answered first would put the replica's label-level top taxon in the
same row as the reanalysis's family-level `top_agrees` the moment the resolutions differ — a row that is true
of neither file. The script prefers the file written at this run's own resolution and records which
one in `published_shared_columns=from_arm1|from_arm2`. Agreement with the published *top taxon* is
only meaningful at the resolution `compare_published.py` was run at, because our roll-up changes our
labels but not theirs.

The published columns are copied from files written from a table this script never sees. If those
are not the same tables, every published number belongs to another comparison, so the read totals
must agree sample by sample before anything is written out. The tolerance is one read:
`compare_published.py` writes 1 for a sample with no reads at all, and the two totals are float sums
of the same integers in a different order.

### Traps worth knowing about

* **A sample with no reads in one arm scores NA, not 0.** The four blanks are kept in every table on
  purpose (phases 4 and 5 want to see them), and a blank that both pipelines correctly left empty must
  not be counted as perfect disagreement.
* **Samples are compared over the UNION of the two pipelines, not the intersection.** `build_tables.py`
  and `dada2_chimera.R` both drop a sample with no surviving reads, so intersecting would delete
  exactly the samples one arm lost — which is a finding, not a nuisance. Such a sample scores NA and
  is visible as a zero read count.
* **`n_scored` sits beside `n_samples`** because the medians and the counts below them are over the
  scored samples only, while `n_samples` counts every sample compared; one denominator for both
  would be wrong. `median_read_count_diff_pct` is likewise over the scored samples only: counting a
  sample one arm lost entirely as −100% would put two numbers with two different denominators in one
  summary row (on the 40-insect dev set against the replica's full run that reads −100.00, which is a
  statement about which samples were run, not about reads). Every median is over all compared
  samples, controls included, while the published comparisons cover the 369 insects only, so the two
  denominators are stated separately in the header.
* **`--out` without a `.tsv` suffix.** The three extra files are derived from a stem, not from
  `a.out.replace(".tsv", ...)`, which does NOTHING when `--out` has no `.tsv` suffix: all four tables
  would land on `--out`, each overwriting the last and leaving only the rank shift, and it would also
  rewrite a `.tsv` occurring in a parent directory name.
* **The rank-shift table joins on family, never on labels.** Labels cannot be joined between the pipelines
  — that is the whole point — and family is the one key that survives relabelling. It is always
  computed from the labels as the pipelines wrote them, never from the rolled-up tables: after a roll-up
  every label is a family by construction, and a table saying so would hide the very difference it
  exists for. `rank_of()` reads a label's rank off its shape (`Genus epithet` / `Genus spp.` /
  `A/B spp.` / a family or higher name), which is valid only because the reanalysis keeps the replica's vocabulary.
* **Results are also reported per site group**, as `docs/results.md` does for the replica: the library gap
  is a garden-bed effect, so a study-wide median hides exactly the thing worth seeing.
* **Sample columns that normalise to one name are refused, not summed.** The duplicate second run of
  `GTF3_919_2` keeps a suffixed name precisely so that two columns do not merge into one sample
  without anyone seeing it.
* **A row that is not exactly as wide as the header stops the run.** Zipping it against the header
  would pad the missing columns with nothing, and "this taxon has no reads in those insects" is a
  composition, not an error anyone would notice.
* **Two rows for one taxon are SUMMED**, which is what `compare_published.py` does with the same
  file. Dropping one would make the two scripts disagree about the table's read total, and the
  cross-check above would then blame the mismatch on the wrong thing. How many rows were merged goes
  in the header, so it is not silent.
* **`family_map()` has three sources in increasing authority**: a lineage element is mapped to the
  first `-aceae` name in its own lineage, but only when every lineage carrying that element agrees
  (which keeps a tribe like *Anthemideae* resolvable while an order like *Asterales* correctly stays
  unresolved); then the genus and family columns; then an arm's own `taxon_lineage.tsv`, keyed by the
  whole label.
* **Dominant-taxon ties break by name** (the candidates are sorted before `max`), so the same data
  cannot name one taxon here and another elsewhere.

## Phase 3 in detail: classification

Reference for `bin/assign_taxonomy_asv.py`, which points here for everything that justifies its
decision procedure.

### Why the reanalysis does not reuse the replica's rule

The replica (`assign_taxonomy.py`) follows the paper: only the hits sharing the EXACT top bit score
vote, and a cluster whose top bit score falls below the 1st percentile is discarded. Both rules are
relative to the data in front of them, which is what the reanalysis changes. An exact tie means one
matching base can separate a species name from a genus one, so the band here is one substitution
wide rather than zero. And a percentile is a property of the sample set, not of the sequences: on
the full replica run the same written rule landed at 719.0 bits for rbcL and 484.0 for ITS2, 81.6%
and 66.8% of each marker's highest top score, and set aside 5,058 and 9,671 clusters. The
reanalysis uses an absolute identity floor per marker and rank instead, so the same sequence gets
the same verdict in any run. Measured from `faithful/07_assignment/*.assignments.tsv`.

### The tie band is stated in substitutions, not per cent

Candidates are the hits within `--band-substitutions` substitutions of the best RAW score. Stating
the band in substitutions is what makes it mean the same thing on a 458 bp rbcL and a 307 bp ITS2
alignment: one substitution is exactly `reward + penalty` raw points. A search run without a raw `score`
column is refused rather than silently falling back, so a table cannot be produced under a band
that was not applied.

### Where the identity floors come from

The floors (`species`/`genus`/`family`/`order`, per marker in `PARAMS`) are set from the published
discrimination of the two markers, not from a calibration run here. Barcode UK (Jones, Twyford et al. 2021)
reports that BLAST assigns 64% of species and 98% of genera correctly on rbcL, against 85% and 100%
on ITS2, and that ITS2 covers only 75% of species where rbcL covers 96%. The floors are stricter on
rbcL at the species rank for that reason: the marker separates species less well, so the library can
support fewer species-level calls.

They are deliberately a judgement, and the repository does not pretend otherwise: **nothing here
calibrates them.** One vsearch self-search over the built library would produce the
identity-against-correctness table, and if the floors are ever to be more than well-reasoned it
should be published beside them.

Every floor is exposed as a flag, so a sensitivity run is one flag away, and the values actually in
force are written to `<marker>.classify_params.tsv` so a table cannot be read under cut-offs it was
not made with.

`--min-asv-len` rests on the same kind of fact: the rbcL insert is fixed at 458 bp and ITS2 is about
307 bp, so a much shorter ASV is a phase-2 artefact rather than a plant. It is reported with rank
`asv_too_short` and its reads stay in the accounting.

### Empty taxonomy fields are never agreement

An empty rank is never read as agreement. No record of the current pinned library is missing an
order, but the guard does not depend on that: step 4 treats an empty field as absence of
information rather than as a match, so a library built with gaps cannot manufacture a call. (An
earlier version of this document counted 191 rbcL and 83 ITS2 such records; that was the old
live-search library and the number no longer applies.) In practice: an ASV whose candidates simply have no order recorded is
written `unassigned` with the note "no order in the reference taxonomy" and accounted under
`__unassigned_no_taxonomy`, not as `order_spanning`. `__order_spanning` must stay what phases 4–5
read it as — a tie across REAL orders that survived chimera removal.

### The subject cap

The search is `blastn -max_target_seqs 500` with **no identity or e-value filter**: everything
blastn finds is returned and the thresholds are applied afterwards, in one place, where they are
recorded. The cap could only bite if a hit the search never returned would have joined the tie
band, and blastn returns hits in descending score order, so the five-hundredth subject is far
outside a band one substitution wide (three raw score points, at --match-reward 1 and
--mismatch-penalty 2). ASSIGN_ASV is passed the same `max_target_seqs` the search
used, so the two cannot drift, and it refuses to run on a search that returned more subjects than
that.

No `-perc_identity` is passed anywhere in this pipeline, so the cap is the only thing bounding the
hit list. A `cap_may_truncate_candidates` warning flag was tried here and removed: it could
essentially never fire, and a flag that cannot fire reassures without checking anything.

### Outputs and what they guarantee

Nothing is dropped silently. Unassigned and order-spanning ASVs keep a row in
`<marker>.asv_taxonomy.tsv` and their reads are counted in `<marker>.taxon_table.full.tsv`;
`<marker>.taxon_table.tsv` holds the assigned taxa only, in the shape `build_tables.py` writes, so
`compare_published.py` runs on it unmodified. `<marker>.rank_summary.tsv` states in ASVs AND in reads
how coarse the answer is, so an agreement number cannot be quoted without its resolution, and
`<marker>.classify_params.tsv` records the thresholds that were in force. Both are written next to
the tables, not to a log. A gate the input could not support is now a refusal rather than a warning; what is
written as `not applied` rather than as its nominal value.

## Phase 3 in detail: two-marker consensus

Reference for `bin/consensus_markers.py`, which points here for everything that justifies the
operation it performs.

### What the published study did, and what this script repeats

The published study combined the markers by "a rule-based, objective, and conservative decision
process" (p.1487) whose rules are in Appendix S2; workbook sheet `consensus_identifications` (366
rows) records its result, mapping each marker-level label to a shared name, and their combined sheet is those labels applied
and the two markers' read counts SUMMED (reconstructing it from the two marker sheets plus that
mapping reproduces it cell for cell). This script repeats the OPERATION from our own reference
lineages and never opens their workbook: reading their mapping would manufacture the agreement we
then measure.

### The rule, and the four things it refuses to do

A label is replaced when the OTHER marker offers exactly one finest label that holds every plant this
one could be — 'Acer pseudoplatanus' when the other marker could only say 'Acer spp.', 'Achillea
spp.' when it could only say the tribe 'Anthemideae'. Then the two tables are relabelled and summed.
Nothing else moves:

* A label is never merged with another label OF ITS OWN MARKER. That is not a consensus, it is a
  coarsening of one marker's own table, and it is how this step runs away (see below).
* A name is never invented. Every consensus label is a label one of the markers actually gave.
* A demotion is never finer, and never a guess: if two labels of the other marker hold this one but
  neither holds the other ('Alcea/Althaea spp.' and 'Alcea/Tilia spp.'), nothing is done — and
  `consensus_map.tsv` records that as `ambiguous`, not as `unchanged`, because a label kept for want
  of an honest choice is not the same as a label nothing else could hold. 8 of the replica's 1,057 labels
  end there.
* A demotion is strict. Two labels can hold each other ('Helianthus sp.' and 'Helianthus spp.'), and
  replacing one with the other is a rename, not a consensus; both rows stay, listed as `equivalent`
  in `consensus_nested_labels.tsv`. The replica's full run emits eight 'Genus sp.' labels, seven of them
  beside the 'Genus spp.' label of the same genus.

Two different binomials of one genus are two plants, and a family-level label is a bin of reads that
could not be resolved rather than a plant, so neither absorbs anything. `--merge-family-labels` lifts
the second rule, which is what the authors did for 'Apiaceae', 'Poaceae' and others.

### The rejected alternative: merging by connected component

Merging by connected component is the obvious way to write this and was tried first. The transitive
closure of "these two labels could be one plant" collapses whole families. Measured on the 1,057
labels of the replica's full run, all four variants of it (linking on genus overlap or on containment, within a
marker or only across markers) put 75–81% of all reads on a family-level label and turn 177
Asteraceae labels into one row called 'Asteraceae': the replica's genus groups overlap each other densely
enough that a whole family becomes one component. Replayed over the 286 marker-level labels of the
authors' own sheet, the closure keeps only 76 of the 158 labels they carried through one to one.

The demotion rule keeps all 158, and every one of the 15 merges it makes is a merge they also made —
13 of them with exactly their membership, 2 with one extra label of the same tribe.

### Why `--merge-family-labels` exists but is not the default

With the flag, the rule reproduces 20 of the authors' 33 merges and still keeps 155 of the 158
labels. How much it moves depends entirely on how many family-level labels the arm produced, so it
must not be described from one arm: on the replica's full tables it moves 0.27% of reads (replica names
almost everything to genus or below), on the reanalysis's 40-insect development set 23.6%, and there it also
takes 78 labels down to a family or a tribe. It is the setting that agrees best with the answer key,
which is exactly why it is not the default and why the pipeline does not pass it: the rule that lets
a family-level bin swallow named genera has to stand on its own, and it does not. It is kept as a
switch so that the sensitivity of the comparison to it can be shown rather than asserted.

### The cost of the demotion rule, stated plainly

It does not reproduce the authors' judgement calls. It misses their merges of sister genera
('Alcea/Malva spp.' + 'Sidalcea spp.' -> 'Alcea/Malva/Sidalcea spp.'), which no rule over names alone
recovers, and it leaves two labels that bracket one plant differently ('Achillea/Cichorium spp.'
against 'Achillea/Taraxacum spp.') in two rows, where that plant is counted twice. Both are written
to `consensus_conflicts.tsv` rather than papered over, and the family roll-up in phase 5 is where the
two resolutions meet.

### The tables, and what the summed table is not

Counts are relabelled and summed (the primary table, because that is what their data are). Two
relative tables follow: per-sample shares of the summed counts, and the mean of the two markers' own
per-sample shares, since rbcL and ITS2 read counts are not on a common scale. A taxon seen by one
marker only is kept, as in their tables; `consensus_detection.tsv` records who saw what and which
rule made the label. Every one of the claims below is easy to make from the summed table and all four
are wrong:

* **Not a count of plants.** A row is a label, and a label can sit inside another label of the same
  table whenever both markers resolved the same genus to different depths, so one plant may hold two
  rows. `consensus_nested_labels.tsv` lists every such pair and `consensus_params.tsv` counts them
  and their reads: 362 of the replica's 784 rows and 81.5% of its reads, where nearly every abundant
  binomial also has a genus row. The number of rows is an upper bound on richness, and the family roll-up in phase
  5 is the level at which a richness statement can be made at all.
* **Not a set of independent observations.** The same pollen grain is counted once by rbcL and again
  by ITS2, with different primer bias, so a count here is reads and nothing more, and a taxon seen
  by both markers effectively votes twice. **Every test in phase 6 is run on this summed table**:
  the manyglm on its counts with `log(total reads)` as the offset, the PERMANOVA on per-sample
  shares of it.
* **Not filtered.** Blanks and the duplicate run stay in, as they must for phases 4 and 5, so a
  column here is not yet a result.

## Phase 6 in detail: the ecology script

`bin/ecology.R` takes one taxon × sample table of read counts and runs the paper's analyses on it.
It is pointed at three tables in turn — the replica's consensus table, the reanalysis's, and the authors' own — so
that the pipelines can be compared on their conclusions and not only on their read counts. The rules for
reading a sample code are in the script's header, where anyone changing the code needs them.
Everything else about it is here.

### The three tests

* **manyglm.** One negative-binomial GLM per plant taxon (`mvabund`), with `log(total reads of the
  sample)` as an offset, so the model is about the composition of a load and not about how deeply
  the sample was sequenced. It is asked for its deviances and no resampled p: restricting the
  resamples to within transect needs `resamp="case"`, which does not finish on a table of this
  width, and the PERMANOVA below already permutes within transect at a fraction of the cost.
* **PERMANOVA.** `vegan::adonis2` with `by = "margin"` on the Bray-Curtis distances between the
  samples' relative read abundances, permuting within transect. It is *not* one of the tests the
  paper's methods name; it is reported because its R2 is a share of the total dissimilarity and so
  stands on one scale whatever the size of the table, because each term is adjusted for the others
  rather than read in sequence, and because it is where the significance of each term comes from.
  The ordination is `bin/ordination.R`, which does all three tables at once.

### What is an option and what is not

`--table`, `--metadata`, `--label`, `--out`, `--perm` and `--formula` are options. The two filters
and the seed are constants at the top of the script. The formula is an option because the paper
compares two models by AIC; it defaults to the one the authors carried forward, season and year and
pollinator order, and `--formula="~ season + year + insect_group"` runs the four-group version.
Nothing else is a matter of which flags a run happened to be given, and
`run_info.tsv` writes all of them out next to the counts of what was used, so a result carries its
own settings.

### Where the metadata comes from

A sample code carries the site, the transect, the month and the year, and hence the season, and the
script derives all of them itself. The one thing it does not carry is which insect the pollen came
off: bumblebees and hoverflies were collected side by side at the same transect on the same day.
That comes from `--metadata`, written by `bin/ecology_metadata.py` from the `Pollen_sample_metadata`
sheet and nothing else, so the ecology never touches the published read counts that phase 5 compares
against. Where both sides know a column the file wins, and every disagreement between the two
derivations is counted and named in `summary.txt` rather than quietly resolved — on the full run
that is one sample, `HBG1_318_1`, whose site the sheet gives as Nature Reserve while every other
HBG sample is Botanic Garden.

### Traps worth knowing about

* **The offset is the whole sample, not the modelled subset.** `total_reads` is the reads on named
  taxa in that sample, including the taxa too rare to enter the model, because that is the
  sequencing effort the composition has to be read against.
* **`__` rows are reads, not plants.** The reanalysis's `*.taxon_table.full.tsv` carries read-accounting rows
  (`__unassigned_no_hit` and friends). They are counted per sample and dropped. A table without them
  — the consensus tables, the published sheet — carries no such record at all, and the summary says
  so rather than letting a fully assigned sample be mistaken for a half-assigned one.
* **The statistics describe the table in front of the script, and nothing wider.** A disagreement
  between pipelines can come from the taxon labels of phase 5 as easily as from the ecology. Insects
  taken at the same transect are not independent replicates, which is why the PERMANOVA permutes
  within transect. Neither test can separate year from the sequencing run, which is the same column
  in this dataset.
