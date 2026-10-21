# Walkthrough: what each stage does, and how to tell whether it matches

The aim of `replica` is to repeat the published analysis of Lowe et al. 2022 (J Appl Ecol) with the
authors' own method, made portable. "Their command" below is read from their code
(github.com/colford/nbgw-plant-illumina-pipeline) and from the paper's Methods and Appendix S2.
After every run, open `results/read_tracking.tsv`: one row per sample, one column per stage.

Numbers to hold the run against. From the paper (p.1488): 40,800,709 reads returned; 22,510,682
after quality control (11,305,697 rbcL; 11,204,985 ITS2); 369 insects successfully sequenced of 382
caught; 191 plant taxa from the two markers combined. From the Dryad workbook, sheet
`consensus_identifications`, which the paper does not print: 180 marker-level rbcL labels and 186
ITS2 labels.

## 0. Input
**What.** Ask the archive (ENA) for the list of runs in project PRJNA763761, download both FASTQ
files of each, verify the MD5 checksums. 374 runs: 369 insects + one insect sequenced twice +
4 blanks (`Extr_ctrl_1/2`, `H2O_1/2`).
**Check.** `raw_pairs` summed over samples should be 40.8 million.
**Exome analogy.** Same as pulling FASTQs for a cohort; nothing new.

## 1. Trim  (`TRIM`, Trimmomatic)
**Their command.** `ILLUMINACLIP:<adapters>:2:30:10 HEADCROP:3 LEADING:20 SLIDINGWINDOW:4:Q MINLEN:200`
**What.** Remove adapter read-through, cut the first 3 bases, cut low-quality ends, drop reads that
end up shorter than 200 bases. Only pairs where both mates survive go on.
**Why it matters here.** MiSeq 2x300 reverse reads decay badly towards the end. Trim too little and
merging fails on mismatches; trim too much and the mates no longer overlap.
**Quality threshold.** `Q` is asked interactively by their script and not reported in the paper,
but the prompt's suggested default is 20 (`Enter quality score [20]:` in `plant-pipeline.py`). We use 20.
**Is `HEADCROP:3` really what they ran?** It is worth asking, because the spacer is six bases and
cropping three is the odd choice: it is the one parameter that shapes everything downstream. Three
things say yes. It is **not a default that could have been left behind**: `Q` is the only prompted
value, and the 3 is written into the shell command the script generates, so running 6 would have meant
editing the source. It has been 3 in **all three commits** of their repository, 2016, 2017 and the
October 2020 release that matches the paper. And the read counts rule 6 out. Dropping the first three
bases from our own cluster centroids collapses 87.4% of rbcL and 87.6% of ITS2 clusters, in groups of
exactly 64, which is the splitting the remaining spacer causes; under `HEADCROP:6` that splitting
disappears and almost all of the reads the study-wide singleton filter removes, 16.9% of rbcL and
29.6% of ITS2, would survive instead. Retained reads would then come out up to 20% (rbcL) and 42%
(ITS2) above what we have now, where our counts already sit a median 1.6% and 2.1% **above** the
published per-insect totals. A run with 6 would overshoot the paper by a fifth to a half; the run with
3 matches it to within two per cent.
**Exome analogy.** Identical to adapter/quality trimming before BWA, but with far more at stake:
an exome tolerates a lost read, an amplicon pair that cannot merge is lost entirely.

**Deviation: the adapter file, and why it does not matter.** Their script hard-codes
`adaptors/NexteraPE-PE.fa`, but these libraries are not Nextera. In the one raw sample still on disk
the Nextera transposase sequence does not occur at all, while the standard `AGATCGGAAGAGC`
read-through occurs in 0.03% of reads, so the default here is `TruSeq3-PE-2.fa`. Both lists are
Trimmomatic's own stock files and carry its GPL-3, so neither is copied into this repository; `TRIM`
takes the one named by `--adapters` out of the installed package, and `--adapters NexteraPE-PE.fa`
runs their recipe exactly as written. The file their `adaptors/` directory holds is byte-identical
to the package copy (md5 `356e426e476e887c83707a484e52a319`). It makes no difference worth having: Trimmomatic 0.33 with their full command on that sample
keeps 63,241 pairs of 74,280 with TruSeq and 63,243 with Nextera, and the two outputs differ in two
reads. The shipped numbers come from the TruSeq default.

## 2. Merge  (`MERGE`, FLASH)
**Their command.** `flash --max-overlap 450 --min-overlap 10 --max-mismatch-density 0.25`
**What.** Overlap the forward and reverse read of each pair into one sequence covering the whole
amplicon (primer to primer).
**Check.** `merged` vs `trimmed_pairs`. rbcL amplicons are long (~500 bp + primers), so overlap is
short and this is where rbcL reads are lost.
**Exome analogy.** None. In resequencing, mates are aligned separately to a reference. Here there
is no reference to align to, so the pair itself must reconstruct the molecule.

## 3. Split by marker and length  (`SPLIT_BY_PRIMER`)
**Their code.** `merged_demultiplex_by_primer.py`, read from their repository on 5 October 2026.
Each read is searched for a forward primer fragment in its **first 50 bases** -- `ACAGAGACTAAAGC`
for rbcL, `ACTTGGTGTGAAT` for ITS2. Only if neither is found is the read searched for a reverse
fragment in its **last 50 bases** -- `TGAACAAGTATGG` for rbcL, `G[ACT]GACC[CT]CA[AG][AG]` for ITS2.
A read matching nothing goes to `unknown`. The match is exact: `re.search` on the fragment, no
mismatch allowance. Primers are **not** removed. Then keep merged reads >= 450 bp (rbcL) or
>= 300 bp (ITS2): `minimum_seq_lengths` in their script says 300, where the paper says 350.
**What.** Both PCR products of an insect were pooled under one index, so each merged read is one
marker or the other. The 6 random bases in front of the forward primer (Appendix S3, Table S1) are
why the search looks inside the first 50 bases rather than at position 1.
**Deviation: primer mismatches.** We allow 2 mismatches per primer fragment
(`--primer_mismatches`) where their `re.search` is exact. `*_primer_exact` in the tracking table
counts the reads their rule would have kept, and the answer is nearly all of them: 13,897,706 of
13,940,743 rbcL reads (99.69%) and 16,232,728 of 16,237,639 ITS2 reads (99.97%), so the allowance
rescues 0.31% and 0.03%. Only 67 merged reads in the whole study match no primer at all.

The deviation is harmless, and it is worth being clear about how little it buys. One base inside
the rbcL forward primer is sometimes miscalled, `...ACTAAAAC...` for `...ACTAAAGC...`, from a
low-diversity cycle, but that variant is in 0.10% of the reads of the one raw sample still on disk,
and the study-wide 99.69% above leaves no room for it to be common anywhere else. The allowance is
kept because two mismatches in a 13-14 base fragment cannot pull in the other marker, not because
exact matching would lose data.
**Deviation: ITS2 minimum length.** The paper says > 350 bp, their script says 300. Following
the paper would remove most of the merged ITS2 reads and could not leave the 11.2 million ITS2 reads
they report, so we follow the script (300). The read-length distribution that shows this is
regenerated by the current pipeline; quote it from there, not from memory.

**Deviation: a read that matches both markers.** Ours takes the first marker that matches, so a read
carrying both forward fragments goes to one file; their loop writes it to both. The fragments are dissimilar enough that this should be rare;
`split_by_primer` reports the counts, so the size of it is measurable rather than assumed.
**Check.** `unknown` should be a fraction of a percent. Note that the paper's "22.5 million after
quality control" corresponds to reads *assigned* at the end (stage 8), not to this stage.
**Exome analogy.** Like demultiplexing, but by amplicon rather than by sample barcode.

## 4. Dereplicate within sample  (`DEREP`, vsearch)
**Their command.** `vsearch --derep_fulllength --sizeout`
**What.** Collapse identical reads, keep the count in the label (`;size=N`).
**Exome analogy.** The opposite of duplicate marking: there, duplicates are artefacts to discard;
here, the number of identical reads *is* the measurement.

## 5. Cluster across samples, drop singletons  (`CLUSTER`, vsearch)
**Their command.** `vsearch --cluster_size <all samples> --id 1.0 --sizein --sizeout --uc ...`, then
keep clusters seen more than once.
**What.** Pool every sample; each distinct sequence in the whole study becomes one cluster. A
sequence seen exactly once anywhere is assumed to be a sequencing error.
**Deviation: exact dereplication instead of alignment clustering.** Their `--cluster_size --id 1.0`
aligns every sequence to every cluster and did not finish on 20 samples in six minutes here; the
cost grows with the square of the number of sequences. By default we pool the samples and collapse
*identical* sequences (`vsearch --derep_fulllength`), which takes seconds. The only difference: a
sequence contained inside a longer one (a truncated read) stays a separate cluster instead of being
absorbed. Their exact command remains available with `--cluster_mode align`.
**This is the step a denoiser (DADA2) replaces.** They use a count threshold; DADA2 uses an error
model. Everything else can stay identical, which makes the comparison clean.
**Exome analogy.** A hard depth filter on variant calls versus a statistical genotype likelihood.

## 6. Reference library  (`REFERENCE`)
**Theirs.** A custom BLAST database of 5,887 UK plant species, attributed by the paper to Jones,
Brennan et al. (2021) and comprising the natives of Stace (2019), the naturalised and alien species
of Preston et al. (2002) and horticultural species from the IRIS BG database at the National
Botanic Garden of Wales (p.1487). Not deposited.
**Ours.** 14,173 sequences, rebuilt from a pinned manifest: rbcL 7,716 and ITS2 6,457. Two
sources, both public, assembled to answer the same question theirs did.

- **Barcode UK** (Jones, Twyford et al. 2021), the vouchered Sanger backbone: 4,466 rbcL and 2,574 ITS2,
  downloaded from GenBank by the accession ranges in that paper. 1,420 species.
- **GenBank fills** for the species Barcode UK does not cover, 3,250 rbcL and 3,883 ITS2. Which
  species those are comes from one list that stands in for the ones their Methods name, not from
  anything they found. Their Methods name Stace (2019) for natives and Preston et al. (2002) for
  naturalised and alien species, both books; the **BSBI Plant Atlas 2020** (Stroh et al. 2024) is
  used here instead, being the BSBI's own successor to the Preston atlas and machine readable.
  At most 3 sequences per species and 10 species per genus, each trimmed to the barcode region
  with vsearch against a Barcode UK template.

**The third source theirs had, and this one does not.** Their Methods also name the IRIS BG
database of the Garden's own living collection, which covers the horticultural species planted in
the beds. That catalogue carries no licence statement, so it is not a source here and nothing in
this repository is derived from it. The library is therefore buildable from public records alone,
and the cost of that is measurable rather than assumed: `bin/reference_gap.py` reaches **257 of the
286** marker-level labels the authors report, against 282 when the garden's collection is used as a
target list. The 29 it cannot reach are ornamentals of the planted beds, among them *Begonia*,
*Camellia*, *Chrysanthemum*, *Aster* and *Cercidiphyllum japonicum*.

**Known gap.** Everything outside the two lists, the depth within a species that a vouchered
collection gives, and the planted ornamentals above. `bin/reference_gap.py` measures it after the
fact against the authors' own results: 257 of their 286 marker-level labels, and 81 of the 83 taxa in
their floral surveys, which were held out of the target lists.
**Exome analogy.** The annotation database. A variant can only be called "known" if the database
has it; a read can only be named if the library has its species.

**Nothing from their results builds the library.** The target lists are Barcode UK and the Plant
Atlas 2020; the authors' consensus names are read only afterwards, by `reference_gap.py`, to measure
what the library can reach. Building the library from the names it is then checked against would
make the agreement circular.

## 7. BLAST and assign  (`BLAST`, `ASSIGN`)
**Their command.** `blastn -max_target_seqs 20`, then a spreadsheet curated by hand.
**Ours.** Same BLAST, then the rules written in the paper (Methods 2.6) as code: among the hits
sharing the top bit score, one species -> species; one genus -> genus; one family -> family (or
"A/B/C spp." for up to three genera, as in their tables); several families -> discard; top bit score
below the 1st percentile -> discard.
**Their actual code, which is more than the paper says.** `vsearch-blast-summary.py` sorts the hits
of a cluster by bit score, keeps every hit sharing the top score, and reduces them to unique
(family, genus, species) triples. One triple is a species call and one shared genus is a genus call,
exactly as we do. Below that they have a rule the paper does not mention: if any single genus holds
**60% or more** of the top-score entries, the cluster is called to that genus and marked with a `%`;
only if no genus reaches 60% does it fall back to the shared family, and to `Various` if the families
disagree too. We do not implement it. Ours calls `A/B spp.` where theirs would call the 2-of-3 genus,
and ours is the more conservative of the two in that case.

**Not reproducible literally.** Their botanist also removed plants that do not grow at the site.
That judgement is not recorded. `assets/published/` holds the floral survey if you want to emulate it.
**Exome analogy.** Annotation plus manual variant review, with the review rules made explicit.

## 8. Tables  (`TABLES`)
**What.** Sum reads per taxon per sample. This matrix is the dataset every ecological analysis
starts from, the counterpart of your genotype matrix.
**Check.** `*_assigned` vs `*_clustered` tells you how much is lost to the reference gap.

## 9. Compare with the publication  (`COMPARE`)
**What.** Against the Dryad workbook (`rbcL_consensus_ID_read_counts`, `ITS2_consensus_ID_read_counts`):
reads per sample, taxa above 1% shared, Spearman correlation of relative read abundance per sample.
**Name translation.** Their sheets use *consensus* names (after combining the two markers, e.g.
"Cirsium/Centaurea/Hypochaeris spp."), ours single-marker names (e.g. "Centaurea nigra"). Each of our
names is translated to their label: identical label if it exists; else the label sharing the genus
(their own mapping table, sheet `consensus_identifications`); else the nearest rank above it in the
plant's lineage that they used (tribe, family). `*.name_translation.tsv` lists every translation.
**Metrics.** Bray-Curtis similarity of relative read abundance per insect (1 = identical
composition), Spearman correlation over the taxa present, and whether the dominant plant agrees.

**Tie-break check.** When one of our names could map to several of their labels, the translation picks
the label most common in their data, a choice that favours agreement. Reversing it (least common label)
or excluding ambiguous names leaves the median unchanged on the 40 test insects: rbcL 0.981 in all
three cases (3.8% of rbcL reads carry an ambiguous name), ITS2 0.973 (no ambiguous names). It only
moves insects that already disagree (worst rbcL insect 0.094 / 0.018 / 0.015).

## What the replica is for: calibrating the assumptions

Some parameters are given in neither the paper nor the deposited script, the two differ in places,
and some steps were a judgement made in the lab rather than a written setting. Every one of those
gaps forces a choice here, and a choice made without a target is a guess.

The replica is what gives those guesses a target. Their published tables are a known answer, so each
forced choice can be made both ways and scored against them, and the settings carried into the
reanalysis are ones whose cost against their tables is known rather than assumed. That is what this
arm is for; getting close to the paper is how it earns that use, not a result in itself.

What it calibrated, and what each choice turned out to be worth:

| the ambiguity | the choice | what it costs, measured against their tables |
|---|---|---|
| Crop 3 of the 6 spacer bases, or all 6? | 3, as hard-coded in their script | 6 would put retained reads up to 20% (rbcL) and 42% (ITS2) above ours, where ours sit a median 1.6% and 2.1% above theirs |
| `--cluster_size --id 1.0` by alignment, or pool identical sequences? | pool, with their command kept behind `--cluster_mode align` | +0.5% rbcL and +1.1% ITS2 reads assigned, largest per-insect Bray-Curtis change 0.003, and hours against seconds |
| ITS2 minimum length 350 (paper) or 300 (script)? | 300 | 350 removes most merged ITS2 reads and cannot leave the 11.2 M they report |
| Exact primer matching, or allow mismatches? | allow 2 | their exact rule would keep 99.69% of rbcL and 99.97% of ITS2 reads anyway, so this one is nearly free in both directions |
| Their 2016-2020 tool versions, or current ones? | theirs, pinned | no difference at all; the ITS2 tables came out byte-identical |
| `NexteraPE-PE.fa`, as their script hard-codes, or TruSeq? | TruSeq, theirs shipped beside it | 2 surviving pairs in 74,280 on one sample, and the two outputs differ in two reads |
| Leave the primers on the BLAST queries, as their pipeline does? | trim them | against a mixed library a GenBank reference gains ~30 aligned bases for free; against their primer-free library it changes nothing |
| Which of their consensus labels does one of our names map to? | the most common, which favours agreement | reversing it or dropping ambiguous names leaves the median unchanged, rbcL 0.981 all three ways |
| Barcode UK alone, or Barcode UK plus GenBank fills? | fills | garden beds rbcL 0.952 with fills against 0.821 without, ITS2 0.935 against 0.890 |

**What it cannot calibrate, and this is the boundary that matters.** Three things have no published
target to be scored against, so for them the repository states a judgement and shows its working
instead of claiming a measurement:

- **Their library.** It was never deposited. The arm bounds the term, the row above being how, but a
  bound is not a calibration.
- **The botanist's filter.** Removing plants that do not grow at the site is unrecorded judgement, and
  nothing here reproduces it. It is the main reason the taxon *names* part while the compositions
  match.
- **Every parameter that exists only in the reanalysis**: the identity floors, the mumu thresholds,
  the DADA2 settings. There is no published answer for a method the paper did not use. The floors come
  from the published discrimination of the two markers, which is a literature value and not a
  calibration run here; see `reanalysis.md`.

## Control runs on 40 test insects (18 grassland, 20 garden beds), 21 Sep 2026

Two of the rows above were measured on this subset rather than the full study. Each run changes one
thing against the faithful configuration (their tool versions, rebuilt library, primer-trimmed
queries, pooled dereplication).

| Change | Reads assigned | Largest per-insect change in Bray-Curtis similarity |
|---|---|---|
| Their exact clustering (`--cluster_mode align`) | +0.5% rbcL, +1.1% ITS2 | 0.003 |
| Current tool versions (Trimmomatic 0.39, vsearch 2.28.1, BLAST+ 2.15.0) | none | 0.000 (ITS2 tables byte-identical) |

Their clustering took 2 h 30 min for ITS2 on these 40 samples (1.1 million sequences); it scales
quadratically, so the full study would take days to weeks. Pooled dereplication is kept, with the
measured difference stated above.

Faithful configuration, median similarity to the published tables:

| | Grassland (18) | Garden beds (20) |
|---|---|---|
| rbcL | 0.996 | 0.952 (Barcode UK only: 0.821) |
| ITS2 | 0.996 | 0.935 (Barcode UK only: 0.890) |

## Primer trimming of BLAST queries
Merged reads keep both primers, as in their pipeline. Barcode UK barcodes had their primers removed;
many GenBank records are full-length genes that include the primer site, so any GenBank reference
gained ~30 aligned bases over a perfect Barcode UK match (a 99%-identical *Aster* outscored a
100%-identical *Centaurea*). Primer bases come from the PCR, not the plant, so `trim_primers.py`
removes them from the queries before BLAST. Against a primer-free library, as their native library
was, this changes nothing.

## The random spacer and singletons
`HEADCROP:3` leaves three of the six random bases in front of the forward primer, so one biological
sequence appears as up to 64 variants. Their study-wide singleton filter therefore removes more reads
than errors alone would explain. We keep this behaviour because it is theirs.

## Software versions
The pipeline pins the authors' own versions: Trimmomatic 0.33, FLASH 1.2.11, vsearch 2.3.2,
BLAST+ 2.2.31, as set in `nextflow.config`. Current releases (0.39, 2.28.1, 2.15.0) were run once as
a control and changed nothing, to the point that the ITS2 tables came out byte-identical; none of the
options used changed meaning between the two sets.

## Data quirks
- `GTF3_919_2` was sequenced in two runs. The published table has it once. Our samplesheet keeps both
  runs, naming the second `GTF3_919_2_<run accession>`; only the first is compared.
- During bulk download the archive intermittently answers with a "403 Forbidden" page instead of the
  file. The MD5 check catches it and the task is retried.

- For runs SRR15911821 (`HBC1_719_2`) and SRR15911818 (`HBC1_719_3`) the European archive served a directory listing instead of the
  read-2 file in October 2026, while listing the correct size. The fault was transient and the
  archive now serves both files, so the recipe in the README runs as written. The run was fetched from NCBI with `fasterq-dump`
  (63,766 and 74,280 pairs, as listed; read 1 identical to the European copy) and supplied through a local
  samplesheet (`--input`, with `file://` URLs and the new checksums).

## Not yet implemented
- The botanist's filter. The authors assessed "whether those plants were present within the study
  site and wider landscape" (p.1487) and removed what was not. Nothing here does that.

The consensus step, the blanks, the statistics and the DADA2 arm were all on this list and are now
built; they live in [`reanalysis.md`](reanalysis.md) and their outputs in `04_consensus`,
`05_controls` and `07_ecology`.
