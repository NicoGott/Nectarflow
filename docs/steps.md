# Every step of both pipelines, and where its parameters live

All parameters are in [`nextflow.config`](../nextflow.config), in three labelled blocks: shared,
replica, reanalysis. The only numbers not there are the classification thresholds of the reanalysis,
which are command-line flags and are written into `*.classify_params.tsv` on every run, so a result
table can never be read under cut-offs it was not produced with.

## replica ([`replica.nf`](../replica.nf))

| # | step | tool | key parameters |
|---|---|---|---|
| 1 | TRIM | Trimmomatic 0.33 | `trim_quality`, `adapters`, their `HEADCROP:3` |
| 2 | MERGE | FLASH 1.2.11 | overlap 10-450 |
| 3 | SPLIT_BY_PRIMER | [`split_by_primer.py`](../bin/split_by_primer.py) | `markers` (primer fragments, `min_len` 450 / 300), `primer_mismatches` |
| 4 | DEREP | vsearch 2.3.2 | identical reads within a sample |
| 5 | CLUSTER | vsearch 2.3.2 | `cluster_mode`, 100% identity, study-wide singletons dropped |
| 6 | REFERENCE | [`build_reference.py`](../bin/build_reference.py), BLAST+ | `ref_manifest`, `ref_sequences`. The library is rebuilt from a pinned manifest, searching nothing; [`resolve_reference.py`](../bin/resolve_reference.py) writes that manifest, by hand |
| 7 | BLAST + ASSIGN | BLAST+ 2.2.31, [`assign_taxonomy.py`](../bin/assign_taxonomy.py) | top bit score wins, 1st-percentile cut-off, `blast_chunk` |
| 8 | TABLES | [`build_tables.py`](../bin/build_tables.py) | reads per taxon per sample |
| 9 | COMPARE | [`compare_published.py`](../bin/compare_published.py) | `published`, `published_sheets`, `--rank` |

Where the original cannot be reproduced literally it is replaced by something explicit: the
unpublished reference library by Barcode UK plus the taxa the published results contain and it lacks,
the manual curation by the rules stated in the paper. Every such deviation is justified with a
measurement in [`replica.md`](replica.md), which also gives the original command for each stage and
how to check the output against the paper.

## reanalysis ([`reanalysis.nf`](../reanalysis.nf))

| # | step | tool | key parameters | what changes |
|---|---|---|---|---|
| 1 | CUTADAPT | cutadapt 4.9 | `primers` (full sequences), `primer_mismatches`, `min_len_after_trim` | the 6 random bases **and** both primers come off, and the markers are split, on the raw pairs |
| 2 | DADA2_DENOISE, DADA2_CHIMERA, CURATE | DADA2 1.30.0, vsearch 2.29.1, mumu 1.1.4 | `dada2`: `trunc_len`, `max_ee`, `min_overlap`; `lulu` | exact sequence variants, error model per sequencing run, chimeras removed, no singleton filter; then the extra rDNA copies of one plant are merged back onto it (LULU) |
| 3 | BLAST_ASV, ASSIGN_ASV, CONSENSUS | BLAST+ 2.16, [`assign_taxonomy_asv.py`](../bin/assign_taxonomy_asv.py), [`consensus_markers.py`](../bin/consensus_markers.py) | `max_target_seqs`, identity floors per rank | a name is only as fine as the identity supports; the two markers are combined by stated rules |
| 4 | CONTROLS | [`controls_report.py`](../bin/controls_report.py) | `controls` | the blanks are reported; nothing is subtracted |
| 5 | COMPARE_PUBLISHED_ARM2, COMPARE_ARMS | [`compare_published.py`](../bin/compare_published.py), [`compare_arms.py`](../bin/compare_arms.py) | `replica_results`, `--rank label\|family` | both pipelines against the paper and against each other |
| 6 | ECOLOGY | [`ecology.R`](../bin/ecology.R) (mvabund, vegan) | `conda_ecology` | one analysis script over all three tables, so the pipelines are compared on their conclusions |

Each change, the evidence behind it and what it did to the numbers: [`reanalysis.md`](reanalysis.md).

## Profiles

`conda`, `micromamba`, `docker`, `singularity`; `dev` (20 samples); and `grasspollen`, which points
the replica at a different study with the same chemistry (Brennan et al. 2019, airborne grass pollen,
PRJNA475478, ITS3R instead of UniPlantR) as a portability check.
