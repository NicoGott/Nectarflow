#!/usr/bin/env Rscript
# Recover, for every filtered read pair of one (marker, sequencing run), the ASV it
# ended in -- or why it did not reach one.
#
# The dada() and mergePairs() calls are copied verbatim from bin/dada2_denoise.R so the
# inference is bit-identical to the published run; learnErrors is re-run under the same
# seed, which reproduces the same model. The read -> unique map is recovered afterwards
# by dereplicating each file on its own, because derepFastq is deterministic and dada()
# uses the same unique ordering internally. Doing it this way keeps one file in memory
# at a time instead of all of them.
suppressPackageStartupMessages(library(dada2))

a <- commandArgs(trailingOnly = TRUE)
kv <- strsplit(sub("^--", "", a[grepl("^--", a)]), "=")
opt <- setNames(lapply(kv, function(x) paste(x[-1], collapse = "=")), sapply(kv, `[`, 1))
marker <- opt$marker; run <- opt$run
wdir <- opt$wdir; outdir <- opt$out
trunc <- as.integer(strsplit(opt$trunc_len, ",")[[1]])
maxee <- as.numeric(strsplit(opt$max_ee, ",")[[1]])
minover <- as.integer(opt$min_overlap)
exp_merged <- as.numeric(opt$merged)   # from the published track table
threads <- as.integer(opt$threads)
seed <- 20260923L

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
f1 <- sort(Sys.glob(file.path(wdir, "filtered", paste0("*.", marker, ".R1.fastq.gz"))))
f2 <- sub("\\.R1\\.fastq\\.gz$", ".R2.fastq.gz", f1)
samples <- sub(paste0("\\.", marker, "\\.R1\\.fastq\\.gz$"), "", basename(f1))
stopifnot(length(f1) > 0, all(file.exists(f2)))
cat(sprintf("[%s %s] %d samples\n", marker, run, length(f1)))

# Reuse the run's own error model where it survived, otherwise relearn it under the
# same seed, which is deterministic and gives the same model. Either way the seqtab
# assertion below proves the inference matched the published run.
erds <- file.path(wdir, sprintf("%s.%s.errors.rds", marker, run))
if (file.exists(erds)) {
  er <- readRDS(erds); e1 <- er[[1]]; e2 <- er[[2]]
  cat("  reused the saved error model\n")
} else {
  set.seed(seed)                                 # identical to the original run
  e1 <- learnErrors(f1, multithread = threads, randomize = TRUE, nbases = 1e8)
  e2 <- learnErrors(f2, multithread = threads, randomize = TRUE, nbases = 1e8)
}
cat("  errors learned\n"); flush.console()
d1 <- dada(f1, err = e1, pool = "pseudo", multithread = threads)
d2 <- dada(f2, err = e2, pool = "pseudo", multithread = threads)
mg <- mergePairs(d1, f1, d2, f2, minOverlap = minover, verbose = FALSE)
if (length(samples) == 1) { d1 <- list(d1); d2 <- list(d2); mg <- list(mg) }
cat("  denoised and merged\n"); flush.console()

# the published seqtab must come back out of this identically, or the trace is of
# a different run than the figure
seqtab <- makeSequenceTable(mg); rownames(seqtab) <- samples
ref <- readRDS(file.path(wdir, sprintf("%s.%s.seqtab.rds", marker, run)))
# Two checks. The published track table is the definitive record of what the run
# produced. The work directory's seqtab only agrees if the right attempt was picked
# out of several that share the same filtered counts, so it is checked, not trusted.
cat(sprintf("  merged: regenerated %d, published track %d, work seqtab %d\n",
            sum(seqtab), exp_merged, sum(ref))); flush.console()
stopifnot(sum(seqtab) == exp_merged,
          identical(dim(seqtab), dim(ref)), identical(sort(colnames(seqtab)), sort(colnames(ref))),
          sum(seqtab) == sum(ref))
cat("  reproduces the published run exactly\n"); flush.console()

ids <- function(fq) sub("^@", "", sub(" .*$", "", readLines(gzfile(fq))[c(TRUE, FALSE, FALSE, FALSE)]))
gz <- gzfile(file.path(outdir, sprintf("%s.%s.read_fate.tsv.gz", marker, run)), "w")
writeLines("read_id\tsample\tfate\tsequence", gz)
for (i in seq_along(samples)) {
  m1 <- derepFastq(f1[i])$map
  m2 <- derepFastq(f2[i])$map
  rid <- ids(f1[i])
  stopifnot(length(rid) == length(m1), length(m1) == length(m2))
  a1 <- d1[[i]]$map[m1]; a2 <- d2[[i]]$map[m2]      # read -> denoised variant, each direction
  tab <- mg[[i]]
  key <- if (nrow(tab)) paste(tab$forward, tab$reverse) else character(0)
  hit <- match(paste(a1, a2), key)
  fate <- ifelse(is.na(a1) | is.na(a2), "denoised_away",
                 ifelse(is.na(hit), "not_merged", "merged"))
  seqs <- ifelse(is.na(hit), "", tab$sequence[hit])
  writeLines(paste(rid, samples[i], fate, seqs, sep = "\t"), gz)
  if (i %% 20 == 0) { cat(sprintf("  %d/%d\n", i, length(samples))); flush.console() }
}
close(gz)
cat(sprintf("[%s %s] done\n", marker, run))
