#!/usr/bin/env Rscript
# Denoise one marker of one sequencing run: quality filter, learn that run's error model,
# infer amplicon sequence variants, merge the pairs. Chimeras are removed later, study-wide.
# learnErrors reads a random subset of the files, so --seed fixes it and ends up in the track file.
if (any(commandArgs(trailingOnly = TRUE) %in% c("--help", "-h"))) {
  cat("usage: dada2_denoise.R --marker=M --run=R --trunc_len=a,b --max_ee=a,b --min_overlap=N --threads=N FASTQ...\n\nDenoise one marker of one sequencing run and merge the pairs. Chimeras are removed later, study-wide.\n"); quit(status = 0)
}
suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
kv   <- strsplit(sub("^--", "", args[grepl("^--", args)]), "=")
opt  <- setNames(lapply(kv, function(x) paste(x[-1], collapse = "=")), sapply(kv, `[`, 1))
files <- args[!grepl("^--", args)]

marker  <- opt$marker
run     <- opt$run
trunc   <- as.integer(strsplit(opt$trunc_len, ",")[[1]])      # 0,0 = no truncation (ITS2)
maxee   <- as.numeric(strsplit(opt$max_ee, ",")[[1]])
minover <- as.integer(opt$min_overlap)
threads <- as.integer(opt$threads)
seed    <- if (is.null(opt$seed)) 20260923L else as.integer(opt$seed)

r1 <- sort(files[grepl("\\.R1\\.fastq\\.gz$", files)])
r2 <- sort(files[grepl("\\.R2\\.fastq\\.gz$", files)])
samples <- sub(paste0("\\.", marker, "\\.R1\\.fastq\\.gz$"), "", basename(r1))
stopifnot(length(r1) == length(r2), length(r1) > 0)

f1 <- file.path("filtered", basename(r1)); f2 <- file.path("filtered", basename(r2))
track <- filterAndTrim(r1, f1, r2, f2, truncLen = trunc, maxN = 0, maxEE = maxee,
                       truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = threads)
rownames(track) <- samples

keep <- file.exists(f1) & file.exists(f2)                     # samples with nothing left
f1 <- f1[keep]; f2 <- f2[keep]; kept_samples <- samples[keep]

set.seed(seed)                                                # which files the error model sees
e1 <- learnErrors(f1, multithread = threads, randomize = TRUE, nbases = 1e8)
e2 <- learnErrors(f2, multithread = threads, randomize = TRUE, nbases = 1e8)
d1 <- dada(f1, err = e1, pool = "pseudo", multithread = threads)
d2 <- dada(f2, err = e2, pool = "pseudo", multithread = threads)
mg <- mergePairs(d1, f1, d2, f2, minOverlap = minover, verbose = FALSE)
# with one sample left both return the bare object, not a list of one; everything below wants a list
if (length(kept_samples) == 1) { d1 <- list(d1); d2 <- list(d2); mg <- list(mg) }

seqtab <- makeSequenceTable(mg)
rownames(seqtab) <- kept_samples
saveRDS(seqtab, sprintf("%s.%s.seqtab.rds", marker, run))

getN <- function(x) if (is.null(x)) 0 else sum(getUniques(x))
denoised <- setNames(rep(0, length(samples)), samples)
merged   <- denoised
denoised[kept_samples] <- sapply(d1, getN)
merged[kept_samples]   <- sapply(mg, getN)
out <- data.frame(sample = samples, marker = marker, seq_run = run, seed = seed,
                  cutadapt_pairs = track[, 1], filtered_pairs = track[, 2],
                  denoised = denoised[samples], merged = merged[samples])
write.table(out, sprintf("%s.%s.dada2_track.tsv", marker, run),
            sep = "\t", quote = FALSE, row.names = FALSE)
