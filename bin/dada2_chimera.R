#!/usr/bin/env Rscript
# Combine the per-run tables of one marker, remove chimeras across the whole study, and write
# the ASV table (samples x ASVs), the ASV sequences, and what each step cost in reads.
if (any(commandArgs(trailingOnly = TRUE) %in% c("--help", "-h"))) {
  cat("usage: dada2_chimera.R --marker=M --threads=N SEQTAB.rds... TRACK.tsv...\n\nCombine one marker's per-run tables, remove chimeras across the whole study, write the ASV table and sequences.\n"); quit(status = 0)
}
suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
kv   <- strsplit(sub("^--", "", args[grepl("^--", args)]), "=")
opt  <- setNames(lapply(kv, function(x) paste(x[-1], collapse = "=")), sapply(kv, `[`, 1))
files <- args[!grepl("^--", args)]
marker  <- opt$marker
threads <- as.integer(opt$threads)

# sorted: Nextflow stages the per-run tables in groupTuple order, which is not fixed, and
# mergeSequenceTables keeps the order it is given, so unsorted input would number the ASVs
# differently on two runs over identical data
tabs <- lapply(sort(files[grepl("seqtab\\.rds$", files)]), readRDS)
seqtab <- if (length(tabs) == 1) tabs[[1]] else mergeSequenceTables(tables = tabs)
nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = threads, verbose = TRUE)

asv <- sprintf("%s_ASV%04d", marker, seq_len(ncol(nochim)))
writeLines(as.vector(rbind(paste0(">", asv), colnames(nochim))), sprintf("%s.asvs.fasta", marker))

out <- cbind(sample = rownames(nochim), as.data.frame(nochim))
colnames(out) <- c("sample", asv)
write.table(out, sprintf("%s.asv_table.tsv", marker), sep = "\t", quote = FALSE, row.names = FALSE)

track <- do.call(rbind, lapply(files[grepl("dada2_track\\.tsv$", files)],
                               read.delim, stringsAsFactors = FALSE))
track$nonchimeric <- rowSums(nochim)[track$sample]
track$nonchimeric[is.na(track$nonchimeric)] <- 0
write.table(track[order(track$sample), ], sprintf("%s.dada2_track.tsv", marker),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("%s: %d ASVs before chimera removal, %d after (%.1f%% of reads kept)\n",
            marker, ncol(seqtab), ncol(nochim), 100 * sum(nochim) / sum(seqtab)))
