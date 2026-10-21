#!/usr/bin/env Rscript
#
# The ordination behind docs/figures/ordination.png, recomputed for the figure rather than taken
# from phase 6.
#
# Phase 6 ordinates every sample it analyses, which is right for the results table. It is wrong for
# a picture: an insect whose whole load is one taxon sits at Bray-Curtis distance 1 from nearly
# every other insect, NMDS cannot place it, and it both flies off the plot and flatters the stress
# of whichever table carries it most. This drops such samples, says how many, and tries harder.
#
#   bin/ordination.R --tables published=a.tsv,replica=b.tsv,reanalysis=c.tsv --out DIR
#
# Each table is taxa x samples, as phase 4 and bin/published_table.py write them. Writes one
# <name>.nmds.tsv of coordinates per table, plus ordination_fit.tsv with the stress and the
# Procrustes agreement between every pair.
suppressMessages(library(vegan))

SEED    <- 20260923L
TRYMAX  <- 400L
MAXDIST <- 0.999      # "shares essentially nothing with anything else"
SHARE   <- 0.95       # ... with at least this share of the other samples

args <- commandArgs(trailingOnly = TRUE)
get <- function(flag, default = NA) {
  i <- grep(paste0("^", flag, "="), args)
  if (length(i)) sub(paste0("^", flag, "="), "", args[i[1]]) else default
}
spec <- get("--tables"); out <- get("--out", ".")
if (is.na(spec)) stop("--tables name=path,name=path,... is required")
dir.create(out, recursive = TRUE, showWarnings = FALSE)   # the ordination is written here
paths <- strsplit(strsplit(spec, ",")[[1]], "=")
names <- vapply(paths, `[`, "", 1); files <- vapply(paths, `[`, "", 2)

set.seed(SEED)
# taxa x samples on disk, samples x taxa for vegan
tabs <- lapply(files, function(f) t(as.matrix(read.delim(f, check.names = FALSE, row.names = 1))))
# phase 4 writes sample codes with a zero-padded tail, the workbook without: match compare_published
for (k in seq_along(tabs)) rownames(tabs[[k]]) <- sub("_0*([0-9]+)$", "_\\1", rownames(tabs[[k]]))
names(tabs) <- names
shared <- Reduce(intersect, lapply(tabs, rownames))

# a sample is dropped if ANY table cannot place it, so all the panels hold the same insects
drop <- character(0)
for (k in names) {
  p <- sweep(tabs[[k]][shared, ], 1, rowSums(tabs[[k]][shared, ]), "/")
  d <- as.matrix(vegdist(p, method = "bray"))
  diag(d) <- NA
  far <- rowMeans(d >= MAXDIST, na.rm = TRUE) >= SHARE
  drop <- union(drop, shared[far])
}
keep <- setdiff(shared, drop)
cat(sprintf("%d samples shared, %d unplaceable and dropped (%s), %d ordinated\n",
            length(shared), length(drop), paste(drop, collapse = ", "), length(keep)))

fits <- list(); co <- list()
for (k in names) {
  p <- sweep(tabs[[k]][keep, ], 1, rowSums(tabs[[k]][keep, ]), "/")
  p <- p[, colSums(p) > 0, drop = FALSE]
  m <- metaMDS(vegdist(p, method = "bray"), k = 2, trymax = TRYMAX, trace = 0, autotransform = FALSE)
  s <- as.data.frame(scores(m, display = "sites"))
  write.table(data.frame(sample = rownames(s), NMDS1 = s$NMDS1, NMDS2 = s$NMDS2),
              file.path(out, paste0(k, ".nmds.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  co[[k]] <- s[keep, 1:2]
  fits[[length(fits) + 1]] <- data.frame(table = k, metric = "stress", value = sprintf("%.5f", m$stress))
  fits[[length(fits) + 1]] <- data.frame(table = k, metric = "converged", value = as.character(isTRUE(m$converged)))
  fits[[length(fits) + 1]] <- data.frame(table = k, metric = "n_taxa", value = as.character(ncol(p)))
  cat(sprintf("%-11s stress %.5f  converged %s  %d taxa\n", k, m$stress, isTRUE(m$converged), ncol(p)))
}
# r is vegan's Procrustes correlation, sqrt(1 - m12): 1 is the same arrangement of the same insects
for (i in seq_along(names)) for (j in seq_along(names)) if (i < j) {
  pr <- protest(co[[names[i]]], co[[names[j]]], permutations = 999)
  r <- sqrt(1 - pr$ss)
  fits[[length(fits) + 1]] <- data.frame(table = paste(names[i], "vs", names[j]),
                                         metric = "procrustes_r", value = sprintf("%.4f", r))
  fits[[length(fits) + 1]] <- data.frame(table = paste(names[i], "vs", names[j]),
                                         metric = "procrustes_p", value = sprintf("%.4f", pr$signif))
  cat(sprintf("%-11s vs %-11s Procrustes r = %.4f  p = %.4f\n", names[i], names[j], r, pr$signif))
}
fits[[length(fits) + 1]] <- data.frame(table = "all", metric = "dropped", value = paste(drop, collapse = " "))
fits[[length(fits) + 1]] <- data.frame(table = "all", metric = "n_samples", value = as.character(length(keep)))
fits[[length(fits) + 1]] <- data.frame(table = "all", metric = "seed", value = as.character(SEED))
fits[[length(fits) + 1]] <- data.frame(table = "all", metric = "trymax", value = as.character(TRYMAX))
write.table(do.call(rbind, fits), file.path(out, "ordination_fit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
