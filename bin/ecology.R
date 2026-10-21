#!/usr/bin/env Rscript
# Phase 6: the ecological analyses of Lowe et al. 2022 -- the authors' manyglm and a PERMANOVA
# beside it -- over any taxon x sample read-count TSV, so that arm 1, arm 2
# and the published table go through one script. See docs/reanalysis.md, phase 6.
# This script never opens the published workbook: the metadata is parsed out of the sample names
# and the one thing they cannot carry, the insect, comes from --metadata (bin/ecology_metadata.py).
#
# PARSING. A sample code is <letters><digit>_<month><yy>_<insect>, e.g. GBC1_518_01 = site code GBC,
# transect replicate 1, May 2018, insect 1. Checked against their Pollen_sample_metadata sheet: the
# month and year read off the name agree with their Month and Year columns for all 369 insects, and
# the season follows from the month (Mar-May Spring, Jun-Aug Summer, Sep-Nov Autumn; only Mar-Sep
# occur). The letter prefix fixes the habitat, the site and their Transect number (SITES below);
# that is exact for 368 of 369 insects -- HBG1_318_01 is listed under Nature Reserve although every
# other HBG sample is Botanic Garden. A --metadata file always wins over the parsed value and the
# number of disagreements is reported, so that single row is not silently overwritten either way.
# The insect group is NOT in the name and cannot be derived from it: bumblebees and hoverflies were
# collected side by side at the same transect on the same day, which is why --metadata is required.
#
# Sample codes are normalised exactly as compare_published.norm does (GBC1_518_01 -> GBC1_518_1),
# so the same metadata file joins to either arm's tables. Samples whose name does not parse -- the
# four blanks and the second sequencing run of GTF3_919_2 -- are dropped and listed.
#
# THE MODEL. The paper (p.1488) gives season, year and the insects as the predictor variables, and
# compares a four-group model against a two-order one by AIC, carrying the order version forward.
# That selected model is the default here. --formula runs another: "~ season + year + insect_group"
# is the four-group version, and habitat can be added, though adding it is not free -- on every
# table here it takes a deviance 14-220 times larger than any other term, which is what a term that
# nearly separates the data looks like, since plenty of taxa occur in the garden beds and nowhere else.

USAGE <- "
usage: ecology.R --table=<TSV> --metadata=<TSV> --out=<directory> [options]

  --table=FILE     taxon x sample counts, first column 'taxon'; '#' lines are ignored and '__'
                   rows (the read-accounting rows of arm 2) are dropped and reported
  --metadata=FILE  TSV with a 'sample' column; its columns override the parsed ones, and it is
                   the only source of insect_group -- see PARSING in the header
  --out=DIR        output directory, created if needed
  --label=TEXT     name for this table in the summary (default: the file name)
  --formula=F      the manyglm and PERMANOVA model (default '~ season + year + insect_class', the
                   model the authors selected by AIC); see THE MODEL in the header
  --perm=N         permutations for the PERMANOVA (default 999)

Options are given as --name=value.
Everything else is a constant below and is written out to run_info.tsv.

Outputs: manyglm_community.tsv, permanova.tsv, samples_used.tsv, taxa_used.tsv, run_info.tsv,
summary.txt. The ordination is bin/ordination.R, which does all three tables at once.
"

DEFAULTS <- list(table = "", metadata = "", out = "", label = "", formula = "",
                 perm = "999")

FORMULA      <- "~ season + year + insect_class"   # the model the authors selected; --formula overrides
MIN_SAMPLES  <- 5L                # a taxon enters the model if it is present in this many samples
MIN_READS    <- 1L                # log(0) is not an offset; an empty sample is out
SEED         <- 20260923L

# letter prefix -> habitat, site, and their Transect number. Read off Pollen_sample_metadata; the
# prefix determines all three for 368 of the 369 insects (the exception is named in the header).
SITES <- data.frame(
  site_code = c("GBC", "GBG", "GTF", "HBC", "HBG", "SCR", "SW", "WBC"),
  habitat   = c("Grassland", "Grassland", "Horticultural",
                rep("Broadleaved woodland and hedgerows", 5)),
  site      = c("Nature Reserve", "Botanic Garden", "Botanic Garden", "Nature Reserve",
                "Botanic Garden", "Botanic Garden", "Botanic Garden", "Nature Reserve"),
  transect  = 1:8, stringsAsFactors = FALSE)

SEASONS <- c("Winter", "Winter", "Spring", "Spring", "Spring", "Summer",
             "Summer", "Summer", "Autumn", "Autumn", "Autumn", "Winter")   # by month number
SEASON_ORDER <- c("Spring", "Summer", "Autumn", "Winter")
BEES <- c("Bumblebee", "Honeybee", "Solitary bee")

# options

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || any(args %in% c("-h", "--help"))) { cat(USAGE); quit(save = "no") }
if (!all(grepl("^--[a-z_]+=", args)))
  stop("options are given as --name=value; got: ",
       paste(args[!grepl("^--[a-z_]+=", args)], collapse = " "), call. = FALSE)
kv  <- strsplit(sub("^--", "", args), "=")
opt <- setNames(lapply(kv, function(x) paste(x[-1], collapse = "=")), sapply(kv, `[`, 1))
if (length(setdiff(names(opt), names(DEFAULTS))))
  stop("unknown option(s): ", paste(setdiff(names(opt), names(DEFAULTS)), collapse = ", "),
       call. = FALSE)
opt <- modifyList(DEFAULTS, opt)
opt$perm    <- as.integer(opt$perm)
if (opt$table == "" || opt$metadata == "" || opt$out == "")
  stop("--table, --metadata and --out are required\n", USAGE, call. = FALSE)
if (!file.exists(opt$table)) stop("no such table: ", opt$table, call. = FALSE)
if (opt$label == "") opt$label <- basename(opt$table)
if (opt$formula == "") opt$formula <- FORMULA

# Both packages up front: a run that dies twenty minutes in because vegan is missing wastes the
# manyglm resampling. mvabund is the one that has to come from CRAN and needs a C++ compiler --
# checked 23 Sep 2026, neither conda-forge nor bioconda has it, unlike vegan and dada2.
absent <- setdiff(c("mvabund", "vegan"),
                  Filter(function(x) requireNamespace(x, quietly = TRUE), c("mvabund", "vegan")))
if (length(absent))
  stop("missing R package(s): ", paste(absent, collapse = ", "),
       "\n  conda install -c conda-forge r-vegan",
       "\n  Rscript -e 'install.packages(\"mvabund\", repos=\"https://cloud.r-project.org\")'",
       call. = FALSE)
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
outfile <- function(name) file.path(opt$out, name)
write_tsv <- function(x, name, comments = character(0)) {
  con <- file(outfile(name), "w")
  on.exit(close(con))
  num <- vapply(x, is.double, logical(1))           # 17 significant digits of a p value is noise
  x[num] <- lapply(x[num], signif, digits = 6)
  if (length(comments)) writeLines(paste0("# ", comments), con)
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

# the table

norm_sample <- function(x) sub("_0*([0-9]+)$", "_\\1", trimws(as.character(x)))

read_counts <- function(path) {
  lines <- readLines(path, warn = FALSE)
  df <- read.delim(text = lines[!grepl("^#", lines)], check.names = FALSE,
                   stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (ncol(df) < 2 || nrow(df) < 1) stop("table has no data: ", path, call. = FALSE)
  m <- as.matrix(df[, -1, drop = FALSE])
  rownames(m) <- as.character(df[[1]])
  suppressWarnings(storage.mode(m) <- "double")
  if (anyNA(m)) stop(sum(is.na(m)), " non-numeric cell(s) in ", path, call. = FALSE)
  t(m)                                                  # in: taxa x samples; out: samples x taxa
}

counts <- read_counts(opt$table)
rownames(counts) <- norm_sample(rownames(counts))
if (anyDuplicated(rownames(counts)))
  stop("duplicate sample name(s) after normalisation: ",
       paste(unique(rownames(counts)[duplicated(rownames(counts))]), collapse = ", "), call. = FALSE)

# arm 2's *.taxon_table.full.tsv carries read-accounting rows (__unassigned_no_hit and friends):
# reads, not plants, so they are counted and dropped. Per sample and not only in total, because the
# share of an insect's reads that never reached a taxon differs between the arms by design. A table
# without __ rows (the consensus tables, the published sheet) carries no such record at all, and
# the summary then says so rather than leaving it to be assumed.
acc <- grepl("^__", colnames(counts))
acc_reads <- sum(counts[, acc, drop = FALSE])
acc_by_sample <- rowSums(counts[, acc, drop = FALSE])
counts <- counts[, !acc, drop = FALSE]
taxa_all <- colnames(counts)
reads_total <- rowSums(counts)                    # sequencing effort, before any taxon filtering

# metadata

parse_samples <- function(samples) {
  m  <- regmatches(samples, regexec("^([A-Z]{2,4})([0-9]+)_([0-9]{1,2})([0-9]{2})_([0-9]+)$",
                                    samples))
  ok <- lengths(m) == 6
  part <- function(i) vapply(seq_along(m), function(j) if (ok[j]) m[[j]][i] else NA_character_,
                             character(1))
  month <- as.integer(part(4))
  d <- data.frame(sample = samples, parses = ok, site_code = part(2),
                  site_replicate = as.integer(part(3)), month = month,
                  year = ifelse(ok, 2000L + as.integer(part(5)), NA_integer_),
                  season = ifelse(is.na(month), NA_character_, SEASONS[month]),
                  stringsAsFactors = FALSE)
  merge(d, SITES, by = "site_code", all.x = TRUE,
        sort = FALSE)[, c(setdiff(names(d), "site_code"), names(SITES))]
}

md <- parse_samples(rownames(counts))
rownames(md) <- md$sample
md <- md[rownames(counts), ]
md$insect_group <- NA_character_

disagree <- character(0)
ext <- read.delim(opt$metadata, check.names = FALSE, stringsAsFactors = FALSE,
                  quote = "", comment.char = "#")
if (!"sample" %in% names(ext)) stop("--metadata has no 'sample' column", call. = FALSE)
ext$sample <- norm_sample(ext$sample)
ext <- ext[!duplicated(ext$sample), ]
rownames(ext) <- ext$sample
hit <- intersect(rownames(md), rownames(ext))
for (col in setdiff(names(ext), "sample")) {
  v <- ext[hit, col]
  if (col %in% names(md)) {
    # both sides claim to know this column: a mismatch means the naming convention moved and one
    # of the two derivations is stale, which is worth a line in the summary rather than silence
    old <- md[hit, col]
    bad <- which(!is.na(old) & !is.na(v) & trimws(as.character(old)) != trimws(as.character(v)))
    if (length(bad))
      disagree <- c(disagree, sprintf("%s: %d sample(s), e.g. %s (parsed '%s', file '%s')",
                                      col, length(bad), hit[bad[1]], old[bad[1]], v[bad[1]]))
    new <- ifelse(is.na(v) | !nzchar(as.character(v)), as.character(old), as.character(v))
    md[hit, col] <- if (is.numeric(md[[col]])) as.numeric(new) else new
  } else {
    md[[col]] <- NA
    md[hit, col] <- v
  }
}
md$metadata_matched <- rownames(md) %in% hit
md$insect_group[!nzchar(as.character(md$insect_group))] <- NA
# the paper's contrast is bees against hoverflies; the four collected groups collapse onto it
md$insect_class <- ifelse(is.na(md$insect_group), NA_character_,
                          ifelse(md$insect_group == "Hoverfly", "Hoverfly",
                                 ifelse(md$insect_group %in% BEES, "Bee", NA_character_)))
# a group that is neither Hoverfly nor one of BEES falls out of every bee-vs-hoverfly test; count it
# here so that it is a reported exclusion and not an empty cell someone has to notice
unclassed_groups <- sort(unique(md$insect_group[!is.na(md$insect_group) & is.na(md$insect_class)]))
n_unclassed <- sum(!is.na(md$insect_group) & is.na(md$insect_class))
md$total_reads <- reads_total[rownames(md)]
md$accounting_reads <- as.integer(acc_by_sample[rownames(md)])
md$n_taxa <- rowSums(counts > 0)
# max.col breaks a tie by column order, i.e. alphabetically, which is an arbitrary answer to
# "which plant dominated this load" and is the unit the chi-square counts. Record where it happened.

# filtering

md$reason <- ifelse(!md$parses, "sample name is not an insect code",
             ifelse(md$total_reads < MIN_READS,
                    sprintf("fewer than %d reads", MIN_READS), ""))
md$used <- md$reason == ""
keep <- rownames(md)[md$used]
if (length(keep) < 3) stop("only ", length(keep), " usable sample(s); nothing to analyse",
                           call. = FALSE)

prevalence <- colSums(counts[keep, , drop = FALSE] > 0)
taxa_keep <- names(prevalence)[prevalence >= MIN_SAMPLES]
taxa_used <- data.frame(taxon = taxa_all,
                        samples_present = as.integer(prevalence[taxa_all]),
                        reads = as.integer(colSums(counts[keep, , drop = FALSE])[taxa_all]),
                        used = taxa_all %in% taxa_keep, stringsAsFactors = FALSE)
taxa_used <- taxa_used[order(-taxa_used$reads, taxa_used$taxon), ]
if (length(taxa_keep) < 2)
  stop("only ", length(taxa_keep), " taxon/taxa present in >= ", MIN_SAMPLES, " samples",
       call. = FALSE)

Y   <- counts[keep, taxa_keep, drop = FALSE]
mdk <- md[keep, ]
off <- log(mdk$total_reads)                       # sequencing effort of the whole sample, not of Y
rra <- sweep(counts[keep, , drop = FALSE], 1, pmax(mdk$total_reads, 1), "/")

# factors, with the baseline levels fixed so a rerun labels the coefficients the same way
lev <- function(x, order = NULL) {
  x <- as.character(x)
  factor(x, levels = if (is.null(order)) sort(unique(x[!is.na(x)]))
                     else c(intersect(order, x), setdiff(sort(unique(x[!is.na(x)])), order)))
}
mdk$season       <- lev(mdk$season, SEASON_ORDER)
mdk$year         <- lev(mdk$year)
mdk$habitat      <- lev(mdk$habitat)
mdk$site         <- lev(mdk$site)
mdk$site_code    <- lev(mdk$site_code)
mdk$insect_group <- lev(mdk$insect_group)
mdk$insect_class <- lev(mdk$insect_class, c("Bee", "Hoverfly"))

usable_terms <- function(form) {
  wanted <- attr(terms(as.formula(form)), "term.labels")
  ok <- vapply(wanted, function(t) {
    v <- all.vars(as.formula(paste("~", t)))
    all(v %in% names(mdk)) && all(vapply(v, function(x)
      sum(!is.na(mdk[[x]])) == nrow(mdk) && length(unique(mdk[[x]][!is.na(mdk[[x]])])) > 1,
      logical(1)))
  }, logical(1))
  list(keep = wanted[ok], dropped = wanted[!ok])
}
terms_used <- usable_terms(opt$formula)

# A study variable that is a function of a technical batch cannot be told apart from it: every
# sample in one sequencing run then shares the same level of the term, so no permutation and no
# resample can say which of the two produced a difference. Nothing in the p value shows this, so
# it is checked here and reported in the summary. BATCH is the list of columns that are facts about
# the sequencing rather than about the insect.
BATCH <- c("seq_run")
confounds <- unlist(lapply(intersect(BATCH, names(mdk)), function(b) {
  if (length(unique(mdk[[b]][!is.na(mdk[[b]])])) < 2) return(NULL)
  unlist(lapply(terms_used$keep, function(t) {
    tab <- table(mdk[[b]], mdk[[t]]) > 0
    if (!all(rowSums(tab) == 1)) return(NULL)
    map <- vapply(rownames(tab), function(r) colnames(tab)[tab[r, ]][1], character(1))
    sprintf("%s is a function of %s (%s), so the two are the same column here and the %s term cannot be read as %s",
            t, b, paste(sprintf("%s -> %s", names(map), map), collapse = "; "), t, t)
  }))
}))

# manyglm

# mvabund and vegan both rename their output columns between versions, so pick them by prefix
grab <- function(d, pat) { i <- grep(pat, colnames(d))[1]
                           if (is.na(i)) rep(NA_real_, nrow(d)) else as.numeric(d[[i]]) }

Ymv <- mvabund::mvabund(Y)
f <- as.formula(paste("Ymv ~", paste(terms_used$keep, collapse = " + ")))
set.seed(SEED)
# offset is a per-sample vector; mvabund applies the same one to every taxon, which is what
# "total reads as offset" means: the model is about the share of the load, not its depth
fit <- mvabund::manyglm(f, family = "negative.binomial", data = mdk, offset = off)
set.seed(SEED)
# nBoot = 0 asks for the deviances and no resampled p. Insects netted at the same transect are
# not independent, and mvabund can only honour that with case resampling, which does not finish on
# a table this wide. The PERMANOVA below permutes within transect and carries the inference; this
# fit contributes the authors' own statistic, which is what it is here for.
an <- anova(fit, nBoot = 0, test = "LR", p.uni = "none", show.time = "none")
tb <- as.data.frame(an$table)
manyglm_community <- data.frame(
  term = rownames(tb), residual_df = grab(tb, "^Res\\.Df"), df = grab(tb, "^Df"),
  deviance = grab(tb, "^Dev"), test = "LR", stringsAsFactors = FALSE)
write_tsv(manyglm_community, "manyglm_community.tsv",
          c("mvabund::manyglm negative binomial, offset = log(total reads)",
            sprintf("formula: Y ~ %s", paste(terms_used$keep, collapse = " + ")),
            "deviances are SEQUENTIAL, in formula order: each term is adjusted only for the ones before it",
            "no p value: the significance of each term is in permanova.tsv, which permutes within transect",
            sprintf("%d samples x %d taxa", nrow(Y), ncol(Y))))

# PERMANOVA

# Bray-Curtis on relative read abundance: read counts per insect vary by orders of magnitude and
# are not a measure of how much pollen was carried, so composition is all that can be compared.
d <- vegan::vegdist(rra, method = "bray")
set.seed(SEED)
# Insects netted at the same transect are not independent, so the permutations are restricted to
# within transect: a term is then tested against reshuffles that keep each transect's insects
# together. A term that does not vary inside a transect cannot be tested this way and is dropped
# here rather than returned as a degenerate p (habitat is the one that does this).
perm_terms <- Filter(function(t) {
  v <- all.vars(as.formula(paste("~", t)))
  all(vapply(v, function(x) any(tapply(mdk[[x]], mdk$transect,
                                       function(z) length(unique(z))) > 1), logical(1)))
}, terms_used$keep)
perm_dropped <- setdiff(terms_used$keep, perm_terms)
if (!length(perm_terms)) stop("no term varies within transect; the PERMANOVA cannot be blocked")

set.seed(SEED)
adf <- as.data.frame(vegan::adonis2(
  as.formula(paste("d ~", paste(perm_terms, collapse = " + "))),
  data = mdk, by = "margin",
  permutations = permute::how(nperm = opt$perm, blocks = factor(mdk$transect))))
permanova <- data.frame(term = rownames(adf), df = grab(adf, "^Df"),
                        sum_of_squares = grab(adf, "^SumOfSqs"), r2 = grab(adf, "^R2"),
                        f = grab(adf, "^F$"), p_value = grab(adf, "^Pr"),
                        stringsAsFactors = FALSE)
write_tsv(permanova, "permanova.tsv",
          c(sprintf("vegan::adonis2 on the same Bray-Curtis distances, by='margin', %d permutations",
                    opt$perm),
            "each term is adjusted for the others, so season and year are not read in sequence",
            "permutations are restricted to within transect, so insects netted at the same place are",
            "not shuffled across places; this removes site-level pseudoreplication and nothing else",
            if (length(perm_dropped))
              sprintf("dropped, constant within transect and so untestable under blocking: %s",
                      paste(perm_dropped, collapse = ", ")) else NULL,
            "NOT one of the tests named in the paper's methods (docs/reanalysis.md, phase 6);",
            "reported because it puts season and year on one scale (R2) next to the manyglm test"))

# what was used

write_tsv(md[order(rownames(md)),
             c("sample", "parses", "site_code", "site_replicate", "transect", "habitat", "site",
               "month", "year", "season", "insect_group", "insect_class", "metadata_matched",
               "total_reads", "accounting_reads", "n_taxa", "used", "reason")],
          "samples_used.tsv",
          c("total_reads: reads on named taxa, the denominator of this sample's relative abundances",
            "accounting_reads: reads in __ rows (unassigned etc); 0 for a table that carries no such rows"))
write_tsv(taxa_used, "taxa_used.tsv",
          c(sprintf("kept if present in >= %d of the %d used samples", MIN_SAMPLES, length(keep)),
            sprintf("the %d taxa not kept carry %s of the %s reads of the used samples (%.2f%%); they are in the distances but not in the manyglm",
                    sum(!taxa_used$used), format(sum(taxa_used$reads[!taxa_used$used]), big.mark = ","),
                    format(sum(taxa_used$reads), big.mark = ","),
                    100 * sum(taxa_used$reads[!taxa_used$used]) / max(sum(taxa_used$reads), 1))))
write_tsv(data.frame(
  parameter = c("table", "label", "metadata", "tests", "formula", "terms_used", "terms_dropped",
                "min_samples", "min_reads", "perm", "seed",
                "samples_in_table", "samples_used", "taxa_in_table",
                "taxa_used", "reads_used", "reads_in_accounting_rows",
                "reads_per_used_sample_min", "reads_per_used_sample_median",
                "used_samples_under_1000_reads",
                "insect_groups_outside_bee_hoverfly", "date", "r_version",
                "mvabund_version", "vegan_version"),
  value = c(normalizePath(opt$table), opt$label, normalizePath(opt$metadata),
            "manyglm,permanova", opt$formula,
            paste(terms_used$keep, collapse = "+"),
            if (length(terms_used$dropped)) paste(terms_used$dropped, collapse = "+") else "<none>",
            MIN_SAMPLES, MIN_READS, opt$perm, SEED,
            nrow(counts), length(keep), length(taxa_all),
            length(taxa_keep), sum(counts[keep, taxa_keep]), acc_reads,
            min(mdk$total_reads), median(mdk$total_reads), sum(mdk$total_reads < 1000),
            if (n_unclassed) sprintf("%d (%s)", n_unclassed,
                                     paste(unclassed_groups, collapse = "; ")) else "0",
            format(Sys.Date()), R.version.string,
            as.character(packageVersion("mvabund")), as.character(packageVersion("vegan"))),
  stringsAsFactors = FALSE), "run_info.tsv")

# summary

p_of <- function(tbl, term) if (is.null(tbl) || !term %in% tbl$term) NA_real_ else
                            tbl$p_value[match(term, tbl$term)]
dev_of <- function(tbl, term) if (is.null(tbl) || !term %in% tbl$term) NA_real_ else
                              tbl$deviance[match(term, tbl$term)]
r2_of <- function(term) if (is.null(permanova) || !term %in% permanova$term) NA_real_ else
                        permanova$r2[match(term, permanova$term)]
fmt <- function(x, d = 4) if (is.na(x)) "NA" else formatC(x, format = "g", digits = d)

L <- c(sprintf("Ecological analyses of %s", opt$label),
       sprintf("%s, seed %d", format(Sys.time(), "%Y-%m-%d %H:%M"), SEED),
       "",
       sprintf("%d samples in the table, %d used (%d dropped%s)", nrow(counts), length(keep),
               nrow(counts) - length(keep),
               if (all(md$used)) "" else
                 paste0(": ", paste(unique(md$reason[!md$used]), collapse = "; "))),
       sprintf("%d taxa in the table, %d kept (present in >= %d samples), %s reads",
               length(taxa_all), length(taxa_keep), MIN_SAMPLES,
               format(sum(counts[keep, taxa_keep]), big.mark = ",")),
       if (acc_reads > 0)
         sprintf("%s reads (%.1f%% of the table) were in __accounting rows and are not modelled; per sample in samples_used.tsv, up to %.1f%% of one insect's reads",
                 format(acc_reads, big.mark = ","),
                 100 * acc_reads / max(acc_reads + sum(counts), 1),
                 100 * max(md$accounting_reads / pmax(md$accounting_reads + md$total_reads, 1)))
       else paste("this table carries no __accounting rows, so what fraction of each insect's reads",
                  "never reached a taxon cannot be seen here; read it from phase 3/4 and state it",
                  "beside any comparison of arms, whose assigned fractions differ by design"),
       # a composition computed from a handful of reads is noise with a taxon name on it
       sprintf("reads per used sample: median %s, minimum %s; %d used sample(s) under 1,000 reads%s",
               format(median(mdk$total_reads), big.mark = ","),
               format(min(mdk$total_reads), big.mark = ","), sum(mdk$total_reads < 1000),
               if (sum(mdk$total_reads < 1000))
                 sprintf(" (%s) -- their relative abundances are sampling noise; raise MIN_READS (now %d), at the top of this script, to exclude them",
                         paste(sprintf("%s:%s", rownames(mdk)[mdk$total_reads < 1000],
                                       format(mdk$total_reads[mdk$total_reads < 1000], trim = TRUE,
                                              big.mark = ",")), collapse = ", "), MIN_READS)
               else ""),
       sprintf("metadata: %s; insect group known for %d of %d used samples%s",
               opt$metadata, sum(!is.na(mdk$insect_group)), length(keep),
               if (sum(!is.na(mdk$insect_group) & is.na(mdk$insect_class)))
                 sprintf("; %d of them are neither Hoverfly nor a bee (%s) and fall out of every bee-vs-hoverfly test",
                         sum(!is.na(mdk$insect_group) & is.na(mdk$insect_class)),
                         paste(unclassed_groups, collapse = ", ")) else ""),
       if (length(disagree)) c("parsed metadata disagrees with --metadata:",
                               paste0("  ", disagree)) else NULL,
       if (length(terms_used$dropped)) sprintf("terms dropped from the model (missing or constant): %s",
                                               paste(terms_used$dropped, collapse = ", ")) else NULL,
       "",
       if (length(confounds)) c("WARNING: a term in the formula is confounded with the sequencing batch:",
                                paste0("  - ", confounds),
                                "    the effect it reports may be the batch and not the variable; say so beside any",
                                "    p value for it, and do not report it as an effect of the variable alone.") else NULL,
       "")

# The two quantities a reader comes here for: how season compares with year, and whether the
# insect groups differ. Numbers only -- no verdict is adjudicated here.
L <- c(L, "Season and year",
       sprintf("    manyglm deviance: season %s, year %s",
               fmt(dev_of(manyglm_community, "season"), 6),
               fmt(dev_of(manyglm_community, "year"), 6)),
       if (!is.null(permanova)) sprintf("    PERMANOVA: season R2=%s p=%s, year R2=%s p=%s",
                                        fmt(r2_of("season"), 3), fmt(p_of(permanova, "season")),
                                        fmt(r2_of("year"), 3), fmt(p_of(permanova, "year"))) else NULL,
       "    the deviances are sequential, in formula order, so season absorbs what it shares with",
       "    year; the order-free comparison is the PERMANOVA R2.",
       "")

gterm <- intersect(c("insect_class", "insect_group"), terms_used$keep)
L <- c(L, "Bees against hoverflies",
       sprintf("    manyglm deviance: %s", if (length(gterm))
               sprintf("%s %s", gterm[1], fmt(dev_of(manyglm_community, gterm[1]), 6))
               else "no insect term in the model"),
       if (!is.null(permanova) && length(gterm))
         sprintf("    PERMANOVA: %s R2=%s p=%s", gterm[1], fmt(r2_of(gterm[1]), 3),
                 fmt(p_of(permanova, gterm[1]))) else NULL,
       "    the groups are unequal in size and not evenly spread over habitats and seasons; see the",
       "    group sizes in samples_used.tsv and the terms in the model before reading a difference",
       "    between bees and hoverflies as a difference in preference rather than in where and when",
       "    they were caught.", "")

L <- c(L, "Files: run_info.tsv (every effective parameter), samples_used.tsv, taxa_used.tsv,",
       "manyglm_community.tsv, permanova.tsv.")

writeLines(L, outfile("summary.txt"))
writeLines(L)
