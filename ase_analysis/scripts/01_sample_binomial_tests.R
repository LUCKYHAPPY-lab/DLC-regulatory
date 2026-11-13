#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript 01_sample_binomial_tests.R COUNT_DIR OUTDIR")
files <- list.files(args[1], pattern = "_ASEReadCounter[.]table$", full.names = TRUE)
if (!length(files)) stop("No ASEReadCounter tables found.")
ids <- sub("_.*$", "", basename(files))
if (anyDuplicated(ids)) stop("Duplicate sample IDs extracted from filenames.")
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(files)) {
  x <- read.delim(files[i], stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("refCount", "altCount", "totalCount")
  if (!all(required %in% names(x))) stop("Missing count columns in ", files[i])
  x$ref_ratio <- x$refCount / x$totalCount
  x$alt_ratio <- x$altCount / x$totalCount
  keep <- with(x, refCount >= 3 & altCount >= 3 & totalCount >= 10 &
                 ref_ratio >= 0.01 & alt_ratio >= 0.01)
  x <- x[which(keep), , drop = FALSE]
  x$p <- vapply(seq_len(nrow(x)), function(j) {
    binom.test(x$refCount[j], x$totalCount[j], p = 0.5)$p.value
  }, numeric(1))
  x$FDR <- p.adjust(x$p, method = "BH")
  write.table(x, file.path(args[2], paste0(ids[i], "_ASEReadCounter_binom_test.txt")),
              row.names = FALSE, sep = "\t", quote = FALSE)
}
