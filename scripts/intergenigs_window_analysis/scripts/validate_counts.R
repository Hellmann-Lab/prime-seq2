#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript validate_counts.R <window_counts.txt> <input.bed>")

# Read files
window_counts <- fread(args[1])
reads <- fread(args[2], col.names = c("chr", "start", "end", "read_id", "umi"))

# Filter reads for chromosomes starting with "chr"
reads <- reads[grep("^chr", chr)]

# Calculate totals
total_reads <- nrow(reads)
total_umis <- length(unique(reads$umi[reads$umi != "NA"]))
window_total_reads <- sum(window_counts$read_count)
window_total_umis <- sum(window_counts$unique_umi_count)

# Print validation results
cat("\nValidation Results:\n")
cat("==================\n")
cat(sprintf("Total reads in input: %d\n", total_reads))
cat(sprintf("Total reads in windows: %d\n", window_total_reads))
cat(sprintf("Total unique UMIs in input: %d\n", total_umis))
cat(sprintf("Total unique UMIs in windows: %d\n", window_total_umis))
cat(sprintf("Number of windows: %d\n", nrow(window_counts)))
cat("\n")

# Check read counts
if (total_reads != window_total_reads | total_umis != window_total_umis) {
    cat("ERROR: Counts don't match!\n")
    cat(sprintf("Difference: %d reads\n", abs(total_reads - window_total_reads)))
    cat(sprintf("Difference: %d UMIs\n", abs(total_umis - window_total_umis)))
    quit(status = 1)
}


cat("Validation passed: Read counts match and UMI counts are as expected.\n") 