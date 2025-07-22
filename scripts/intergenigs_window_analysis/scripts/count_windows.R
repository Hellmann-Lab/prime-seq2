#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript count_windows.R <input.bed>")

# Read BED file
reads <- fread(args[1], col.names = c("chr", "start", "end", "read_id", "umi"))

# Filter for chromosomes starting with "chr"
reads <- reads[grep("^chr", chr)]

# Create GRanges object for reads
gr <- GRanges(
    seqnames = reads$chr,
    ranges = IRanges(start = reads$start, end = reads$start) # just use the starting position
)

# Create 10kb windows
chr_lengths <- reads[, .(max_end = max(end)), by = chr]
chr_lengths <- chr_lengths[order(chr)]

# Create windows for each chromosome separately
windows <- GRanges()
for (i in 1:nrow(chr_lengths)) {
    chr <- chr_lengths$chr[i]
    max_end <- chr_lengths$max_end[i]
    chr_windows <- tileGenome(
        setNames(max_end, chr),
        tilewidth = 10000,
        cut.last.tile.in.chrom = TRUE
    )
    windows <- c(windows, chr_windows)
}

# Count reads and unique UMIs per window
read_counts <- countOverlaps(windows, gr)

# Count unique UMIs per window with fractional counting
overlaps <- findOverlaps(windows, gr)
umi_counts <- integer(length(windows))  # Initialize with zeros

# Get unique UMIs and their window appearances
umi_window_mapping <- data.table(
    window_id = queryHits(overlaps),
    umi = reads$umi[subjectHits(overlaps)]
)

# Calculate fractional counts
umi_window_mapping <- umi_window_mapping[umi != "NA"]  # Remove NA UMIs
umi_window_mapping[, window_count := .N, by = umi]  # Count windows per UMI
umi_window_mapping[, fractional_count := 1/window_count]  # Calculate fractional count

# Sum fractional counts per window
umi_counts <- umi_window_mapping[, .(unique_umi_count = sum(fractional_count)), by = window_id]

# Create results
results <- data.frame(
    chr = seqnames(windows),
    start = start(windows),
    end = end(windows),
    read_count = read_counts,
    unique_umi_count = 0  # Initialize with zeros
)

# Fill in UMI counts for windows that have them
results$unique_umi_count[umi_counts$window_id] <- umi_counts$unique_umi_count

# Save results
output_file <- sub("\\.bed$", "_window_counts.txt", args[1])
write.table(results, output_file, sep = "\t", row.names = FALSE, quote = FALSE) 