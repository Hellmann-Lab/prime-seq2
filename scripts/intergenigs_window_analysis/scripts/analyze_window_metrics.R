#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))

# Function to print formatted messages
print_step <- function(msg, start_time = NULL) {
    if (!is.null(start_time)) {
        elapsed <- difftime(Sys.time(), start_time, units = "secs")
        message(sprintf("\n%s (%.1f seconds)", msg, as.numeric(elapsed)))
    } else {
        message(sprintf("\n%s", msg))
    }
}

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript analyze_window_metrics.R <input.bam>")

# Read BAM file using samtools
bam_file <- args[1]
print_step("Reading BAM file...")
start_time <- Sys.time()
cmd <- sprintf("samtools view %s | awk '{print $3 \"\t\" $4 \"\t\" $4+length($10) \"\t\" $1 \"\t\" $5 \"\t\" $2}'", bam_file)
reads <- fread(cmd = cmd, col.names = c("chr", "start", "end", "read_id", "mapq", "flag"))
print_step("BAM file read complete", start_time)

# Filter for chromosomes starting with "chr"
print_step("Filtering chromosomes...")
start_time <- Sys.time()
reads <- reads[grep("^chr", chr)]
print_step(sprintf("Filtered to %d reads", nrow(reads)), start_time)

# Create GRanges object for reads
print_step("Creating GRanges object...")
start_time <- Sys.time()
gr <- GRanges(
    seqnames = reads$chr,
    ranges = IRanges(start = reads$start, end = reads$start)
)
print_step("GRanges object created", start_time)

# Create 10kb windows
print_step("Creating genomic windows...")
start_time <- Sys.time()
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
print_step(sprintf("Created %d windows", length(windows)), start_time)

# Find overlaps between reads and windows
print_step("Finding read-window overlaps...")
start_time <- Sys.time()
overlaps <- findOverlaps(windows, gr)
print_step(sprintf("Found %d overlaps", length(overlaps)), start_time)

# Convert to data.table for faster processing
print_step("Processing overlaps...")
start_time <- Sys.time()
overlaps_dt <- data.table(
    window_id = queryHits(overlaps),
    read_id = subjectHits(overlaps)
)

# Add window information
windows_dt <- data.table(
    window_id = 1:length(windows),
    chr = as.character(seqnames(windows)),
    start = start(windows),
    end = end(windows)
)

# Pre-calculate strand information to ensure consistent types
reads[, is_reverse := as.integer(bitwAnd(flag, 16) > 0)]
reads[, is_forward := as.integer(!is_reverse)]

# Process MAPQ and strand information separately to avoid type issues
# First, get MAPQ statistics - convert to numeric explicitly
print_step("Calculating MAPQ statistics...")
results_mapq <- overlaps_dt[, .(
    mean_mapq = as.numeric(mean(reads$mapq[read_id])),
    median_mapq = as.numeric(median(reads$mapq[read_id]))
), by = window_id]

# Then, get strand counts - ensure integer type
print_step("Calculating strand counts...")
results_strand <- overlaps_dt[, .(
    forward_strand_count = as.integer(sum(reads$is_forward[read_id])),
    reverse_strand_count = as.integer(sum(reads$is_reverse[read_id]))
), by = window_id]

# Merge the results
print_step("Merging results...")
results <- merge(results_mapq, results_strand, by = "window_id", all = TRUE)

# Merge with window information
results <- merge(windows_dt, results, by = "window_id", all.x = TRUE)

# Fill NA values with 0 for counts
results[is.na(mean_mapq), `:=`(
    mean_mapq = 0,
    median_mapq = 0,
    forward_strand_count = 0,
    reverse_strand_count = 0
)]

# Calculate strand ratio
results[, strand_ratio := forward_strand_count / (forward_strand_count + reverse_strand_count)]
results[is.nan(strand_ratio), strand_ratio := NA]

# Remove window_id column and reorder columns
results[, window_id := NULL]
setcolorder(results, c("chr", "start", "end", "mean_mapq", "median_mapq", 
                     "forward_strand_count", "reverse_strand_count", "strand_ratio"))
print_step("Overlap processing complete", start_time)

# Print summary statistics
print_step("Summary statistics:")
message(sprintf("Total windows: %d", nrow(results)))
message(sprintf("Windows with reads: %d (%.1f%%)", 
                sum(results$forward_strand_count + results$reverse_strand_count > 0),
                mean(results$forward_strand_count + results$reverse_strand_count > 0) * 100))
message(sprintf("Mean MAPQ: %.2f", mean(results$mean_mapq)))
message(sprintf("Median MAPQ: %.2f", median(results$median_mapq)))
message(sprintf("Mean strand ratio: %.2f", mean(results$strand_ratio, na.rm = TRUE)))

# Save results
print_step("Saving results...")
start_time <- Sys.time()

# Get the base name of the BAM file without extension
bam_basename <- basename(bam_file)
bam_basename <- sub("\\.bam$", "", bam_basename)

# Create output directory if it doesn't exist
output_dir <- "output"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Create output file path in the output directory
output_file <- file.path(output_dir, paste0(bam_basename, "_window_metrics.txt"))

write.table(results, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
print_step(sprintf("Results saved to: %s", output_file), start_time) 