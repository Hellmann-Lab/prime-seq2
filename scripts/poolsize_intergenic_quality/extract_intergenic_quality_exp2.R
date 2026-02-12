#!/usr/bin/env Rscript
# Extract quality metrics for intergenic reads from experiment 2
# Metrics: MAPQ, aligned length, soft-clipping fraction, mismatch rate

library(Rsamtools)
library(GenomicAlignments)
library(here)
library(dplyr)

here::i_am("scripts/poolsize_intergenic_quality/extract_intergenic_quality_exp2.R")
setwd(here("scripts/poolsize_intergenic_quality"))

# Function to parse CIGAR string and extract aligned length and soft-clipping
parse_cigar <- function(cigar_string) {
  # Parse CIGAR operations
  # M = match/mismatch (counts toward aligned length on reference)
  # I = insertion (doesn't count toward aligned length on reference, but counts in query)
  # D = deletion (counts toward aligned length on reference, not in query)
  # S = soft clipping (doesn't count toward aligned length, but counts in query)
  # H = hard clipping (doesn't count toward aligned length or query)
  # N = skipped region (counts toward aligned length on reference)
  
  if (is.na(cigar_string) || cigar_string == "*") {
    return(list(aligned_length = 0, soft_clip_start = 0, soft_clip_end = 0, total_length = 0, soft_clip_fraction = 0))
  }
  
  # Use GenomicAlignments to get aligned length on reference
  aligned_length <- cigarWidthAlongReferenceSpace(cigar_string)
  
  # Get total query length (read length including soft-clipping)
  total_length <- cigarWidthAlongQuerySpace(cigar_string)
  
  # Parse CIGAR string to extract soft-clipping
  # Split by operation type
  ops <- strsplit(cigar_string, "(?<=[A-Z])", perl = TRUE)[[1]]
  
  soft_clip_start <- 0
  soft_clip_end <- 0
  
  if (length(ops) > 0) {
    # Check first operation
    first_op <- ops[1]
    if (grepl("S$", first_op)) {
      soft_clip_start <- as.numeric(gsub("[A-Z]", "", first_op))
    }
    
    # Check last operation
    last_op <- ops[length(ops)]
    if (grepl("S$", last_op)) {
      soft_clip_end <- as.numeric(gsub("[A-Z]", "", last_op))
    }
  }
  
  total_soft_clip <- soft_clip_start + soft_clip_end
  soft_clip_fraction <- ifelse(total_length > 0, total_soft_clip / total_length, 0)
  
  return(list(
    aligned_length = aligned_length,
    soft_clip_start = soft_clip_start,
    soft_clip_end = soft_clip_end,
    total_soft_clip = total_soft_clip,
    total_length = total_length,
    soft_clip_fraction = soft_clip_fraction
  ))
}

# Function to extract quality metrics from BAM file
extract_intergenic_quality <- function(bam_file, barcode_file, sample_name) {
  cat("Processing sample:", sample_name, "\n")
  cat("  BAM file:", bam_file, "\n")
  
  # Read barcode whitelist if provided
  barcode_whitelist <- NULL
  if (!is.null(barcode_file) && file.exists(barcode_file)) {
    barcode_whitelist <- readLines(barcode_file)
    cat("  Loaded", length(barcode_whitelist), "barcodes from whitelist\n")
  }
  
  # Read BAM file with intergenic filter
  cat("  Reading BAM file...\n")
  param <- ScanBamParam(
    tag = c("ES", "IS", "BC", "nM"),
    what = c("qname", "mapq", "cigar", "qwidth")
  )
  
  bam_data <- scanBam(bam_file, param = param)[[1]]
  
  if (length(bam_data$qname) == 0) {
    cat("  No reads found\n")
    return(NULL)
  }
  
  # Handle missing tags - replace NULL with NA
  n_reads <- length(bam_data$qname)
  
  # Check and handle each tag
  ES_tag <- if (is.null(bam_data$tag$ES) || length(bam_data$tag$ES) == 0) {
    rep(NA_character_, n_reads)
  } else {
    bam_data$tag$ES
  }
  
  IS_tag <- if (is.null(bam_data$tag$IS) || length(bam_data$tag$IS) == 0) {
    rep(NA_character_, n_reads)
  } else {
    bam_data$tag$IS
  }
  
  BC_tag <- if (is.null(bam_data$tag$BC) || length(bam_data$tag$BC) == 0) {
    rep(NA_character_, n_reads)
  } else {
    bam_data$tag$BC
  }
  
  # nM is STAR's tag for number of mismatches (lowercase n)
  nM_tag <- if (is.null(bam_data$tag$nM) || length(bam_data$tag$nM) == 0) {
    rep(NA_integer_, n_reads)
  } else {
    bam_data$tag$nM
  }
  
  # Ensure all vectors have the same length
  if (length(ES_tag) != n_reads || length(IS_tag) != n_reads || 
      length(BC_tag) != n_reads || length(nM_tag) != n_reads) {
    stop("Tag vectors have different lengths. This should not happen.")
  }
  
  # Convert to data frame
  reads_df <- data.frame(
    qname = bam_data$qname,
    mapq = bam_data$mapq,
    cigar = bam_data$cigar,
    qwidth = bam_data$qwidth,
    ES = ES_tag,
    IS = IS_tag,
    BC = BC_tag,
    nM = nM_tag,
    stringsAsFactors = FALSE
  )
  
  cat("  Total reads:", nrow(reads_df), "\n")
  
  # Filter for intergenic reads
  intergenic <- reads_df$ES == "Unassigned_NoFeatures" & 
                reads_df$IS == "Unassigned_NoFeatures"
  reads_df <- reads_df[intergenic, ]
  
  cat("  Intergenic reads:", nrow(reads_df), "\n")
  
  if (nrow(reads_df) == 0) {
    cat("  No intergenic reads found\n")
    return(NULL)
  }
  
  # Filter by barcode whitelist if provided
  if (!is.null(barcode_whitelist)) {
    reads_df <- reads_df[reads_df$BC %in% barcode_whitelist, ]
    cat("  After barcode filtering:", nrow(reads_df), "reads\n")
  }
  
  if (nrow(reads_df) == 0) {
    cat("  No reads after barcode filtering\n")
    return(NULL)
  }
  
  # Parse CIGAR strings
  cat("  Parsing CIGAR strings...\n")
  cigar_info <- lapply(reads_df$cigar, parse_cigar)
  
  reads_df$aligned_length <- sapply(cigar_info, function(x) x$aligned_length)
  reads_df$soft_clip_fraction <- sapply(cigar_info, function(x) x$soft_clip_fraction)
  reads_df$total_length <- sapply(cigar_info, function(x) x$total_length)
  
  # Calculate mismatch rate
  # nM tag is STAR's number of mismatches (lowercase n)
  # Mismatch rate = nM / aligned_length
  # Handle missing nM tags (set to NA)
  reads_df$mismatch_rate <- ifelse(!is.na(reads_df$nM) & reads_df$aligned_length > 0,
                                   reads_df$nM / reads_df$aligned_length,
                                   NA_real_)
  
  # Extract poolsize from sample name (e.g., "80_1_1" -> "80", "920_1_4" -> "920")
  poolsize <- gsub("_.*", "", sample_name)
  reads_df$poolsize <- poolsize
  reads_df$sample <- sample_name
  
  # Select relevant columns
  result <- reads_df %>%
    select(poolsize, sample, mapq, aligned_length, soft_clip_fraction, mismatch_rate)
  
  cat("  Extracted", nrow(result), "intergenic reads\n")
  
  return(result)
}

# Experiment 2: FC2025_06_01_poolsize2
base_dir <- here("data/FC2025_06_01_poolsize2/03_zUMIs")

# Find all sample directories
sample_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
sample_dirs <- sample_dirs[sample_dirs != ""]

cat("Found sample directories:", paste(sample_dirs, collapse = ", "), "\n\n")

all_data <- data.frame()

for (sample_dir in sample_dirs) {
  # Find BAM file in this directory
  bam_pattern <- paste0("poolsize2_", sample_dir, ".filtered.Aligned.GeneTagged.sorted.bam")
  bam_file <- file.path(base_dir, sample_dir, bam_pattern)
  
  if (!file.exists(bam_file)) {
    cat("Warning: BAM file not found:", bam_file, "\n")
    next
  }
  
  # Find corresponding barcode file
  # Try to match poolsize (80 or 920)
  poolsize <- gsub("_.*", "", sample_dir)
  barcode_file <- file.path(base_dir, paste0(poolsize, "_BCs.txt"))
  
  # If that doesn't exist, set to NULL
  if (!file.exists(barcode_file)) {
    cat("  Warning: Barcode file not found:", barcode_file, "\n")
    barcode_file <- NULL
  }
  
  data <- extract_intergenic_quality(bam_file, barcode_file, sample_dir)
  if (!is.null(data)) {
    all_data <- rbind(all_data, data)
  }
}

# Save results
output_file <- "intergenic_quality_exp2.rds"
saveRDS(all_data, output_file)
cat("\nSaved results to:", output_file, "\n")
cat("Total intergenic reads:", nrow(all_data), "\n")

# Print summary statistics
if (nrow(all_data) > 0) {
  summary_stats <- all_data %>%
    group_by(poolsize) %>%
    summarise(
      n_reads = n(),
      mean_mapq = mean(mapq, na.rm = TRUE),
      median_mapq = median(mapq, na.rm = TRUE),
      mean_aligned_length = mean(aligned_length, na.rm = TRUE),
      median_aligned_length = median(aligned_length, na.rm = TRUE),
      mean_soft_clip_fraction = mean(soft_clip_fraction, na.rm = TRUE),
      median_soft_clip_fraction = median(soft_clip_fraction, na.rm = TRUE),
      mean_mismatch_rate = mean(mismatch_rate, na.rm = TRUE),
      median_mismatch_rate = median(mismatch_rate, na.rm = TRUE),
      .groups = "drop"
    )
  
  print(summary_stats)
}

