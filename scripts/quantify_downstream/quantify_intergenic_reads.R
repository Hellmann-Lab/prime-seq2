#!/usr/bin/env Rscript

#########################################################
#### Quantify Intergenic Reads Downstream of Genes ####
#########################################################

# This script quantifies intergenic reads that map within 10kb downstream 
# of genes for multiple BAM files, following the approach from the manuscript.

#### Load required libraries ####
#################################

# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_packages <- c("GenomicAlignments", "GenomicRanges", "rtracklayer", "stringr", "dplyr", "ggplot2")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
  library(pkg, character.only = TRUE)
}

#### Configuration ####
#######################

# Path to the GTF annotation file
GTF_FILE <- "/data/share/htp/Felix_genotyping/zUMIs_tests/own_genomes/mus_musculus/gencode.vM34.primary_assembly.annotation.gtf"

# Distance threshold for downstream reads (10kb)
DOWNSTREAM_DISTANCE <- 10000

# Output directory
OUTPUT_DIR <- "/data/share/htp/prime-seq_NextGen/scripts/quantify_downstream/results"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

#### Function to create gene ranges BED file ####
#################################################

create_gene_ranges_bed <- function(gtf_file, output_dir) {
  cat("Creating gene ranges BED file...\n")
  
  # Create temporary files for bash processing
  temp_gtf <- file.path(output_dir, "temp_gene_ranges.gtf")
  temp_gtf2 <- file.path(output_dir, "temp_gene_ranges1.gtf")
  bed_file <- file.path(output_dir, "gene_ranges.bed")
  
  # Extract gene entries from GTF
  system(paste("grep -P '\\tgene\\t'", gtf_file, ">", temp_gtf))
  
  # Add transcript_id column if missing
  system(paste("awk '{ if ($0 ~ \"transcript_id\") print $0; else print $0\" transcript_id \\\"\\\";\"; }'", 
               temp_gtf, ">", temp_gtf2))
  
  # Convert to BED format using bedtools
  system(paste("gtf2bed <", temp_gtf2, ">", bed_file))
  
  # Read and process the BED file in R
  gene_ranges <- read.table(bed_file, sep = "\t", stringsAsFactors = FALSE)
  
  # Extract gene names from the last column
  for (i in 1:nrow(gene_ranges)) {
    gene_info <- gene_ranges[i, 10]
    # Extract gene_name using regex
    gene_name_match <- stringr::str_match(gene_info, "gene_name\\s*(.*?)\\s*;")
    if (!is.na(gene_name_match[1, 2])) {
      gene_ranges[i, 10] <- gene_name_match[1, 2]
    } else {
      gene_ranges[i, 10] <- paste0("gene_", i)
    }
  }
  
  # Save processed gene ranges
  write.table(gene_ranges, bed_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Clean up temporary files
  unlink(c(temp_gtf, temp_gtf2))
  
  cat("Gene ranges BED file created:", bed_file, "\n")
  return(bed_file)
}

#### Function to extract intergenic reads from BAM ####
#######################################################

extract_intergenic_reads <- function(bam_file, output_dir, sample_name, barcode_whitelist_path) {
  cat("Processing", sample_name, "...\n")
  
  # Read barcode whitelist
  cat("  Reading barcode whitelist...\n")
  if (file.exists(barcode_whitelist_path)) {
    barcode_whitelist <- readLines(barcode_whitelist_path)
    cat("  Loaded", length(barcode_whitelist), "barcodes from whitelist\n")
  } else {
    cat("  WARNING: Barcode whitelist not found at", barcode_whitelist_path, "\n")
    cat("  Proceeding without barcode filtering\n")
    barcode_whitelist <- NULL
  }
  
  # Read BAM file with all reads to check ES and IS tags
  cat("  Reading BAM file...\n")
  loaded_data <- readGAlignments(
    bam_file, 
    index = bam_file,
    param = ScanBamParam(
      tag = c("ES", "IS", "BC"),
      what = "flag"
    )
  )
  
  if (length(loaded_data) == 0) {
    cat("  No reads found in", sample_name, "\n")
    return(list(bed_file = NULL, total_intergenic = 0))
  }
  
  # Convert to data frame
  tags_df <- data.frame(loaded_data)
  
  # Print tag distribution for debugging
  cat("  Tag distribution:\n")
  cat("    ES tag values:", table(tags_df$ES, useNA = "ifany"), "\n")
  cat("    IS tag values:", table(tags_df$IS, useNA = "ifany"), "\n")
  
  # Filter for intergenic reads: neither ES nor IS is "Unassigned_NoFeatures"
  cat("  Filtering for intergenic reads...\n")
  intergenic_reads <- tags_df$ES == "Unassigned_NoFeatures" & tags_df$IS == "Unassigned_NoFeatures"
  
  tags_df <- tags_df[intergenic_reads, ]
  
  if (nrow(tags_df) == 0) {
    cat("  No intergenic reads found in", sample_name, "\n")
    return(list(bed_file = NULL, total_intergenic = 0))
  }
  
  cat("  Found", nrow(tags_df), "intergenic reads\n")
  
  # Filter by barcode whitelist if available
  if (!is.null(barcode_whitelist)) {
    cat("  Filtering by barcode whitelist...\n")
    valid_barcodes <- tags_df$BC %in% barcode_whitelist
    tags_df <- tags_df[valid_barcodes, ]
    cat("  Found", nrow(tags_df), "intergenic reads with valid barcodes\n")
  }
  
  # Store total count of intergenic reads with valid barcodes
  total_intergenic <- nrow(tags_df)
  
  # Convert to BED format
  gr_tags <- makeGRangesFromDataFrame(tags_df)
  ga_tags <- as(gr_tags, "GAlignments")
  bed_data <- asBED(ga_tags)
  
  # Save intergenic reads BED file
  bed_file <- file.path(output_dir, paste0(sample_name, "_intergenic_reads.bed"))
  export.bed(bed_data, con = bed_file)
  
  return(list(bed_file = bed_file, total_intergenic = total_intergenic))
}

#### Function to find reads within distance of gene ends ####
#############################################################

find_downstream_reads <- function(intergenic_bed, gene_ranges_bed, output_dir, sample_name, distance = 10000) {
  cat("  Finding reads within", distance, "bp downstream of genes...\n")
  
  # Sort BED files
  sorted_intergenic <- file.path(output_dir, paste0(sample_name, "_intergenic_sorted.bed"))
  sorted_genes <- file.path(output_dir, paste0(sample_name, "_genes_sorted.bed"))
  
  system(paste("sortBed -i", intergenic_bed, ">", sorted_intergenic))
  system(paste("sortBed -i", gene_ranges_bed, ">", sorted_genes))
  
  # Find closest genes using bedtools
  results_file <- file.path(output_dir, paste0(sample_name, "_closest_results.txt"))
  system(paste("bedtools closest -a", sorted_intergenic, "-b", sorted_genes, 
               "-s -D a -fu >", results_file))
  
  # Read results
  if (file.exists(results_file) && file.size(results_file) > 0) {
    results <- read.table(results_file, sep = "\t", stringsAsFactors = FALSE)
    
    # Filter for reads within specified distance downstream (negative values = downstream)
    downstream_reads <- results[results$V23 >= -distance & results$V23 < 0, ]
    
    cat("  Found", nrow(downstream_reads), "reads within", distance, "bp downstream\n")
    
    # Clean up temporary files
    unlink(c(sorted_intergenic, sorted_genes))
    
    return(downstream_reads)
  } else {
    cat("  No results found for", sample_name, "\n")
    unlink(c(sorted_intergenic, sorted_genes))
    return(NULL)
  }
}

#### Function to analyze and summarize results ####
##################################################

analyze_downstream_reads <- function(results, sample_name, output_dir) {
  if (is.null(results) || nrow(results) == 0) {
    return(NULL)
  }
  
  cat("  Analyzing results for", sample_name, "...\n")
  
  # Basic statistics
  total_reads <- nrow(results)
  unique_genes <- length(unique(results$V22))
  
  # Distance distribution
  distances <- results$V23
  mean_distance <- mean(distances)
  median_distance <- median(distances)
  
  # Create distance histogram with mean and median lines
  p <- ggplot(data.frame(distance = distances), aes(x = distance)) +
    geom_histogram(binwidth = 100, fill = "steelblue", alpha = 0.7) +
    # Add mean line
    geom_vline(xintercept = mean_distance, color = "red", linetype = "dashed", size = 1) +
    # Add median line
    geom_vline(xintercept = median_distance, color = "orange", linetype = "dashed", size = 1) +
    # Add labels for mean and median
    annotate("text", x = mean_distance, y = Inf, 
             label = paste("Mean:", round(mean_distance, 1), "bp"), 
             vjust = 2, hjust = -0.1, color = "red", fontface = "bold") +
    annotate("text", x = median_distance, y = Inf, 
             label = paste("Median:", round(median_distance, 1), "bp"), 
             vjust = 4, hjust = -0.1, color = "orange", fontface = "bold") +
    theme_classic() +
    labs(
      title = paste("Intergenic reads downstream of genes -", sample_name),
      x = "Distance from 3' gene end (bp)",
      y = "Number of reads"
    ) +
    scale_x_continuous(limits = c(-10000, 0))
  
  plot_file <- file.path(output_dir, paste0(sample_name, "_distance_histogram.pdf"))
  ggsave(plot_file, p, width = 10, height = 6)
  
  # Gene-level summary
  gene_summary <- results %>%
    group_by(V22) %>%
    summarise(
      read_count = n(),
      mean_distance = mean(V23),
      min_distance = min(V23),
      max_distance = max(V23)
    ) %>%
    arrange(desc(read_count))
  
  # Save gene summary
  summary_file <- file.path(output_dir, paste0(sample_name, "_gene_summary.txt"))
  write.table(gene_summary, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Return summary statistics and distance data for combined plot
  summary_stats <- list(
    sample = sample_name,
    total_reads = total_reads,
    unique_genes = unique_genes,
    mean_distance = mean_distance,
    median_distance = median_distance,
    gene_summary = gene_summary,
    distance_data = data.frame(sample = sample_name, distance = distances)
  )
  
  return(summary_stats)
}

#### Function to create combined histogram plot ####
##################################################

create_combined_histogram <- function(all_results, output_dir) {
  if (length(all_results) == 0) {
    return(NULL)
  }
  
  # Combine all distance data
  all_distance_data <- do.call(rbind, lapply(all_results, function(x) {
    if (!is.null(x$distance_data)) {
      return(x$distance_data)
    } else {
      return(NULL)
    }
  }))
  
  if (is.null(all_distance_data) || nrow(all_distance_data) == 0) {
    return(NULL)
  }
  
  # Calculate mean and median for each sample
  sample_stats <- all_distance_data %>%
    group_by(sample) %>%
    summarise(
      mean_dist = mean(distance),
      median_dist = median(distance)
    )
  
  # Create combined histogram with facets and mean/median lines
  p_combined <- ggplot(all_distance_data, aes(x = distance)) +
    geom_histogram(binwidth = 100, fill = "steelblue", alpha = 0.7) +
    # Add mean lines for each sample
    geom_vline(data = sample_stats, aes(xintercept = mean_dist), 
               color = "red", linetype = "dashed", size = 0.8) +
    # Add median lines for each sample
    geom_vline(data = sample_stats, aes(xintercept = median_dist), 
               color = "orange", linetype = "dashed", size = 0.8) +
    # Add labels for mean and median
    geom_text(data = sample_stats, 
              aes(x = mean_dist, y = Inf, 
                  label = paste("Mean:", round(mean_dist, 1))), 
              vjust = 2, hjust = -0.1, color = "red", fontface = "bold", size = 3) +
    geom_text(data = sample_stats, 
              aes(x = median_dist, y = Inf, 
                  label = paste("Median:", round(median_dist, 1))), 
              vjust = 4, hjust = -0.1, color = "orange", fontface = "bold", size = 3) +
    facet_wrap(~sample, scales = "free_y", ncol = 2) +
    theme_classic() +
    theme(
      strip.background = element_rect(fill = "lightgray"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = "Intergenic reads downstream of genes - All samples",
      x = "Distance from 3' gene end (bp)",
      y = "Number of reads"
    ) +
    scale_x_continuous(limits = c(-10000, 0))
  
  # Save combined plot
  combined_plot_file <- file.path(output_dir, "combined_distance_histograms.pdf")
  ggsave(combined_plot_file, p_combined, width = 12, height = 8)
  
  cat("Combined histogram plot saved:", combined_plot_file, "\n")
  
  return(combined_plot_file)
}

#### Main execution ####
#######################

# Read BAM file list (tab-separated)
bam_files <- read.table("/data/share/htp/prime-seq_NextGen/scripts/quantify_downstream/bam_files.txt", 
                        header = TRUE, 
                        sep = "\t", 
                        stringsAsFactors = FALSE,
                        strip.white = TRUE)

# Create gene ranges BED file
gene_ranges_bed <- create_gene_ranges_bed(GTF_FILE, OUTPUT_DIR)

# Process each BAM file
all_results <- list()

for (i in 1:nrow(bam_files)) {
  sample_name <- bam_files$project[i]
  bam_file <- bam_files$bam_path[i]
  barcode_whitelist_path <- bam_files$BC_WL_path[i]
  
  cat("\n=== Processing", sample_name, "===\n")
  
  # Extract intergenic reads
  intergenic_data <- extract_intergenic_reads(bam_file, OUTPUT_DIR, sample_name, barcode_whitelist_path)
  
  if (!is.null(intergenic_data$bed_file) && intergenic_data$total_intergenic > 0) {
    # Find downstream reads
    downstream_results <- find_downstream_reads(
      intergenic_data$bed_file, gene_ranges_bed, OUTPUT_DIR, sample_name, DOWNSTREAM_DISTANCE
    )
    
    # Analyze results
    analysis_results <- analyze_downstream_reads(downstream_results, sample_name, OUTPUT_DIR)
    
    if (!is.null(analysis_results)) {
      # Add total intergenic count and calculate fraction
      analysis_results$total_intergenic <- intergenic_data$total_intergenic
      analysis_results$downstream_fraction <- analysis_results$total_reads / intergenic_data$total_intergenic
      all_results[[sample_name]] <- analysis_results
    }
  }
}

# Create summary report
cat("\n=== SUMMARY REPORT ===\n")

if (length(all_results) > 0) {
  # Combine all results
  summary_df <- do.call(rbind, lapply(all_results, function(x) {
    data.frame(
      Sample = x$sample,
      Total_Intergenic = x$total_intergenic,
      Downstream_Reads = x$total_reads,
      Downstream_Fraction = round(x$downstream_fraction * 100, 2),
      Unique_Genes = x$unique_genes,
      Mean_Distance = round(x$mean_distance, 1),
      Median_Distance = round(x$median_distance, 1),
      stringsAsFactors = FALSE
    )
  }))
  
  # Save summary
  write.table(summary_df, file.path(OUTPUT_DIR, "summary_report.txt"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Print summary
  print(summary_df)
  
  # Create combined histogram plot
  cat("\nCreating combined histogram plot...\n")
  create_combined_histogram(all_results, OUTPUT_DIR)
  
  # Create comparison plot
  if (nrow(summary_df) > 1) {
    p_compare <- ggplot(summary_df, aes(x = Sample, y = Downstream_Fraction)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "Fraction of intergenic reads within 10kb downstream",
        x = "Sample",
        y = "Fraction (%)"
      ) +
      geom_text(aes(label = paste0(Downstream_Fraction, "%")), 
                vjust = -0.5, size = 3)
    
    ggsave(file.path(OUTPUT_DIR, "comparison_plot.pdf"), p_compare, width = 10, height = 6)
  }
  
} else {
  cat("No results found for any samples.\n")
}

cat("\nAnalysis complete! Results saved in:", OUTPUT_DIR, "\n") 