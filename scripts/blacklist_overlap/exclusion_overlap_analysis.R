#!/usr/bin/env Rscript
# Script: exclusion_overlap_analysis.R
# Purpose: Analyze overlaps between BAM files and genomic exclusion ranges
# Author: 
# Date: 

# 1. Setup and Configuration ----
# Install and load required packages
library("BiocManager")
library("GenomicRanges")
library("AnnotationHub")
library("Rsamtools")
library("ggplot2")
library("tibble")
library("dplyr")
library("tidyr")

# Define colors for conditions
pool_colors <- c("80ng" = "#C4C082",
                 "320ng" = "#7D7A3B",
                 "920ng" = "#454321")

# Define BAM files
bam_files <- c(
    "/data/share/htp/prime-seq_NextGen/data/FC2024_08_01_poolsize/03_zUMIs/80ng/poolsize_80ng.filtered.Aligned.GeneTagged.sorted.bam",
    "/data/share/htp/prime-seq_NextGen/data/FC2024_08_01_poolsize/03_zUMIs/320ng/poolsize_320ng.filtered.Aligned.GeneTagged.sorted.bam",
    "/data/share/htp/prime-seq_NextGen/data/FC2024_08_01_poolsize/03_zUMIs/920ng/poolsize_920ng.filtered.Aligned.GeneTagged.sorted.bam"
)

# 2. Data Loading and Validation ----
# Check if BAM files exist
for (bam_file in bam_files) {
    if (!file.exists(bam_file)) {
        stop(sprintf("BAM file not found: %s", bam_file))
    }
}

# Get exclusion ranges using AnnotationHub
ah <- AnnotationHub()
# Query for mouse GRCm39/mm39 exclusion ranges
query_data <- query(ah, c("excluderanges", "mm39"))
# Get the mouse exclusion ranges
exclusion_ranges <- query_data[["AH107321"]]  # GRCm39/mm39 exclusion ranges
# Sort and keep standard chromosomes
exclusion_ranges <- sort(exclusion_ranges)

# 3. Helper Functions ----
# Function to extract condition name from BAM file path
get_condition_name <- function(bam_file) {
    # Extract the condition (e.g., "80ng") from the path
    condition <- gsub(".*/([0-9]+ng)/.*", "\\1", bam_file)
    return(condition)
}

# Function to classify read type based on ES and IS tags
classify_read_type <- function(ES, IS) {
    # Convert to character and handle NA values
    ES <- as.character(ES)
    IS <- as.character(IS)
    
    # Handle NA values
    if (is.na(ES) || is.na(IS)) {
        return("unknown")
    }
    
    # Classify based on tags
    if (ES != "Unassigned_NoFeatures") {
        return("exonic")
    } else if (IS != "Unassigned_NoFeatures") {
        return("intronic")
    } else {
        return("intergenic")
    }
}

# 4. BAM Analysis Function ----
# Function to analyze BAM file
analyze_bam <- function(bam_file) {
    condition_name <- get_condition_name(bam_file)
    cat(sprintf("Processing: %s\n", condition_name))
    
    # Read BAM file
    tryCatch({
        # Create BamFile object
        bf <- BamFile(bam_file)
        
        # Define what to read from BAM file
        what <- c("rname", "pos", "qwidth")
        tag <- c("ES", "IS")
        
        # Read all alignments with specified tags
        bam <- scanBam(bf, param = ScanBamParam(what = what, tag = tag))[[1]]
        total_reads <- length(bam$pos)
        cat(sprintf("Successfully read %d reads from %s\n", total_reads, condition_name))
        
        # Extract ES and IS tags
        ES <- bam$tag$ES
        IS <- bam$tag$IS
        
        # Create tibble with all necessary information
        temp <- tibble(rname = bam$rname,
                       pos = bam$pos, 
                       width = bam$qwidth,
                       ES = ES,
                       IS = IS) %>%
            filter(!is.na(rname) & !is.na(pos) & !is.na(width))
        
        # Classify read types
        temp <- temp %>%
            rowwise() %>%
            mutate(read_type = classify_read_type(ES, IS)) %>%
            ungroup()
        
        # Create GRanges directly from BAM data
        bam_gr <- GRanges(seqnames = as.character(temp$rname),
                          ranges = IRanges(start = temp$pos, 
                                           width = temp$width))
        
        # Calculate overlaps with exclusion ranges
        overlaps <- findOverlaps(bam_gr, exclusion_ranges)
        
        # Create detailed overlap annotation
        overlap_annotation <- tibble(
            read_id = queryHits(overlaps),
            exclusion_id = subjectHits(overlaps),
            read_chr = as.character(seqnames(bam_gr[queryHits(overlaps)])),
            read_start = start(bam_gr[queryHits(overlaps)]),
            read_end = end(bam_gr[queryHits(overlaps)]),
            exclusion_chr = as.character(seqnames(exclusion_ranges[subjectHits(overlaps)])),
            exclusion_start = start(exclusion_ranges[subjectHits(overlaps)]),
            exclusion_end = end(exclusion_ranges[subjectHits(overlaps)]),
            exclusion_name = (exclusion_ranges[subjectHits(overlaps)])$name,
            exclusion_strand = as.character(strand(exclusion_ranges[subjectHits(overlaps)])),
            read_type = temp$read_type[queryHits(overlaps)]
        )
        
        # Calculate statistics
        valid_total_reads <- length(bam_gr)
        overlapping_reads <- length(unique(queryHits(overlaps)))
        overlap_percentage <- (overlapping_reads / valid_total_reads) * 100
        
        cat(sprintf("Found %d reads (%.2f%%) overlapping exclusion ranges\n", 
                    overlapping_reads, overlap_percentage))
        
        return(list(
            condition = condition_name,
            total_reads = total_reads,
            valid_reads = valid_total_reads,
            overlapping_reads = overlapping_reads,
            overlap_percentage = overlap_percentage,
            filtered_reads = total_reads - valid_total_reads,
            overlap_annotation = overlap_annotation
        ))
    }, error = function(e) {
        cat(sprintf("Error processing %s: %s\n", condition_name, e$message))
        return(NULL)
    })
}

# 5. Main Analysis ----
# Analyze each BAM file
results <- lapply(bam_files, analyze_bam)

# Remove NULL results
results <- results[!sapply(results, is.null)]

# Create summary statistics tibble
summary_stats <- tibble(
    condition = sapply(results, function(x) x$condition),
    total_reads = sapply(results, function(x) x$total_reads),
    valid_reads = sapply(results, function(x) x$valid_reads),
    overlapping_reads = sapply(results, function(x) x$overlapping_reads),
    overlap_percentage = sapply(results, function(x) x$overlap_percentage),
    filtered_reads = sapply(results, function(x) x$filtered_reads)
)

# Create combined overlap annotation tibble
combined_annotations <- bind_rows(
    lapply(results, function(x) {
        x$overlap_annotation %>%
            mutate(condition = x$condition) %>%
            select(condition, everything())
    })
)

# 6. Basic Overlap Analysis Plots ----
# Plot 1: Fraction of excluded reads by exclusion type and read type
p1 <- combined_annotations %>%
    count(condition, exclusion_name, read_type, name="excluded_reads") %>%
    left_join(summary_stats, by="condition") %>%
    mutate(excluded_reads_fract = excluded_reads/valid_reads) %>%
    mutate(condition = factor(condition, levels = c("80ng", "320ng", "920ng"))) %>%
    mutate(read_type = factor(read_type, levels = c("exonic", "intronic", "intergenic", "unknown"))) %>%
    ggplot(aes(y=excluded_reads_fract, x=exclusion_name, fill=condition))+
    geom_col(position = "dodge")+
    coord_flip()+
    facet_wrap(~read_type)+
    scale_fill_manual(values = pool_colors)+
    labs(title = "Fraction of Excluded Reads by Exclusion Type",
         subtitle = "Grouped by Read Type",
         y = "Fraction of Excluded Reads",
         x = "Exclusion Type",
         fill = "Condition")+
    theme_bw()+
    theme(plot.title = element_text(size=14, face="bold"),
          plot.subtitle = element_text(size=12),
          axis.title = element_text(size=12),
          strip.text = element_text(size=12),
          legend.title = element_text(size=12),
          legend.text = element_text(size=10))

# Save plot 1
ggsave("/data/share/htp/prime-seq_NextGen/blacklist_overlap/excluded_reads_by_type.pdf", 
       p1, width = 12, height = 6, device = "pdf")

# Plot 2: Fraction of excluded reads by chromosome, condition, and read type
p2 <- combined_annotations %>%
    count(condition, exclusion_name, read_chr, read_type, name="excluded_reads") %>%
    left_join(summary_stats, by="condition") %>%
    mutate(excluded_reads_fract = excluded_reads/valid_reads) %>%
    mutate(condition = factor(condition, levels = c("80ng", "320ng", "920ng"))) %>%
    mutate(read_chr = factor(read_chr, 
                             levels = c(paste0("chr", 1:19), "chrX", "chrY", "chrM"),
                             ordered = TRUE)) %>%
    mutate(read_type = factor(read_type, levels = c("exonic", "intronic", "intergenic", "unknown"))) %>%
    ggplot(aes(y=excluded_reads_fract, x=exclusion_name, fill=condition))+
    geom_col(position = "dodge")+
    coord_flip()+
    facet_grid(read_type~read_chr, scales = "free_y")+
    scale_fill_manual(values = pool_colors)+
    labs(title = "Fraction of Excluded Reads by Chromosome",
         subtitle = "Grouped by Read Type and Condition",
         y = "Fraction of Excluded Reads",
         x = "Exclusion Type",
         fill = "Condition")+
    theme_bw()+
    theme(plot.title = element_text(size=14, face="bold"),
          plot.subtitle = element_text(size=12),
          axis.title = element_text(size=12),
          strip.text = element_text(size=10),
          legend.title = element_text(size=12),
          legend.text = element_text(size=10))

# Save plot 2
ggsave("/data/share/htp/prime-seq_NextGen/blacklist_overlap/excluded_reads_by_chromosome_and_type.pdf", 
       p2, width = 18, height = 8, device = "pdf")

# 7. 10kb Windows Analysis ----
# Function to create 10kb windows
create_windows <- function(granges) {
    # Get chromosome lengths from the GRanges object
    chr_lengths <- seqlengths(granges)
    # Create windows
    windows <- tileGenome(chr_lengths, tilewidth = 10000, cut.last.tile.in.chrom = TRUE)
    return(windows)
}

# Create windows
genome_windows <- create_windows(exclusion_ranges)

# Convert combined_annotations to GRanges
read_granges <- GRanges(
    seqnames = combined_annotations$read_chr,
    ranges = IRanges(
        start = combined_annotations$read_start,
        end = combined_annotations$read_end
    ),
    condition = combined_annotations$condition,
    read_type = combined_annotations$read_type
)

# Find overlaps between windows and reads
window_read_overlaps <- findOverlaps(genome_windows, read_granges)

# Create a data frame with window information and read overlap counts
window_data <- tibble(
    chr = as.character(seqnames(genome_windows)),
    start = start(genome_windows),
    end = end(genome_windows),
    window_id = paste0(as.character(seqnames(genome_windows)), ":", 
                       start(genome_windows), "-", 
                       end(genome_windows))
) %>%
    mutate(chr = factor(chr, 
                        levels = c(paste0("chr", 1:19), "chrX", "chrY", "chrM"),
                        ordered = TRUE))

# Count reads per window by condition and read type
window_counts <- tibble(
    window_id = window_data$window_id[queryHits(window_read_overlaps)],
    condition = read_granges$condition[subjectHits(window_read_overlaps)],
    read_type = read_granges$read_type[subjectHits(window_read_overlaps)]
) %>%
    count(window_id, condition, read_type, name = "read_count") %>%
    complete(window_id, condition, read_type, fill = list(read_count = 0))

# Join with window data
window_data <- window_data %>%
    left_join(window_counts, by = "window_id")

# 8. Window Analysis Plots ----
# Plot 3: Top 10 windows by mean excluded reads across conditions
p3 <- window_data  %>% 
    filter(!is.na(condition)) %>%
    pivot_wider(names_from = "condition", values_from = "read_count") %>%
    mutate(mean = rowMeans(select(., `320ng`, `80ng`, `920ng`), na.rm = TRUE)) %>%
    arrange(read_type, desc(mean)) %>%
    group_by(read_type) %>%
    mutate(rank = row_number()) %>%
    pivot_longer(cols = c(6:8), names_to = "condition", values_to = "read_count") %>%
    ungroup() %>%
    mutate(condition = factor(condition, levels = c("80ng", "320ng", "920ng"))) %>%
    filter(rank<10) %>%
    ggplot(aes(y=read_count, x=rank, color=condition))+
    geom_line(size=1.2)+
    geom_point(size=2)+
    facet_wrap(~read_type, scales = "free_y")+
    scale_color_manual(values = pool_colors)+
    labs(title = "Top 10 Windows with Highest Mean Excluded Reads",
         subtitle = "Grouped by Read Type",
         x = "Window Rank (1 = highest mean excluded reads)",
         y = "Number of Excluded Reads",
         color = "Condition")+
    theme_bw()+
    theme(plot.title = element_text(size=14, face="bold"),
          plot.subtitle = element_text(size=12),
          axis.title = element_text(size=12),
          strip.text = element_text(size=12),
          legend.title = element_text(size=12),
          legend.text = element_text(size=10))

# normalize to valid reads first
p3 <- window_data  %>% 
    filter(!is.na(condition)) %>%
    pivot_wider(names_from = "condition", values_from = "read_count") %>%
    mutate(mean = rowMeans(select(., `320ng`, `80ng`, `920ng`), na.rm = TRUE)) %>%
    arrange(read_type, desc(mean)) %>%
    group_by(read_type) %>%
    mutate(rank = row_number()) %>%
    pivot_longer(cols = c(6:8), names_to = "condition", values_to = "read_count") %>%
    ungroup() %>%
    mutate(condition = factor(condition, levels = c("80ng", "320ng", "920ng"))) %>%
    filter(rank<10) %>%
    ggplot(aes(y=read_count, x=rank, color=condition))+
    geom_line(size=1.2)+
    geom_point(size=2)+
    facet_wrap(~read_type, scales = "free_y")+
    scale_color_manual(values = pool_colors)+
    labs(title = "Top 10 Windows with Highest Mean Excluded Reads",
         subtitle = "Grouped by Read Type",
         x = "Window Rank (1 = highest mean excluded reads)",
         y = "Number of Excluded Reads",
         color = "Condition")+
    theme_bw()+
    theme(plot.title = element_text(size=14, face="bold"),
          plot.subtitle = element_text(size=12),
          axis.title = element_text(size=12),
          strip.text = element_text(size=12),
          legend.title = element_text(size=12),
          legend.text = element_text(size=10))

# Save plot 3
ggsave("/data/share/htp/prime-seq_NextGen/blacklist_overlap/top10_windows_by_read_type.pdf", 
       p3, width = 12, height = 8, device = "pdf")

# Plot 4: Top 10 windows overall by mean excluded reads
p4 <- window_data  %>% 
    filter(!is.na(condition)) %>%
    pivot_wider(names_from = "condition", values_from = "read_count") %>%
    mutate(mean = rowMeans(select(., `320ng`, `80ng`, `920ng`), na.rm = TRUE)) %>%
    arrange(read_type, desc(mean)) %>%
    group_by(read_type) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    filter(rank<10) %>%
    ggplot(aes(y=mean, x=reorder(window_id, mean), fill=read_type))+
    geom_col()+
    coord_flip()+
    scale_fill_brewer(palette = "Set2")+
    labs(title = "Top 10 Windows with Highest Mean Excluded Reads",
         subtitle = "Across All Read Types",
         y = "Mean Number of Excluded Reads",
         x = "Genomic Window",
         fill = "Read Type")+
    theme_bw()+
    theme(plot.title = element_text(size=14, face="bold"),
          plot.subtitle = element_text(size=12),
          axis.title = element_text(size=12),
          axis.text.y = element_text(size=8),
          legend.title = element_text(size=12),
          legend.text = element_text(size=10))

p4

# Save plot 4
ggsave("/data/share/htp/prime-seq_NextGen/blacklist_overlap/top10_windows_overall.pdf", 
       p4, width = 10, height = 8, device = "pdf")


# facetted
window_data  %>% 
    filter(!is.na(condition)) %>%
    pivot_wider(names_from = "condition", values_from = "read_count") %>%
    mutate(mean = rowMeans(select(., `320ng`, `80ng`, `920ng`), na.rm = TRUE)) %>%
    arrange(read_type, desc(mean)) %>%
    group_by(read_type) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    filter(rank<10) %>%
    pivot_longer(cols=c(6:8), values_to = "excl_reads", names_to = "condition") %>%
    mutate(condition = factor(condition, levels = c("80ng", "320ng", "920ng"))) %>% View()
    ggplot(aes(y=excl_reads, x=reorder(window_id, mean), fill=read_type))+
    geom_col()+
    coord_flip()+
    scale_fill_brewer(palette = "Set2")+
    labs(title = "Top 10 Windows with Highest Mean Excluded Reads",
         subtitle = "Across All Read Types",
         y = "Mean Number of Excluded Reads",
         x = "Genomic Window",
         fill = "Read Type")+
    theme_bw()+
    facet_wrap(~condition)+
    theme(plot.title = element_text(size=14, face="bold"),
          plot.subtitle = element_text(size=12),
          axis.title = element_text(size=12),
          axis.text.y = element_text(size=8),
          legend.title = element_text(size=12),
          legend.text = element_text(size=10))

    