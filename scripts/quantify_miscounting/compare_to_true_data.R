# Load required libraries
library(data.table)
library(Rsamtools)
library(dplyr)

# Set your working directory
setwd("/home/felix/prime-seq_NextGen/scripts/quantify_miscounting/STAR_out")

# Read BAM file mapping
bam_map <- fread("../true_bam_files.txt") # adjust path if needed
# bam_map should have columns: project, bam_path

# List all *_counts_detail.txt files
count_files <- list.files(pattern = "*_counts_detail.txt$")

# Function to read counts detail file
read_counts_detail <- function(file) {
  dt <- fread(file)
  # Columns: read, exon, intron, anyIntergenic
  return(dt)
}

# Function to extract ES:Z: and IS:Z: tags from BAM
extract_bam_tags <- function(bam_file) {
  param <- ScanBamParam(tag = c("ES", "IS"), what = "qname")
  bam <- scanBam(bam_file, param = param)[[1]]
  df <- data.frame(
    read = bam$qname,
    ES = bam$tag$ES,
    IS = bam$tag$IS,
    stringsAsFactors = FALSE
  )
  return(df)
}

# Main comparison loop
results <- list()

for (count_file in count_files) {
  # Extract project/sample name from file name
  # e.g. STAR_out/PoP64_PTO_allbest_counts_detail.txt -> PoP64_PTO
  sample_name <- sub("_allbest_counts_detail.txt$", "", basename(count_file))
  sample_name <- sub("_counts_detail.txt$", "", sample_name) # fallback

  # Find BAM path for this sample
  bam_path <- bam_map$bam_path[bam_map$project == sample_name]
  if (length(bam_path) == 0) {
    warning(paste("No BAM file found for", sample_name))
    next
  }

  # Read counts detail
  counts <- read_counts_detail(count_file)

  # Read BAM tags
  bam_tags <- extract_bam_tags(bam_path)

  # Merge by read
  merged <- left_join(counts, bam_tags, by = "read")

  # Summarize: for each combination of (exon, intron, ES, IS)
  summary <- merged %>%
    mutate(
      assignment = case_when((exon == TRUE | intron == TRUE) & anyIntergenic == TRUE ~ "InEx & Intergenic",
                                                  exon == TRUE | intron == TRUE ~ "InEx only",
                                                  anyIntergenic == TRUE ~ "Intergenic only",
                                                  T ~ "other")) %>%
    mutate(true_assignment = case_when(ES == "Assigned3" | IS == "Assigned3" ~ "InEx",
                                       T ~ "Intergenic")) %>%
    # group_by(assignment, true_assignment) %>%
    summarise(n = n(), .by=c(assignment, true_assignment)) %>%
    arrange(desc(n))

  results[[sample_name]] <- summary
  print(paste("Summary for", sample_name))
  print(summary)
}

# Optionally, save results
lapply(names(results), function(n) write.table(results[[n]], paste0(n, "_comparison_summary.txt"), sep = "\t", row.names = FALSE, quote = FALSE))
