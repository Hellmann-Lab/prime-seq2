# FP 29/12/2025
# Script: count_antisense_intergenics.R
# Purpose: This script recounts intergenic reads using reverse strandedness to identify intergenic reads
#          mapping antisense of genes.

# Load required libraries
library(tibble)  # For data frame manipulation
library(here)    # For consistent file path handling
library(dplyr)   # For data manipulation
library(data.table)  # For efficient data processing
library(Matrix)  # For sparse matrices

# Set working directory using here package
here::i_am("scripts/antisense_reads/Antisense_power/count_antisense_intergenics.R")
setwd(here("scripts/antisense_reads/Antisense_power"))

# Source custom functions from zUMIs
source("/data/share/htp/zUMIs-Prime/zUMIs/runfeatureCountFUN.R")
source("/data/share/htp/zUMIs-Prime/zUMIs/misc/featureCounts.R")

.runFeatureCount<-function(abamfile,RG,saf,strand,type,primaryOnly,cpu,mem,fcounts_clib,multi_overlap_var,fraction_overlap){
  print(paste0("Assigning reads to features (",type,")"))
  #  fc.stat<-Rsubread::featureCounts(files=abamfile,
  fc.stat <- featureCounts(files=abamfile,
                           annot.ext=saf,
                           isGTFAnnotationFile=F,
                           primaryOnly=primaryOnly,
                           allowMultiOverlap=multi_overlap_var,
                           countMultiMappingReads=primaryOnly,
                           nthreads=cpu,
                           reportReads="BAM",
                           strandSpecific=strand,
                           isPairedEnd=T,
                           countChimericFragments=F,
                           fcounts_clib = fcounts_clib,
                           largestOverlap = TRUE,
                           fracOverlap = fraction_overlap,
                           isIntronInput = ifelse(type == "in", 1, 0))#$stat
  # fn<-paste0(abamfile,".featureCounts.bam")
  # nfn<-paste0(abamfile,".",type,".featureCounts.bam")
  # 
  # system(paste0("mv ",fn," ",nfn,".tmp"))
  # 
  # invisible(suppressWarnings(suppressMessages(gc(verbose=F))))
  # return(nfn)
  return(fc.stat)
}

# Define the path to the featureCounts library
fcounts_clib <- "/data/share/htp/zUMIs-Prime/zUMIs/misc/fcountsLib2"

# Read BAM files from the configuration file
bam_files <- read.table(here("scripts/quantify_downstream/bam_files.txt"), 
                       header = TRUE, 
                       sep = "\t", 
                       stringsAsFactors = FALSE)

# Extract project names, BAM paths, and barcode file paths
projects <- bam_files$project
bam_paths <- bam_files$bam_path
barcode_paths <- bam_files$BC_WL_path

# Step 1: Create intergenic BAM files with barcode filtering
# This section filters the original BAM files to extract reads that are unassigned to any features
# (intergenic reads) and contain project-specific barcodes using samtools
for (i in seq_along(projects)) {
  project_name <- projects[i]
  bam_path <- bam_paths[i]
  barcode_file <- barcode_paths[i]
  
  # Extract intergenic reads and filter for project barcodes in one command
  cmd <- paste0(
    "\"/data/home/felix/samtools-1.21/samtools\" view -h -e '[ES] == \"Unassigned_NoFeatures\" && [IS] == \"Unassigned_NoFeatures\"' -x ES -x IS -D BC:\"", barcode_file, "\" ",
    "\"", bam_path, "\" >  \"",
    here("scripts/antisense_reads/Antisense_power"),
    "/",
    project_name,
    "_intergenic.bam\""
  )
  cat("Processing", project_name, ":", cmd, "\n")  # Print command for debugging
  system(cmd)
}

# Step 2: Process intergenic reads
# This section uses featureCounts to analyze the intergenic reads for both exonic and intronic regions
for (i in seq_along(projects)) {
  project_name <- projects[i]
  
  # Define input BAM file path
  abamfile <- paste0(here("scripts/antisense_reads/Antisense_power/"), 
                     project_name, 
                     "_intergenic.bam"
                     )
  
  # Create SAF (Simplified Annotation Format) from GTF file
  saf<-.makeSAF(gtf = "/data/share/htp/Felix_genotyping/zUMIs_tests/own_genomes/mus_musculus/gencode.vM34.primary_assembly.annotation.gtf",
                samtoolsexc = "samtools")

  # Count reads in exonic regions
  fnex<-.runFeatureCount(abamfile,
                         saf=saf$exons,
                         strand=2,  # reverse
                         type="ex",
                         primaryOnly = "yes",
                         cpu = 5,
                         mem = 80,
                         fcounts_clib = fcounts_clib,
                         multi_overlap_var = FALSE,
                         fraction_overlap = 0)
  
  # Define path for intron analysis
  intron_in <- here("scripts/antisense_reads/Antisense_power",
                    paste0(project_name, 
                    "_intergenic.bam.featureCounts.bam"))

  # Count reads in intronic regions
  fnin  <-.runFeatureCount(intron_in,
                           saf=saf$introns,
                           strand=2,  # reverse
                           type="in",
                           primaryOnly = "yes",
                           cpu = 5,
                           mem = 80,
                           fcounts_clib = fcounts_clib,
                           multi_overlap_var = FALSE,
                           fraction_overlap = 0)
}

# Step 3: Extract tags and create count matrices
# Define tag names
CELL_TAG <- "BC"
UMI_TAG <- "UB"
GENE_EXON_TAG <- "GE"
GENE_INTRON_TAG <- "GI"

# Function to extract barcode, gene, and UMI tags from BAM file
extract_tags_from_bam <- function(bam_file) {
  cat("Extracting tags from", bam_file, "\n")
  
  # Extract BC, UB, GE, GI tags using samtools and awk
  # Use inex mode: prefer GE, fall back to GI if GE is empty
  cmd <- paste0(
    "\"/data/home/felix/samtools-1.21/samtools\" view -F 2308 \"", bam_file, "\" | ",
    "awk 'BEGIN{OFS=\"\\t\"} {",
    "bc=\"\"; ub=\"\"; ge=\"\"; gi=\"\"; ",
    "for(i=12;i<=NF;i++){",
    "  if(index($i, \"", CELL_TAG, ":Z:\")==1) bc=substr($i,6);",
    "  else if(index($i, \"", UMI_TAG, ":Z:\")==1) ub=substr($i,6);",
    "  else if(index($i, \"", GENE_EXON_TAG, ":Z:\")==1) ge=substr($i,6);",
    "  else if(index($i, \"", GENE_INTRON_TAG, ":Z:\")==1) gi=substr($i,6);",
    "} ",
    "gene=(ge!=\"\") ? ge : gi; ",
    "if(bc!=\"\" && ub!=\"\" && gene!=\"\") print bc, gene, ub;",
    "}'"
  )
  
  # Read the output into a data.table
  tags_dt <- fread(cmd = cmd, sep = "\t", header = FALSE, 
                   col.names = c("barcode", "gene", "umi"),
                   stringsAsFactors = FALSE)
  
  return(tags_dt)
}

# Function to create read count matrix (genes x barcodes)
create_read_count_matrix <- function(tags_dt) {
  cat("Creating read count matrix...\n")
  
  # Count all reads per (gene, barcode) pair
  # Note: .N and .() are data.table syntax, not global variables
  read_counts <- tags_dt[, .N, by = .(gene, barcode)]  # nolint: object_usage_linter
  
  # Convert to wide format (genes as rows, barcodes as columns)
  read_matrix <- dcast(read_counts, gene ~ barcode, value.var = "N", fill = 0)
  
  # Convert to matrix with genes as rownames
  gene_names <- read_matrix$gene
  read_matrix <- as.matrix(read_matrix[, -1, with = FALSE])
  rownames(read_matrix) <- gene_names
  
  return(read_matrix)
}

# Function to create UMI count matrix (genes x barcodes)
create_umi_count_matrix <- function(tags_dt) {
  cat("Creating UMI count matrix...\n")
  
  # Count unique (gene, UMI) pairs per barcode
  # This is equivalent to counting unique UMIs per (gene, barcode)
  # Note: .() and uniqueN are data.table syntax, not global variables
  umi_counts <- tags_dt[, .(umi_count = uniqueN(umi)), by = .(gene, barcode)]  # nolint: object_usage_linter
  
  # Convert to wide format (genes as rows, barcodes as columns)
  umi_matrix <- dcast(umi_counts, gene ~ barcode, value.var = "umi_count", fill = 0)
  
  # Convert to matrix with genes as rownames
  gene_names <- umi_matrix$gene
  umi_matrix <- as.matrix(umi_matrix[, -1, with = FALSE])
  rownames(umi_matrix) <- gene_names
  
  return(umi_matrix)
}

# Process each project and create count matrices
all_read_matrices <- list()
all_umi_matrices <- list()
summary_counts <- c()

for (i in seq_along(projects)) {
  project_name <- projects[i]
  cat("\n=== Processing project:", project_name, "===\n")
  
  # Path to final featureCounts BAM file
  final_bam <- here("scripts/antisense_reads/Antisense_power", 
                    paste0(project_name, "_intergenic.bam.featureCounts.bam.featureCounts.bam"))
  
  if (!file.exists(final_bam)) {
    cat("WARNING: BAM file not found:", final_bam, "\n")
    next
  }
  
  # Extract tags
  tags_dt <- extract_tags_from_bam(final_bam)
  
  if (nrow(tags_dt) == 0) {
    cat("WARNING: No reads with valid tags found for", project_name, "\n")
    summary_counts <- c(summary_counts, 0)
    next
  }
  
  cat("Found", nrow(tags_dt), "reads with valid tags\n")
  summary_counts <- c(summary_counts, nrow(tags_dt))
  
  # Create count matrices
  read_matrix <- create_read_count_matrix(tags_dt)
  umi_matrix <- create_umi_count_matrix(tags_dt)
  
  # Store matrices
  all_read_matrices[[project_name]] <- read_matrix
  all_umi_matrices[[project_name]] <- umi_matrix
  
  cat("Read matrix dimensions:", dim(read_matrix), "\n")
  cat("UMI matrix dimensions:", dim(umi_matrix), "\n")
}

# Create summary table
df_counts <- tibble(Project = projects, Count = summary_counts)

# Save summary results
write.table(df_counts, "counts_antisense_intergenics_all_projects.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Display summary
cat("\n=== Summary ===\n")
print(df_counts)

# Step 4: Save count matrices in zUMIs-like format
# Create a structure similar to dgecounts.rds
cat("\n=== Saving count matrices ===\n")

for (i in seq_along(projects)) {
  project_name <- projects[i]
  
  if (!project_name %in% names(all_read_matrices)) {
    next
  }
  
  # Create output directory
  output_dir <- here("scripts/antisense_reads/Antisense_power", project_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create zUMIs-like structure
  dgecounts <- list(
    readcount = list(
      inex = list(
        all = all_read_matrices[[project_name]]
      )
    ),
    umicount = list(
      inex = list(
        all = all_umi_matrices[[project_name]]
      )
    )
  )
  
  # Save as RDS file
  output_file <- file.path(output_dir, paste0(project_name, "_antisense.dgecounts.rds"))
  saveRDS(dgecounts, output_file)
  cat("Saved count matrices to:", output_file, "\n")
  
  # Convert to data.frame with gene_id as first column for text output
  read_df <- as.data.frame(all_read_matrices[[project_name]])
  read_df <- cbind(gene_id = rownames(read_df), read_df)
  write.table(read_df, read_output, sep = "\t", row.names = FALSE, quote = FALSE)
  
  umi_df <- as.data.frame(all_umi_matrices[[project_name]])
  umi_df <- cbind(gene_id = rownames(umi_df), umi_df)
  write.table(umi_df, umi_output, sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("Saved text files:\n")
  cat("  Read counts:", read_output, "\n")
  cat("  UMI counts:", umi_output, "\n")
}

# Step 5: Optional cleanup
# Note: BAM files are kept by default as they may be needed for debugging
# Uncomment the following section if you want to remove intermediate BAM files

for (i in seq_along(projects)) {
  project_name <- projects[i]

  # Remove the intergenic BAM file
  intergenic_bam <- here("scripts/antisense_reads/Antisense_power", paste0(project_name, "_intergenic.bam"))
  if (file.exists(intergenic_bam)) {
    file.remove(intergenic_bam)
  }

  # Remove the intermediate featureCounts BAM file (keep the final one)
  fc_bam1 <- here("scripts/antisense_reads/Antisense_power", paste0(project_name, "_intergenic.bam.featureCounts.bam"))
  if (file.exists(fc_bam1)) {
    file.remove(fc_bam1)
  }
  
  fc_bam2 <- here("scripts/antisense_reads/Antisense_power", paste0(project_name, "_intergenic.bam.featureCounts.bam.featureCounts.bam"))
  if (file.exists(fc_bam2)) {
    file.remove(fc_bam2)
  }
}

cat("\n=== Analysis complete ===\n")
cat("Count matrices saved in:", here("scripts/antisense_reads/Antisense_power"), "\n")

