# FP 08/04/2025
# Script: count_antisense_intergenics.R
# Purpose: This script recounts intergenic reads using reverse strandedness to identify intergenic reads
#          mapping antisense of genes.

# Load required libraries
library(tibble)  # For data frame manipulation
library(here)    # For consistent file path handling

# Set working directory using here package
here::i_am("scripts/intergenic_analysis/count_antisense_intergenics.R")
setwd(here("scripts/intergenic_analysis"))

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

# Define experimental conditions
conditions <- c("80ng", "320ng", "920ng")

# Step 1: Create intergenic BAM files
# This section filters the original BAM files to extract reads that are unassigned to any features
# (intergenic reads) using samtools
for (c in conditions) {
  cmd <- paste0(
    "\"/data/home/felix/samtools-1.21/samtools\" view -h -e '[ES] == \"Unassigned_NoFeatures\" && [IS] == \"Unassigned_NoFeatures\"' -x ES -x IS ",
    "\"", 
    # here("data/FC2024_08_01_poolsize/03_zUMIs/"),
    here("data/FC2024_08_01_poolsize/03b_zUMIs_downsampled/"),
    c,
    "/poolsize_",
    c,
    ".filtered.Aligned.GeneTagged.sorted.bam\" >  \"",
    here("scripts/intergenic_analysis"),
    "/",
    c,
    "_intergenic.bam\""
  )
  cat(cmd, "\n")  # Print command for debugging
  system(cmd)
}

# Step 2: Process intergenic reads
# This section uses featureCounts to analyze the intergenic reads for both exonic and intronic regions
for (c in conditions) {
  # Define input BAM file path
  abamfile <- paste0(here("scripts/intergenic_analysis/"), 
                     c, 
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
  intron_in <- here("scripts",
                    paste0(
                    "intergenic_analysis/",
                     c, 
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

# Step 3: Count and summarize results
# Initialize vector to store counts
counts <- c()

# Count reads with gene information (GI:Z or GE:Z tags) for each condition
for (c in conditions) {
  cmd <- paste0(
    "samtools view \"", here("scripts/intergenic_analysis"), "/", c, "_intergenic.bam.featureCounts.bam.featureCounts.bam\" | awk '/GI:Z:/ || /GE:Z:/' | wc -l"
  )
  cat(cmd, "\n")  # Print command for debugging
  
  result <- system(cmd, intern = TRUE)
  counts <- c(counts, as.numeric(result))
}

# Create a tibble with the results
df_counts <- tibble(Condition = conditions, Count = counts)

# Save results to a tab-separated file
write.table(df_counts, "counts_antisense_intergenics_downsampled.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Display results
print(df_counts)

# Step 4: Clean up temporary files
# Remove all BAM files created during the analysis
for (c in conditions) {
  # Remove the intergenic BAM file
  intergenic_bam <- here("scripts/intergenic_analysis", paste0(c, "_intergenic.bam"))
  if (file.exists(intergenic_bam)) {
    file.remove(intergenic_bam)
  }
  
  # Remove the featureCounts BAM files
  fc_bam1 <- here("scripts/intergenic_analysis", paste0(c, "_intergenic.bam.featureCounts.bam"))
  fc_bam2 <- here("scripts/intergenic_analysis", paste0(c, "_intergenic.bam.featureCounts.bam.featureCounts.bam"))
  
  if (file.exists(fc_bam1)) {
    file.remove(fc_bam1)
  }
  if (file.exists(fc_bam2)) {
    file.remove(fc_bam2)
  }
}

