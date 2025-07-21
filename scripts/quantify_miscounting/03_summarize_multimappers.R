#!/usr/bin/env Rscript
library(tidyverse)

OUTDIR <- "/home/felix/prime-seq_NextGen/scripts/quantify_miscounting/STAR_out"
# PROJECT <- "PoP64_PTO"
PROJECT <- Sys.getenv("PROJECT")
if (PROJECT == "") stop("PROJECT environment variable not set!")
PREFIX  <- file.path(OUTDIR, paste0(PROJECT, "_allbest"))

# read in lists
all_reads  <- read_lines(paste0(PREFIX, "_reads_all.txt"))
exon_reads <- read_lines(paste0(PREFIX, "_reads_exon.txt"))
intn_reads <- read_lines(paste0(PREFIX, "_reads_intron.txt"))
anyint_reads <- read_lines(paste0(PREFIX, "_reads_any_intergenic.txt"))

df <- tibble(read = all_reads) %>%
  mutate(
    exon       = read %in% exon_reads,
    intron     = read %in% intn_reads,
    anyIntergenic = read %in% anyint_reads
  )

# Calculate total mapped reads using system command
bam_file <- file.path(OUTDIR, paste0(PROJECT, "_allbest_Aligned.sortedByCoord.out.bam"))
total_mapped_reads <- as.integer(system(
  paste0(
    "samtools view ", bam_file, " | awk '{print $1}' | sort | uniq | wc -l"
  ),
  intern = TRUE
))

summary_df <- df %>%
  count(exon, intron, anyIntergenic) %>%
  arrange(desc(n)) %>%
  mutate(fraction_of_total_mapped_reads = n / total_mapped_reads)

summary_df %>%
  write_tsv(paste0(PREFIX, "_counts.txt"))

df %>%
  write_tsv(paste0(PREFIX, "_counts_detail.txt"))

# Delete all files with the prefix except for the _counts.txt file
files_to_delete <- list.files(
  path = OUTDIR,
  pattern = paste0(basename(PREFIX), "_"),
  full.names = TRUE
)
# files_to_delete <- files_to_delete[!grepl(paste0("^", PREFIX, ".*(_counts\\.txt)?$"), files_to_delete)]
# # Try to remove files, and if not, try to remove as directories
# for (f in files_to_delete) {
#   if (file.info(f)$isdir) {
#     unlink(f, recursive = TRUE)
#   } else {
#     file.remove(f)
#   }
# }
