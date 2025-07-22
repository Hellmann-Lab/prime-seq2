# FP 02/05/2025
# Script to analyze length bias in gene expression data
# Compares UMI and read counts across different transcript length bins
# for different pool sizes (80ng, 320ng, 920ng)

library(tidyverse)
library(rtracklayer)
library(purrr)
library(GenomicFeatures)
library(here)

# Set working directory using here package
here::i_am("scripts/intergenic_analysis/length_bias/expression_length_bins.R")

pool_colors <- c("80ng" = "#C4C082",
                 "320ng" = "#7D7A3B",
                 "920ng" = "#454321"
)


# Load GTF annotation file
gtf <- readGFF(
  "/data/share/htp/Felix_genotyping/zUMIs_tests/own_genomes/mus_musculus/gencode.vM34.primary_assembly.annotation.gtf"
)

# Load GenomicFeatures for transcript length calculation
txdb <- makeTxDbFromGFF(
  "/data/share/htp/Felix_genotyping/zUMIs_tests/own_genomes/mus_musculus/gencode.vM34.primary_assembly.annotation.gtf"
)

# Calculate total exonic length for each transcript
exon_lengths <- exonsBy(txdb, by="tx", use.names=TRUE)
transcript_lengths <- sum(width(reduce(exon_lengths)))

# Select best transcript for each gene based on multiple criteria:
# 1. Presence of tag (CCDS, etc.)
# 2. Transcript support level (lower is better)
# 3. Length -> decide:(longer is better) / (shorter is better)
length_prio <- "long" # long or short

length_per_gene <- gtf %>% 
  filter(type == "transcript") %>% 
  dplyr::select(type, tag, gene_id, transcript_id, transcript_support_level) %>%
  left_join(
    as_tibble(transcript_lengths, rownames="transcript_id"), 
    by="transcript_id"
  ) %>%
  group_by(gene_id) %>%
  mutate(
    has_tag = !is.na(tag) & tag != "NA",
    n_transcripts = n()
  ) %>%
  arrange(
    gene_id,
    desc(has_tag),  # First by presence of tag
    transcript_support_level,  # Then by transcript support level (lower is better)
    if(length_prio == "short") value else -value  # Finally by length (negative for descending)
  ) %>%
  slice_head(n=1) %>%  # Take the first row after arranging
  ungroup() %>%
  dplyr::rename("length"="value") %>%
  dplyr::select(gene_id, length)

# Define input amounts to analyze
names <- c("80ng", "320ng", "920ng")

# Function to process each dataset and count type
process_dataset <- function(name, count_type) {
  # Load count matrix and join with gene lengths
  file_path <- here(
    "data/FC2024_08_01_poolsize/03b_zUMIs_downsampled",
    name,
    "zUMIs_output/expression",
    paste0("poolsize_",
    name,
    ".dgecounts.rds")
  )
  
  counts <- as_tibble(
    as.matrix(readRDS(file_path)[[count_type]][["inex"]][["all"]]), 
    rownames="gene_id"
  ) %>%
    left_join(length_per_gene, by="gene_id")
  
  # Normalize counts by total reads/UMIs per sample
  counts_norm <- counts %>%
    mutate(across(-c(gene_id, length), ~.x/sum(.x)))
  
  # Create three different binning schemes for transcript lengths
  # Very Coarse bins: <1.6, >1.6
  length_bins_only2 <- counts_norm %>%
    mutate(
      length_bin = cut(
        length, 
        breaks = c(0, 1600, #5000, 
                   Inf),
        labels = c("<1.6 kb", #"2-5 kb", 
                   ">1.6 kb")
      )
    ) %>%
    group_by(length_bin) %>%
    summarise(across(-c(gene_id, length), sum)) %>%
    arrange(length_bin) %>%
    pivot_longer(
      cols=-c(length_bin), 
      names_to = "BC", 
      values_to = "proportion"
    ) %>%
    mutate(dataset = name, binning = "only2", count_type = count_type)
  
  # Coarse bins: 0-2 kb, 2-5 kb, >5 kb
  length_bins_coarse <- counts_norm %>%
    mutate(
      length_bin = cut(
        length, 
        breaks = c(0, 2000, 5000, 
                   Inf),
        labels = c("<2 kb", "2-5 kb", 
                   ">5 kb")
      )
    ) %>%
    group_by(length_bin) %>%
    summarise(across(-c(gene_id, length), sum)) %>%
    arrange(length_bin) %>%
    pivot_longer(
      cols=-c(length_bin), 
      names_to = "BC", 
      values_to = "proportion"
    ) %>%
    mutate(dataset = name, binning = "coarse", count_type = count_type)
  
  # Fine bins: 0-1 kb, 1-2 kb, 2-3 kb, 3-4 kb, 4-5 kb, 5-7.5 kb, >7.5 kb
  length_bins_fine <- counts_norm %>%
    mutate(
      length_bin = cut(
        length, 
        breaks = c(0, 1000, 2000, 3000, 4000, 5000, 7500, Inf),
        labels = c("0-1 kb", "1-2 kb", "2-3 kb", "3-4 kb", "4-5 kb", "5-7.5 kb", ">7.5 kb")
      )
    ) %>%
    group_by(length_bin) %>%
    summarise(across(-c(gene_id, length), sum)) %>%
    arrange(length_bin) %>%
    pivot_longer(
      cols=-c(length_bin), 
      names_to = "BC", 
      values_to = "proportion"
    ) %>%
    mutate(dataset = name, binning = "fine", count_type = count_type)
  
  return(bind_rows(length_bins_coarse, length_bins_fine, length_bins_only2))
}

# Process all combinations of datasets and count types
all_length_bins <- map_dfr(names, function(name) {
  bind_rows(
    process_dataset(name, "umicount"),
    process_dataset(name, "readcount")
  )
}) %>%
  mutate(
    dataset = factor(dataset, levels=c("80ng", "320ng", "920ng")),
    count_type = factor(count_type, levels=c("umicount", "readcount")),
    binning = factor(binning, levels=c("only2", "coarse", "fine"))
  )

write.table(all_length_bins, 
            here("scripts/intergenic_analysis/length_bias/", paste0("counts_per_bin.txt"))
            , sep = "\t", row.names = FALSE, quote = FALSE)

# Create visualization with both count types
p1 <- all_length_bins %>%
  ggplot(aes(y=proportion, x=length_bin, color=dataset)) +
  geom_boxplot() +
  facet_grid(count_type~binning, scales="free", space = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    y="Fraction of counts per sample", 
    x="Transcript length", 
    title=paste("Transcript choice: CCDS, support level,", if(length_prio == "long") "longest" else "shortest")
  )+
  scale_color_manual(values = pool_colors)


# Create visualization with only UMI counts
p2 <- all_length_bins %>%
  filter(count_type == "umicount") %>%
  ggplot(aes(y=proportion, x=length_bin, color=dataset)) +
  geom_boxplot() +
  facet_grid(.~binning, scales="free", space = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    y="Fraction of UMIs per sample", 
    x="Transcript length", 
    title=paste("Transcript choice: CCDS, support level,", if(length_prio == "long") "longest" else "shortest")
  )+
  scale_color_manual(values = pool_colors)


# Create visualization with only UMI counts
p3 <- all_length_bins %>%
  filter(count_type == "umicount") %>%
  filter(binning == "only2") %>%
  ggplot(aes(y=proportion, x=length_bin, color=dataset)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    y="Fraction of UMIs per sample", 
    x="Transcript length", 
    title=paste("Transcript choice: CCDS, support level,", if(length_prio == "long") "longest" else "shortest")
  )+
  scale_color_manual(values = pool_colors)


# Save plots as PDFs
ggsave(
  here("scripts/intergenic_analysis/length_bias", 
       paste0("length_bias_all_counts_", 
              if(length_prio == "long") "longest" else "shortest", 
              ".pdf")), 
  p1, 
  width = 10, 
  height = 8
)

ggsave(
  here("scripts/intergenic_analysis/length_bias", 
       paste0("length_bias_umi_counts_", 
              if(length_prio == "long") "longest" else "shortest", 
              ".pdf")), 
  p2, 
  width = 10, 
  height = 4
)

ggsave(
  here("scripts/intergenic_analysis/length_bias", 
       paste0("length_bias_umi_counts_only2_", 
              if(length_prio == "long") "longest" else "shortest", 
              ".pdf")), 
  p3, 
  width = 5, 
  height = 4
)

