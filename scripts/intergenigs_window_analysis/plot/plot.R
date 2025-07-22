suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

# Define consistent colors for pool sizes
pool_colors <- c("80ng" = "#C4C082",
                 "320ng" = "#7D7A3B",
                 "920ng" = "#454321"
)

# 1 get data ----
## 1.1 window counts ----
window_files <- list.files("/data/share/htp/prime-seq_NextGen/scripts/intergenics_analysis/output", pattern = "_window_counts.txt$", full.names = TRUE)
all_data <- lapply(window_files, function(f) {
  # Extract condition name from filename
  condition_name <- basename(f)
  condition_name <- gsub("_window_counts.txt$", "", condition_name)
  condition_name <- gsub("_intergenic$", "", condition_name)
  
  condition_name <- gsub("\\.filtered.*", "", condition_name)
  condition_name <- gsub("poolsize_", "", condition_name)
  
  # Read data
  dt <- fread(f)
  # Remove "chr" prefix from chromosome column
  dt[, chr := gsub("chr", "", chr)]
  dt[, condition := condition_name]
  return(dt)
})

# Combine all data
combined_data <- rbindlist(all_data) %>%
  mutate(condition = factor(condition, levels=c("80ng", "320ng", "920ng")))


## 1.2 reads for normalizations ----
names <-c("80ng", "320ng", "920ng")

# Create empty list to store dataframes
all_reads <- list()

# Loop through each condition size
for (condition_size in names) {
  file_path <- paste0("/data/share/htp/prime-seq_NextGen/data/FC2024_08_01_poolsize/03_zUMIs/", 
                     condition_size, "/zUMIs_output/stats/poolsize_", condition_size, ".readspercell.txt")
  
  # Read and process each file
  reads_data <- read.delim(file_path) %>%
    filter(type %in% c("Intergenic", "Exon", "Intron") & RG != "bad") %>%
    summarise(total_reads = sum(N), .by="RG") %>%
    mutate(condition = condition_size)  # Add condition size column
  
  all_reads[[condition_size]] <- reads_data
}

# Combine all dataframes
combined_reads <- bind_rows(all_reads) %>%
  mutate(condition = factor(condition, levels=c("80ng", "320ng", "920ng")))



combined_reads %>% 
  mutate(sample_comb = paste0(condition, RG)) %>% 
  ggplot(aes(y=total_reads, x=reorder(sample_comb, total_reads), fill=condition)) + 
  geom_col() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "Sample and Read Group", 
       y = "Total Reads", 
       title = "Total Reads per Sample and Read Group",
       subtitle = "Includes Intergenic, Exon, and Intron reads") +
  scale_fill_manual(values = pool_colors)

# Save the first plot
p0 <- combined_reads %>% 
  mutate(sample_comb = paste0(condition, RG)) %>% 
  ggplot(aes(y=total_reads, x=reorder(sample_comb, total_reads), fill=condition)) + 
  geom_col() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "Sample and Read Group", 
       y = "Total Reads", 
       title = "Total Reads per Sample and Read Group",
       subtitle = "Includes Intergenic, Exon, and Intron reads") +
  scale_fill_manual(values = pool_colors)
ggsave("/data/share/htp/prime-seq_NextGen/scripts/intergenics_analysis/plot/00_total_reads_per_condition_and_rg.pdf", 
       p0, width = 10, height = 6)


combined_reads <- combined_reads %>%
  summarise(total_reads = sum(total_reads), .by="condition")

# 2 plot normalized by intergenic reads ----
## 2.1 reads ----
# Sort data by count
sorted_data <- combined_data %>% 
  # First calculate total counts per condition for normalization
  group_by(condition) %>%
  mutate(
    normalized_reads = read_count / sum(read_count),
    normalized_umis = unique_umi_count / sum(unique_umi_count)
  ) %>%
  ungroup() %>%
  mutate(
    # Create numeric version of chr for finding max number
    chr_num = as.numeric(chr),
    # Get maximum numeric chromosome
    max_chr = max(chr_num[!is.na(chr_num)]),
    # Create properly ordered factor
    chr = factor(chr, 
                 levels = c(as.character(1:max_chr), "M", "X", "Y"),
                 ordered = TRUE)
  ) %>% 
  dplyr::select(-chr_num, -max_chr) %>%
  arrange(desc(normalized_reads)) %>%
  # Group by condition to assign tile_by_count within each group
  group_by(condition) %>%
  # Arrange by normalized_reads within each group (ascending)
  arrange(desc(normalized_reads), .by_group = TRUE) %>%
  # Add row numbers as tile_by_count (starting from 1 for lowest count)
  mutate(
    # Use row_number() to handle ties randomly
    tile_by_count = row_number(),
    # Calculate cumulative sum of normalized counts
    cumsum_reads = cumsum(normalized_reads),
    cumsum_umis = cumsum(normalized_umis)
  ) %>%
  # Restore original sorting
  ungroup() %>%
  arrange(desc(normalized_reads))

# Create line plot with facets by condition and chromosome
p1 <- ggplot(sorted_data %>% filter(tile_by_count <=50)
             , aes(y=normalized_reads, x=tile_by_count, color=condition)) +
  geom_line() +
  # facet_wrap(. ~ chr, scales = "free_x") +
  theme_bw() +
  labs(x = "Window rank", 
       y = "Normalized Reads", 
       title = "Window Counts (Normalized by Condition Intergenics)") +
  # scale_y_log10() +
  scale_color_manual(values = pool_colors)
p1
ggsave("/data/share/htp/prime-seq_NextGen/scripts/intergenics_analysis/plot/01_window_counts_normalized_by_intergenic.pdf", 
       p1, width = 8, height = 6)

# Create cumulative sum plot
p2 <- ggplot(sorted_data %>% filter(tile_by_count <=50)
             , aes(y=cumsum_reads, x=tile_by_count, color=condition)) +
  geom_line() +
  # facet_wrap(. ~ chr, scales = "free_x") +
  theme_bw() +
  labs(x = "Window rank", y = "Cumulative Normalized Count", 
       title = "Cumulative Window Counts (Normalized by Condition Total Intergenics)") +
  scale_color_manual(values = pool_colors)
p2
ggsave("/data/share/htp/prime-seq_NextGen/scripts/intergenics_analysis/plot/02_cumulative_window_counts_normalized_by_intergenic.pdf", 
       p2, width = 8, height = 6)


## 2.2 UMIs ----
# Sort data by count
sorted_data_umis <- combined_data %>% 
  # First calculate total counts per condition for normalization
  group_by(condition) %>%
  mutate(
    normalized_reads = read_count / sum(read_count),
    normalized_umis = unique_umi_count / sum(unique_umi_count)
  ) %>%
  ungroup() %>%
  mutate(
    # Create numeric version of chr for finding max number
    chr_num = as.numeric(chr),
    # Get maximum numeric chromosome
    max_chr = max(chr_num[!is.na(chr_num)]),
    # Create properly ordered factor
    chr = factor(chr, 
                 levels = c(as.character(1:max_chr), "M", "X", "Y"),
                 ordered = TRUE)
  ) %>% 
  dplyr::select(-chr_num, -max_chr) %>%
  arrange(desc(normalized_umis)) %>%
  # Group by condition to assign tile_by_count within each group
  group_by(condition) %>%
  # Arrange by normalized_umis within each group (ascending)
  arrange(desc(normalized_umis), .by_group = TRUE) %>%
  # Add row numbers as tile_by_count (starting from 1 for lowest count)
  mutate(
    # Use row_number() to handle ties randomly
    tile_by_count = row_number(),
    # Calculate cumulative sum of normalized counts
    cumsum_reads = cumsum(normalized_reads),
    cumsum_umis = cumsum(normalized_umis)
  ) %>%
  # Restore original sorting
  ungroup() %>%
  arrange(desc(normalized_umis))

# Create line plot with facets by condition and chromosome
p3 <- ggplot(sorted_data_umis %>% filter(tile_by_count <=50)
             , aes(y=normalized_umis, x=tile_by_count, color=condition)) +
  geom_line() +
  # facet_wrap(. ~ chr, scales = "free_x") +
  theme_bw() +
  labs(x = "Window rank", 
       y = "Normalized UMIs", 
       title = "Window Counts (Normalized by Condition Intergenics)") +
  # scale_y_log10() +
  scale_color_manual(values = pool_colors)
p3
ggsave("/data/share/htp/prime-seq_NextGen/scripts/intergenics_analysis/plot/03_window_counts_umis_normalized_by_intergenic.pdf", 
       p3, width = 8, height = 6)

# Create cumulative sum plot
p4 <- ggplot(sorted_data_umis %>% filter(tile_by_count <=50)
             , aes(y=cumsum_umis, x=tile_by_count, color=condition)) +
  geom_line() +
  # facet_wrap(. ~ chr, scales = "free_x") +
  theme_bw() +
  labs(x = "Window rank", y = "Cumulative Normalized UMIs", 
       title = "Cumulative Window Counts (Normalized by Condition Total Intergenics)") +
  scale_color_manual(values = pool_colors)
p4
ggsave("/data/share/htp/prime-seq_NextGen/scripts/intergenics_analysis/plot/04_cumulative_window_counts_umis_normalized_by_intergenic.pdf", 
       p4, width = 8, height = 6)


# 3 plot normalized by intergenic reads ----
## 3.1 reads ----
# Sort data by count
sorted_data2 <- combined_data %>% 
  left_join(combined_reads, by="condition") %>%
  mutate(
    normalized_reads = read_count / total_reads,
    normalized_umis = unique_umi_count / total_reads
  ) %>%
  ungroup() %>%
  mutate(
    # Create numeric version of chr for finding max number
    chr_num = as.numeric(chr),
    # Get maximum numeric chromosome
    max_chr = max(chr_num[!is.na(chr_num)]),
    # Create properly ordered factor
    chr = factor(chr, 
                 levels = c(as.character(1:max_chr), "M", "X", "Y"),
                 ordered = TRUE)
  ) %>% 
  dplyr::select(-chr_num, -max_chr) %>%
  arrange(desc(normalized_reads)) %>%
  # Group by condition to assign tile_by_count within each group
  group_by(condition) %>%
  # Arrange by normalized_reads within each group (ascending)
  arrange(desc(normalized_reads), .by_group = TRUE) %>%
  # Add row numbers as tile_by_count (starting from 1 for lowest count)
  mutate(
    # Use row_number() to handle ties randomly
    tile_by_count = row_number(),
    # Calculate cumulative sum of normalized counts
    cumsum_reads = cumsum(normalized_reads),
    cumsum_umis = cumsum(normalized_umis)
  ) %>%
  # Restore original sorting
  ungroup() %>%
  arrange(desc(normalized_reads))

# Create line plot with facets by condition and chromosome
p5 <- ggplot(sorted_data2 %>% filter(tile_by_count <=50)
             , aes(y=normalized_reads, x=tile_by_count, color=condition)) +
  geom_line() +
  # facet_wrap(. ~ chr, scales = "free_x") +
  theme_bw() +
  labs(x = "Window rank", 
       y = "Normalized Reads", 
       title = "Window Counts (Normalized by Condition In+Ex+Inter)") +
  # scale_y_log10() +
  scale_color_manual(values = pool_colors)
p5
ggsave("/data/share/htp/prime-seq_NextGen/scripts/intergenics_analysis/plot/05_window_counts_normalized_by_total.pdf", 
       p5, width = 8, height = 6)

# Create cumulative sum plot
p6 <- ggplot(sorted_data2 %>% filter(tile_by_count <=50)
             , aes(y=cumsum_reads, x=tile_by_count, color=condition)) +
  geom_line() +
  # facet_wrap(. ~ chr, scales = "free_x") +
  theme_bw() +
  labs(x = "Window rank", y = "Cumulative Normalized Count", 
       title = "Cumulative Window Counts (Normalized by Condition In+Ex+Inter)") +
  scale_color_manual(values = pool_colors)
p6
ggsave("/data/share/htp/prime-seq_NextGen/scripts/intergenics_analysis/plot/06_cumulative_window_counts_normalized_by_total.pdf", 
       p6, width = 8, height = 6)


## 3.2 UMIs ----
# Sort data by count
sorted_data2_umis <- combined_data %>% 
  left_join(combined_reads, by="condition") %>%
  mutate(
    normalized_reads = read_count / total_reads,
    normalized_umis = unique_umi_count / total_reads
  ) %>%
  ungroup() %>%
  mutate(
    # Create numeric version of chr for finding max number
    chr_num = as.numeric(chr),
    # Get maximum numeric chromosome
    max_chr = max(chr_num[!is.na(chr_num)]),
    # Create properly ordered factor
    chr = factor(chr, 
                 levels = c(as.character(1:max_chr), "M", "X", "Y"),
                 ordered = TRUE)
  ) %>% 
  dplyr::select(-chr_num, -max_chr) %>%
  arrange(desc(normalized_umis)) %>%
  # Group by condition to assign tile_by_count within each group
  group_by(condition) %>%
  # Arrange by normalized_umis within each group (ascending)
  arrange(desc(normalized_umis), .by_group = TRUE) %>%
  # Add row numbers as tile_by_count (starting from 1 for lowest count)
  mutate(
    # Use row_number() to handle ties randomly
    tile_by_count = row_number(),
    # Calculate cumulative sum of normalized counts
    cumsum_reads = cumsum(normalized_reads),
    cumsum_umis = cumsum(normalized_umis)
  ) %>%
  # Restore original sorting
  ungroup() %>%
  arrange(desc(normalized_umis))

# Create line plot with facets by condition and chromosome
p7 <- ggplot(sorted_data2_umis %>% filter(tile_by_count <=50)
             , aes(y=normalized_umis, x=tile_by_count, color=condition)) +
  geom_line() +
  # facet_wrap(. ~ chr, scales = "free_x") +
  theme_bw() +
  labs(x = "Window rank", 
       y = "Normalized UMIs", 
       title = "Window Counts (Normalized by Condition In+Ex+Inter)") +
  # scale_y_log10() +
  scale_color_manual(values = pool_colors)
p7
ggsave("/data/share/htp/prime-seq_NextGen/scripts/intergenics_analysis/plot/07_window_counts_umis_normalized_by_total.pdf", 
       p7, width = 8, height = 6)

# Create cumulative sum plot
p8 <- ggplot(sorted_data2_umis %>% filter(tile_by_count <=50)
             , aes(y=cumsum_umis, x=tile_by_count, color=condition)) +
  geom_line() +
  # facet_wrap(. ~ chr, scales = "free_x") +
  theme_bw() +
  labs(x = "Window rank", y = "Cumulative Normalized UMIs", 
       title = "Cumulative Window Counts (Normalized by Condition In+Ex+Inter)") +
  scale_color_manual(values = pool_colors)
p8
ggsave("/data/share/htp/prime-seq_NextGen/scripts/intergenics_analysis/plot/08_cumulative_window_counts_umis_normalized_by_total.pdf", 
       p8, width = 8, height = 6)


# 4 top 100 tiles per condition (by UMIs) ----
## 4.1 for each separate ----
top_each <- sorted_data2_umis %>%
  arrange(tile_by_count, condition)

# Save top 100 each table
write.table(top_each, 
            "/data/share/htp/prime-seq_NextGen/scripts/intergenics_analysis/plot/top_100_windows_per_condition.txt",
            sep="\t", quote=FALSE, row.names=FALSE)
   
## 4.2 for all ----
top_all <- sorted_data2_umis %>%
  select(c(chr, start, end, condition, normalized_umis)) %>%
  pivot_wider(names_from = condition, values_from = normalized_umis) %>%
  mutate(mean = rowMeans(across(c(`920ng`, `320ng`, `80ng`)))) %>%
  arrange(desc(mean)) %>%
  # Create a location field by merging chr, start, and end
  mutate(location = paste(chr, start, sep=":")) %>%
  mutate(location = paste(location, end, sep="-"))

# Save top all table
write.table(top_all %>%
              slice_head(n=100), 
            "/data/share/htp/prime-seq_NextGen/scripts/intergenics_analysis/plot/top_100_windows_mean_across_conditions.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

# Save the top 25 windows plot
p9 <- top_all %>%
  slice_head(n=25) %>%
  ggplot(aes(y=mean, x=reorder(location, mean), fill=chr))+
  geom_col() +
  coord_flip()+
  theme_bw() +
  labs(x = "Location", y = "Mean Normalized UMIs", 
       title = "Top Windows by Mean Normalized Count",
       subtitle = "Normalized by Condition In+Ex+Inter")
ggsave("/data/share/htp/prime-seq_NextGen/scripts/intergenics_analysis/plot/09_top_25_windows_mean_by_mean.pdf", 
       p9, width = 12, height = 8)

# Save the top 100 windows plot
p10 <- top_all %>%
  slice_head(n=50) %>%
  mutate(tile=c(1:50)) %>%
  pivot_longer(cols=4:6, names_to = "condition", values_to = "norm_count") %>%
  mutate(condition = factor(condition, levels=c("80ng", "320ng", "920ng")))%>%
  ggplot(aes(y=norm_count, x=tile, color=condition))+
  geom_line() +
  theme_bw() +
  labs(x = "Window rank (equal across conditions)", y = "Mean Normalized UMIs", 
       title = "Top Windows by Mean Normalized Count",
       subtitle = "Normalized by Condition In+Ex+Inter") +
  scale_color_manual(values = pool_colors)
ggsave("/data/share/htp/prime-seq_NextGen/scripts/intergenics_analysis/plot/10_top_50_windows_by_mean.pdf", 
       p10, width = 8, height = 6)
