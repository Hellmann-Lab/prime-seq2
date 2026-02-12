# respond to reviewers' comments on poolsize

# Load libraries
library(tidyverse)
library(here)
library(ggpubr)
library(ggbeeswarm)
library(rstatix)
library(patchwork)
library(readxl)
library(cowplot)

# Color scheme
color_vector <- c("80ng" = "#C4C082",
                  "320ng" = "#7D7A3B",
                  "920ng" = "#454321",
                  "80" = "#C4C082",
                  "320" = "#7D7A3B",
                  "920" = "#454321",
                  "80_1_1" = "#C4C082",
                  "920_1_1" = "#8EDEE1",
                  "920_1_4" = "#096581",
                  "920_4_4" = "#1C333F")
here::i_am("scripts/poolsize_figure.R")
# 1 funnel ----
## poolsize 1 (FC2024_08_01_poolsize) ----
# thats flowcell /data/sequencing/illumina/240819_VL00272_298_AACMFLNHV/
# the raw data is on the SFB server. unsure if the deML data is still somewhere ...
# read in data

### Mapping Fractions ----
file_paths <- list.files(path = here("data/FC2024_08_01_poolsize/03_zUMIs"),
                         pattern = "\\.readspercell\\.txt$", recursive=T)
file_paths

# Extract project name from directory structure (e.g., "80ng/zUMIs_output/stats/poolsize_80ng.readspercell.txt" -> "80ng")
projects <- str_extract(file_paths, "(?<=^|/)(\\d+)")  # keep only numeric part
names <- projects

rpc_all_smpl_ps1 <- map_df(seq_along(file_paths), function(i) {
  read.csv(here("data/FC2024_08_01_poolsize/03_zUMIs", file_paths[i]), sep= "\t") %>%
    filter(RG != "bad") %>%
    mutate(project = projects[i]) %>%
    group_by(RG) %>%
    mutate(sum = sum(N)) %>%
    ungroup() %>%
    mutate(fraction = N/sum) %>%
    group_by(project) %>%
    mutate(sum_project = sum(N)) %>%
    ungroup() %>%
    group_by(type) %>%
    mutate(sum_type = sum(N)) %>%
    ungroup() %>%
    mutate(fraction_type_project = sum_type/sum_project) %>%
    dplyr::select(c(RG, N, type, project, fraction, fraction_type_project))
}) %>%
  #filter(type != "Unmapped") %>%
  mutate(type = factor(type, levels=c("Exon", "Intron", "Intergenic", "Ambiguity", "Unmapped")))

### deML ----
ass_reads_ps1 <- read.delim(file=here("data/FC2024_08_01_poolsize/deML_summary_to_read_in.txt")) %>%
  dplyr::select(RG, assigned, total) %>%
  dplyr::rename("project"="RG") %>%
  filter(str_detect(project, "primeseq_poolsize")) %>%
  mutate(project = str_extract(project, "(?<=_poolsize_).*")) %>% 
  mutate(non_assigned = total-assigned, 
         fract_non_assigned = paste(round(non_assigned/total, digits=3)*100, "%")) %>%
  pivot_longer(cols=c(2,4), names_to = "category", values_to = "reads") %>%
  mutate(fract_non_assigned = ifelse(category == "non_assigned", fract_non_assigned, NA))

### Trimming data ----
trim_df_ps1 <- map_df(projects, function(n) {
  read_delim(here(paste0("data/FC2024_08_01_poolsize/02_trimming/", n, "ng.txt")),
             col_names = c("category", "reads")) %>%
    mutate(reads = str_remove_all(reads, "[ ,%\\)]")) %>%
    separate(reads, into = c("reads", "percentage"), sep = "\\(", convert = TRUE) %>%
    mutate(condition = n)
}) %>%
  mutate(reads = as.numeric(reads))

### Barcodes assignment ----
file_paths1 <- list.files(path = here("data/FC2024_08_01_poolsize/03_zUMIs"),
                          pattern = "kept_barcodes_binned\\.txt$", recursive = TRUE)

nreads_barcodes_ps1 <- map_df(seq_along(file_paths1), function(i) {
  barcodes <- read.csv(here("data/FC2024_08_01_poolsize/03_zUMIs", file_paths1[i]))
  barcodes %>%
    mutate(project = projects[i])
}) %>%
  bind_rows()

nreads_ps1 <- nreads_barcodes_ps1 %>%
  group_by(project) %>%
  summarise('Barcode assigned' = sum(n))

### Read UMI data ----
read_umi_raw_ps1 <- map_df(projects, function(n) {
  gene_counts <- read_rds(here(paste0("data/FC2024_08_01_poolsize/03_zUMIs/", n, "ng/zUMIs_output/stats/poolsize_", n, "ng.bc.READcounts.rds"))) %>%
    dplyr::rename(read_count = N, SampleID = RG)
  umi_counts <- read_delim(here(paste0("data/FC2024_08_01_poolsize/03_zUMIs/", n, "ng/zUMIs_output/stats/poolsize_", n, "ng.UMIcounts.txt"))) %>%
    dplyr::rename(umi_count = Count)
  
  full_join(gene_counts, umi_counts, by = c("SampleID", "type")) %>%
    mutate(project = n)
}) 
read_umi_ps1 <- read_umi_raw_ps1 %>%
  filter(type %in% c("Exon", "Intron") & SampleID != "bad") %>%
  mutate(umi_fraction = umi_count / read_count)

read_umi_summary_ps1 <- read_umi_ps1 %>%
  group_by(project) %>%
  summarize(UMI = sum(umi_count), .groups = 'drop')


##  poolsize 2 (FC2025_06_01_poolsize2) ----
# read in data

# only demultiplexed data

### Mapping Fractions ----
file_paths <- list.files(path = here("data/FC2025_06_01_poolsize2/03_zUMIs"),
                         pattern = "\\.readspercell\\.txt$", recursive=T)
file_paths

# Extract project name from directory structure (e.g., "80_1_1/zUMIs_output/stats/poolsize2_80_1_1.readspercell.txt" -> "80_1_1")
projects <- str_extract(file_paths, "^[^/]+")
names <- projects

rpc_all_smpl_ps2 <- map_df(seq_along(file_paths), function(i) {
  read.csv(here("data/FC2025_06_01_poolsize2/03_zUMIs", file_paths[i]), sep= "\t") %>%
    filter(RG != "bad") %>%
    mutate(project = projects[i]) %>%
    group_by(RG) %>%
    mutate(sum = sum(N)) %>%
    ungroup() %>%
    mutate(fraction = N/sum) %>%
    group_by(project) %>%
    mutate(sum_project = sum(N)) %>%
    ungroup() %>%
    group_by(type) %>%
    mutate(sum_type = sum(N)) %>%
    ungroup() %>%
    mutate(fraction_type_project = sum_type/sum_project) %>%
    dplyr::select(c(RG, N, type, project, fraction, fraction_type_project))
}) %>%
  #filter(type != "Unmapped") %>%
  mutate(type = factor(type, levels=c("Exon", "Intron", "Intergenic", "Ambiguity", "Unmapped")))

### Trimming data ----
trim_df_ps2 <- map_df(projects, function(n) {
  read_delim(here(paste0("data/FC2025_06_01_poolsize2/02_trimming/", n, ".txt")),
             delim = ":", 
             col_names = c("category", "reads_raw"),
             trim_ws = TRUE) %>%
    mutate(reads = str_extract(reads_raw, "[0-9,]+")) %>%
    mutate(reads = str_remove_all(reads, ",")) %>%
    mutate(reads = as.numeric(reads)) %>%
    mutate(percentage = str_extract(reads_raw, "\\([0-9.]+%\\)")) %>%
    dplyr::select(category, reads, percentage) %>%
    mutate(condition = n)
}) %>%
  mutate(reads = as.numeric(reads))

### Barcodes assignment ----
file_paths2 <- list.files(path = here("data/FC2025_06_01_poolsize2/03_zUMIs"),
                          pattern = "kept_barcodes_binned\\.txt$", recursive = TRUE)

nreads_barcodes_ps2 <- map_df(seq_along(file_paths2), function(i) {
  barcodes <- read.csv(here("data/FC2025_06_01_poolsize2/03_zUMIs", file_paths2[i]))
  barcodes %>%
    mutate(project = projects[i])
}) %>%
  bind_rows()

nreads_ps2 <- nreads_barcodes_ps2 %>%
  group_by(project) %>%
  summarise('Barcode assigned' = sum(n))

### Read UMI data ---- 
read_umi_raw_ps2 <- map_df(projects, function(n) {
  gene_counts <- read_rds(here(paste0("data/FC2025_06_01_poolsize2/03_zUMIs/", n, "/zUMIs_output/stats/poolsize2_", n, ".bc.READcounts.rds"))) %>%
    dplyr::rename(read_count = N, SampleID = RG)
  umi_counts <- read_delim(here(paste0("data/FC2025_06_01_poolsize2/03_zUMIs/", n, "/zUMIs_output/stats/poolsize2_", n, ".UMIcounts.txt"))) %>%
    dplyr::rename(umi_count = Count)
  
  full_join(gene_counts, umi_counts, by = c("SampleID", "type")) %>%
    mutate(project = n)
}) 
read_umi_ps2 <- read_umi_raw_ps2 %>%
  filter(type %in% c("Exon", "Intron") & SampleID != "bad") %>%
  mutate(umi_fraction = umi_count / read_count)

read_umi_summary_ps2 <- read_umi_ps2 %>%
  group_by(project) %>%
  summarize(UMI = sum(umi_count), .groups = 'drop')

## combine funnel data ----
### poolsize 1 ----
read_funnel_ps1  <-
  ass_reads_ps1 %>%
  filter(category == "assigned") %>%
  dplyr::select(project,  'Total Reads' = total, 'Index assigned' = reads) %>%
  full_join(
    trim_df_ps1 %>%
      dplyr::select(-percentage) %>%
      tidyr::pivot_wider(names_from = category, values_from = reads) %>%
      group_by(condition) %>%
      mutate('Trimmed & Filtered' = `Total read pairs processed` - `Pairs that were too short`) %>%
      dplyr::select(project = condition, 'Trimmed & Filtered'), 
    by="project"
  ) %>%
  full_join(
    nreads_ps1,
    by="project"
  ) %>%
  pivot_longer(cols=-1, names_to = "step", values_to = "reads") %>%
  bind_rows(rpc_all_smpl_ps1 %>% 
              dplyr::select(BC = RG, N, type, project) %>%
              ungroup() %>% 
              group_by(project, BC) %>% 
              filter(type %in% c("Intron", "Exon")) %>% 
              summarise('Exon & Intron selected'=sum(N), .groups = 'drop') %>%
              full_join(read_umi_ps1 %>% 
                          dplyr::select(BC = SampleID, UMI= umi_count, type, project) %>%
                          group_by(project, BC) %>%
                          filter(type %in% c("Intron", "Exon")) %>% 
                          summarise('UMI collapsed'=sum(UMI), .groups = 'drop'),
                        by=c("project", "BC")) %>%
              pivot_longer(cols=-c(BC, project), names_to = "step", values_to = "reads")
  )  %>%
  ungroup()%>%
  mutate(dataset = "Experiment 1",
         poolsize = project) # keep original project name as poolsize

### poolsize 2 (no deML data available) ----
# Extract pool size from project names (e.g., "80_1_1" -> "80ng", "920_1_1" -> "920ng")
read_funnel_ps2  <-
  trim_df_ps2 %>%
  dplyr::select(-percentage) %>%
  tidyr::pivot_wider(names_from = category, values_from = reads) %>%
  group_by(condition) %>%
  mutate('Trimmed & Filtered' = `Total read pairs processed` - `Pairs that were too short`) %>%
  mutate('Index assigned' = `Total read pairs processed`) %>%
  dplyr::select(project = condition, 'Trimmed & Filtered', 'Index assigned') %>%
  full_join(
    nreads_ps2,
    by="project"
  ) %>%
  pivot_longer(cols=-1, names_to = "step", values_to = "reads") %>%
  bind_rows(rpc_all_smpl_ps2 %>% 
              dplyr::select(BC = RG, N, type, project) %>%
              ungroup() %>% 
              group_by(project, BC) %>% 
              filter(type %in% c("Intron", "Exon")) %>% 
              summarise('Exon & Intron selected'=sum(N), .groups = 'drop') %>%
              full_join(read_umi_ps2 %>% 
                          dplyr::select(BC = SampleID, UMI= umi_count, type, project) %>%
                          group_by(project, BC) %>%
                          filter(type %in% c("Intron", "Exon")) %>% 
                          summarise('UMI collapsed'=sum(UMI), .groups = 'drop'),
                        by=c("project", "BC")) %>%
              pivot_longer(cols=-c(BC, project), names_to = "step", values_to = "reads")
  )  %>%
  ungroup()%>%
  mutate(dataset = "Experiment 2",
         poolsize = str_extract(project, "^\\d+")) # keep original project name as poolsize

### combine both datasets ----
read_funnel_all_ps <- bind_rows(read_funnel_ps1, read_funnel_ps2) %>%
  mutate(step = factor(step, levels = c("Total Reads", "Index assigned", "Trimmed & Filtered", "Barcode assigned", "Exon & Intron selected", "UMI collapsed"))) %>%
  filter(step != "Total Reads") %>%
  group_by(project, dataset, poolsize) %>%
  mutate(max_dataset = max(reads, na.rm = TRUE),
         reads_norm = reads / max_dataset) %>% 
  mutate(reads_norm = reads_norm * ifelse(!is.na(BC), case_when(project %in% c("80", "80_1_1") ~ 8, 
                                                                project %in% c("320", "920", "920_1_1", "920_1_4", "920_4_4") ~ 16)
                                          , 1)) %>% # adjust if needed
  ungroup() %>%
  mutate(poolsize = factor(poolsize, levels=c("80", "320", "920"))) 

### calculate average (using poolsize for grouping instead of project) ----
# read_funnel_avg_all_ps <- read_funnel_all_ps %>%
# ungroup() %>%
# group_by(poolsize, step, dataset) %>%
# summarise(reads = sum(reads, na.rm = TRUE), .groups = 'drop') %>% 
# group_by(dataset, poolsize) %>%
# mutate(reads_rel_dataset = reads / max(reads, na.rm = TRUE)) %>%
# group_by(poolsize, step) %>%
# summarise(reads_rel = mean(reads_rel_dataset), .groups = 'drop')

read_funnel_avg_all_ps <- read_funnel_all_ps %>%
  ungroup() %>%
  group_by(project, step, dataset) %>%
  summarise(reads = sum(reads, na.rm = TRUE), .groups = 'drop') %>% 
  group_by(dataset, project) %>%
  mutate(reads_rel_dataset = reads / max(reads, na.rm = TRUE)) %>%
  group_by(dataset, project, step) %>%
  summarise(reads_rel = mean(reads_rel_dataset), .groups = 'drop')


## Funnel change ----
read_funnel_dif1 <- bind_rows(read_funnel_all_ps, 
                              bind_rows(nreads_barcodes_ps1 %>% select(BC = XC, reads = n, project) %>%
                                          mutate(dataset = "Experiment 1",
                                                 step = "Barcode assigned"),
                                        nreads_barcodes_ps2 %>% select(BC = XC, reads = n, project) %>%
                                          mutate(dataset = "Experiment 2",
                                                 step = "Barcode assigned")
                              ) %>%
                                mutate(step = factor(step, levels = c("Total Reads", "Index assigned", "Trimmed & Filtered", "Barcode assigned", "Exon & Intron selected", "UMI collapsed")))
) %>%
  arrange(dataset, step, BC)

read_funnel_dif <- read_funnel_dif1 %>%
  group_by(project, dataset, BC) %>%
  mutate(reads_before = lag(reads)) %>%
  left_join(read_funnel_all_ps %>% filter(step == "Trimmed & Filtered") %>% select(project, dataset, trimmed_reads = reads),
            by = c("project", "dataset"),
            relationship = "many-to-one") %>%
  filter(!is.na(max_dataset)) %>%
  mutate(step_change = ifelse(!is.na(reads_before), reads/reads_before, reads_norm))


read_funnel_dif_avg <- read_funnel_dif %>%
  group_by(project, dataset,step) %>%
  summarise(mean = mean(step_change-1))

funnel_change_ps1 <- ggplot()+
  geom_col(read_funnel_dif_avg %>% filter(dataset == "Experiment 1") %>%
             mutate(project = factor(project, levels=c("80", "320", "920"))),
           mapping = aes(y=mean, x=project, fill=project))+
  stat_summary(data=read_funnel_dif %>% filter(dataset == "Experiment 1")%>%
                 mutate(project = factor(project, levels=c("80", "320", "920"))), 
               aes(y=step_change-1, x=project),
               color = "black",
               fun.data = "mean_se",
               fun.args = list(mult = 1),
               geom="errorbar", 
               #position=position_dodge(0.2), 
               width=.2,
               size=.2) + 
  coord_flip()+
  facet_grid(step~.)+
  theme_pubr(legend = "none")+
  scale_fill_manual(values = color_vector)+
  scale_x_discrete(limits=rev)+
  labs(x = "", 
       y = "Fraction of previous step")+  
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y =element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank())+
  scale_y_continuous(position = "right", breaks = c(0.0, -0.2, -0.4), limits = c(-0.35, 0))+
  geom_hline(yintercept = 0, linetype="dashed")

funnel_change_ps2 <- ggplot()+
  geom_col(read_funnel_dif_avg %>% filter(dataset == "Experiment 2"),
           mapping = aes(y=mean, x=project, fill=project))+
  stat_summary(data=read_funnel_dif %>% filter(dataset == "Experiment 2"), 
               aes(y=step_change-1, x=project),
               color = "black",
               fun.data = "mean_se",
               fun.args = list(mult = 1),
               geom="errorbar", 
               #position=position_dodge(0.2), 
               width=.2,
               size=.2) + 
  coord_flip()+
  facet_grid(step~.)+
  theme_pubr(legend = "none")+
  scale_fill_manual(values = color_vector)+
  scale_x_discrete(limits=rev)+
  labs(x = "", 
       y = "Fraction of previous step")+  
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y =element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank())+
  scale_y_continuous(position = "right", breaks = c(0.0, -0.2, -0.4), limits = c(-0.35, 0))+
  geom_hline(yintercept = 0, linetype="dashed")

funnel_change_ps1+funnel_change_ps2
## create funnel plot ----
funnel_plot_ps1 <- ggplot() +
  stat_summary(data = read_funnel_all_ps %>% filter(dataset == "Experiment 1")%>%
                 mutate(project = factor(project, levels=c("80", "320", "920"))), 
               aes(y = step, x = reads_norm, color = project), 
               fun.data = "mean_se",
               fun.args = list(mult = 1),
               geom = "errorbar", 
               position = position_dodge(0.2), 
               width = .2,
               size = .5) + 
  geom_line(data = read_funnel_avg_all_ps %>% filter(dataset == "Experiment 1")%>%
              mutate(project = factor(project, levels=c("80", "320", "920"))), 
            aes(y = step, x = reads_rel, color = project, group = project), 
            size = 0.9) + 
  scale_x_continuous(position = "top")+ 
  coord_cartesian(xlim = c(0.35, 1))+
  ylab("Processing Steps \n \u27F5") +
  xlab("Relative Reads") +
  scale_y_discrete(limits = rev(levels(read_funnel_all_ps$step))) +
  labs(color = "Pool size") +
  theme_pubr(legend = "bottom") +
  scale_color_manual(values = color_vector) #+
# theme(axis.text.y = element_text(size = 13),
#       legend.text = element_text(size = 13),
#       legend.title = element_text(size = 13))
funnel_plot_ps2 <- ggplot() +
  stat_summary(data = read_funnel_all_ps %>% filter(dataset == "Experiment 2"), 
               aes(y = step, x = reads_norm, color = project), 
               fun.data = "mean_se",
               fun.args = list(mult = 1),
               geom = "errorbar", 
               position = position_dodge(0.2), 
               width = .2,
               size = .5) + 
  geom_line(data = read_funnel_avg_all_ps %>% filter(dataset == "Experiment 2"), 
            aes(y = step, x = reads_rel, color = project, group = project), 
            size = 0.9) + 
  scale_x_continuous(position = "top") +
  coord_cartesian(xlim = c(0.35, 1))+
  ylab("Processing Steps \n \u27F5") +
  xlab("Relative Reads") +
  scale_y_discrete(limits = rev(levels(read_funnel_all_ps$step))) +
  labs(color = "Pool size") +
  theme_pubr(legend = "bottom") +
  scale_color_manual(values = color_vector) #+
# theme(axis.text.y = element_text(size = 13),
#       legend.text = element_text(size = 13),
#       legend.title = element_text(size = 13))

funnel_plot_ps1+funnel_plot_ps2

# # 2 UMI chimera ----
# # Load chimera data for Experiment 1
# pool_sizes <- c("80", "320", "920")
# ds_status <- c("no_ds", "ds")
# 
# all_data_chimera <- data.frame()
# 
# for (ps in pool_sizes) {
#   for (ds in ds_status) {
#     if (ds == "ds") {
#       filename <- paste0("poolsize_", ps, "ng_ds.UMI_chimera_per_bc.rds")
#     } else {
#       filename <- paste0("poolsize_", ps, "ng.UMI_chimera_per_bc.rds")
#     }
#     
#     filepath <- here("scripts/poolsize_chimera", filename)
#     
#     if (file.exists(filepath)) {
#       data <- readRDS(filepath)
#       data$poolsize <- ps
#       data$downsampling <- ds
#       all_data_chimera <- rbind(all_data_chimera, data)
#     }
#   }
# }
# 
# # Convert to factors for proper ordering
# all_data_chimera$poolsize <- factor(as.character(all_data_chimera$poolsize), levels = c("80", "320", "920"))
# all_data_chimera$downsampling <- factor(all_data_chimera$downsampling, levels = c("no_ds", "ds"), 
#                                         labels = c("No downsampling", "With downsampling"))
# 
# # Create beeswarm plot for Experiment 1
# plot_beeswarm_chimera1 <- ggplot(all_data_chimera, aes(x = poolsize, y = fract, color = poolsize)) +
#   geom_quasirandom(alpha = 0.5) +
#   stat_summary(fun.data = mean_se, geom = "errorbar",
#                width = 0) +
#   stat_summary(fun = mean, geom = "point",
#                shape = 23, size = 2, fill = "white") +
#   stat_compare_means(method = "wilcox.test", 
#                      label = "p.signif",
#                      comparisons = list(c("80", "320"), c("80", "920"), c("320", "920")),
#                      tip.length = 0) +
#   facet_wrap(~ downsampling, ncol = 2) +
#   xlab("Pool size (ng)") +
#   ylab("Chimera fraction per barcode") +
#   theme_pubr() +
#   theme(#axis.text.x = element_text(size=13),
#     #       strip.text.x = element_text(size=13),
#     #       plot.title = element_text(size=14, face="bold"),
#     legend.position = "none") +
#   scale_color_manual(values = color_vector)
# plot_beeswarm_chimera1
# 
# # Load chimera data for Experiment 2
# poolsize2_files <- c("poolsize2_80_1_1", "poolsize2_920_1_1", "poolsize2_920_1_4", "poolsize2_920_4_4")
# 
# all_data_chimera2 <- data.frame()
# 
# for (filename_base in poolsize2_files) {
#   filename <- paste0(filename_base, ".UMI_chimera_per_bc.rds")
#   filepath <- here("scripts/poolsize_chimera", filename)
#   
#   if (file.exists(filepath)) {
#     data <- readRDS(filepath)
#     # Extract pool size from filename
#     poolsize_match <- regmatches(filename_base, regexpr("_(80|920)_", filename_base))
#     poolsize_val <- as.numeric(gsub("_", "", poolsize_match))
#     data$poolsize <- poolsize_val
#     data$sample <- filename_base
#     all_data_chimera2 <- rbind(all_data_chimera2, data)
#   }
# }
# 
# # Convert to factors for proper ordering
# all_data_chimera2$poolsize <- factor(all_data_chimera2$poolsize, levels = c(80, 920))
# all_data_chimera2$sample <- factor(all_data_chimera2$sample)
# 
# # Create beeswarm plot for Experiment 2
# plot_beeswarm_chimera2 <- ggplot(all_data_chimera2, aes(x = poolsize, y = fract, color = poolsize)) +
#   geom_quasirandom(alpha = 0.5) +
#   stat_summary(fun.data = mean_se, geom = "errorbar",
#                width = 0) +
#   stat_summary(fun = mean, geom = "point",
#                shape = 23, size = 2, fill = "white") +
#   stat_compare_means(method = "wilcox.test", 
#                      label = "p.signif",
#                      comparisons = list(c("80", "920")),
#                      tip.length = 0) +
#   xlab("Pool size (ng)") +
#   ylab("Chimera fraction per barcode") +
#   theme_pubr() +
#   theme(#axis.text.x = element_text(size=13),
#     #       strip.text.x = element_text(size=13),
#     #       plot.title = element_text(size=14, face="bold"),
#     legend.position = "none") +
#   scale_color_manual(values = color_vector)
# plot_beeswarm_chimera2

# 3 Mapping stats ----
# Read STAR logs
star_logs <- read_excel(here("scripts/poolsize_intergenic_quality/STARlogs.xlsx")) %>%
  mutate(project = factor(project, levels = c("80", "320", "920", "80_1_1", "920_1_1", "920_1_4", "920_4_4"))) %>%
  mutate(dataset = ifelse(dataset == "poolsize1", "Experiment 1", "Experiment 2"))

# Read mismatch/quality/softclip data
# Experiment 1 (Experiment 1)
quality_exp1 <- readRDS(here("scripts/poolsize_intergenic_quality/intergenic_quality_exp1.rds")) %>%
  mutate(dataset = "Experiment 1") %>%
  mutate(poolsize = str_remove(as.character(poolsize), "ng")) %>%
  mutate(sample = poolsize)

# Experiment 2 (Experiment 2)
quality_exp2 <- readRDS(here("scripts/poolsize_intergenic_quality/intergenic_quality_exp2.rds")) %>%
  mutate(dataset = "Experiment 2")

# Prepare data for mapping stats plots
star_short_plot <- ggplot(star_logs %>% filter(type=="short"), aes(x=project, y=fraction, fill=project)) +
  geom_col()+
  # stat_summary(fun = mean, geom = "crossbar")+
  theme_pubr()+
  labs(x = "", color = "", y="Fraction of reads too short\n(STAR)")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #strip.text.x = element_text(size=13),
        #plot.title = element_text(size=14, face="bold"),
        legend.position = "none") +
  # scale_x_continuous(
  #   breaks = c(80, 320, 920),
  #   labels = c("80", "320", "920")
  # )+
  facet_grid(.~dataset, scales = "free", space = "free")+
  scale_fill_manual(values = color_vector)
star_short_plot

softclip <- bind_rows(quality_exp1, quality_exp2) %>%
  group_by(dataset, sample, poolsize) %>%
  mutate(softclipyesno = ifelse(soft_clip_fraction > 0, "yes", "no")) %>%
  dplyr::count(softclipyesno) %>%
  mutate(softclip_fraction = n / sum(n)) %>%
  filter(softclipyesno == "yes") %>%
  mutate(poolsize = as.numeric(poolsize))%>%
  mutate(sample = factor(sample, levels = c("80", "320", "920", "80_1_1", "920_1_1", "920_1_4", "920_4_4")))

softclip_plot <- ggplot(softclip, aes(x=sample, y=softclip_fraction, fill=sample)) +
  geom_col()+
  # stat_summary(fun = mean, geom = "crossbar")+
  theme_pubr()+
  labs(x = "", color = "", y="Fraction of intergenic reads\nwith soft-clipping")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #strip.text.x = element_text(size=13),
        #plot.title = element_text(size=14, face="bold"),
        legend.position = "none") +
  # scale_x_continuous(
  #   breaks = c(80, 320, 920),
  #   labels = c("80", "320", "920")
  # )+
  facet_grid(.~dataset, scales = "free", space = "free")+
  scale_fill_manual(values = color_vector)
softclip_plot

mismatch <- bind_rows(quality_exp1, quality_exp2) %>%
  group_by(dataset, sample, poolsize) %>%
  mutate(mismatchyesno = ifelse(mismatch_rate > 0, "yes", "no")) %>%
  dplyr::count(mismatchyesno) %>%
  mutate(mismatch_fraction = n / sum(n))%>%
  filter(mismatchyesno == "yes") %>%
  mutate(poolsize = as.numeric(poolsize))%>%
  mutate(sample = factor(sample, levels = c("80", "320", "920", "80_1_1", "920_1_1", "920_1_4", "920_4_4")))

mismatch_plot <- ggplot(mismatch, aes(x=sample, y=mismatch_fraction, fill=sample)) +
  geom_col()+
  theme_pubr()+
  labs(x = "", color = "", y="Fraction of intergenic reads\nwith mismatches")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #strip.text.x = element_text(size=13),
        #plot.title = element_text(size=14, face="bold"),
        legend.position = "none") +
  # scale_x_continuous(
  #   breaks = c(80, 320, 920),
  #   labels = c("80", "320", "920")
  # )+
  facet_grid(.~dataset, scales = "free", space = "free")+
  scale_fill_manual(values = color_vector)
mismatch_plot



# 4 Example ----
# todo

# Combine all plots into multipanel figure ----
# Extract legends before removing them from plots
# funnel_legend <- as_ggplot(ggpubr::get_legend(funnel_plot_ps))
# mapping_legend <- as_ggplot(ggpubr::get_legend(star_short_plot + theme(legend.position = "bottom")))
# 
# # Remove legends from individual plots for cleaner multipanel
# funnel_plot_ps <- funnel_plot_ps + theme(legend.position = "none")

# Combine all plots with labels
# Layout:
# Row 1: Funnel plot (full width)
# Row 2: Chimera plots (Experiment 1 and Experiment 2 side by side)
# Row 3: Mapping stats (short reads, softclip, mismatch)
layout <- c("ABC\nADE\nFGH")
spacer_labeled <- ggplot() +
  annotate("text", x = 0, y = 0, label = "", size = 5) +
  theme_void()
combined_multipanel <- spacer_labeled +
  funnel_plot_ps1 + 
  funnel_change_ps1 +
  funnel_plot_ps2 + 
  funnel_change_ps2 +
  star_short_plot + 
  softclip_plot + 
  mismatch_plot +
  plot_layout(heights = c(1, 1, 1),
              design = layout) +
  plot_annotation(tag_levels = list(c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')))
combined_multipanel

ggsave(combined_multipanel,
       filename=here("figures/fig_S13_poolsize.pdf"),
       height=10,
       width=13,
       device = cairo_pdf)
