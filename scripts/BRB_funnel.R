library(tidyverse)
library(data.table)
library(ggbeeswarm)
library(purrr)
library(ggpubr)
library(readxl)
library(patchwork)
library(cowplot)
library(ShortRead)


# read in data
#------------------------------------------------------------------------------------------------------------------------------------------------------------
# data1 ----
## Mapping Fractions
setwd("/data/share/htp/prime-seq_NextGen/other_methods/BRBseq/03_zUMIs")
file_paths <- list.files(pattern = "\\.readspercell\\.txt$", recursive=T)
file_paths

projects <- str_extract(file_paths, "(......)(?=\\.readspercell\\.txt)")
names <- projects

rpc_all_smpl <- map_df(seq_along(file_paths), function(i) {
  read.csv(file_paths[i], sep= "\t") %>%
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
    select(c(RG, N, type, project, fraction, fraction_type_project))
}) %>%
  #filter(type != "Unmapped") %>%
  mutate(type = factor(type, levels=c("Exon", "Intron", "Intergenic", "Ambiguity", "Unmapped")))

## deML (no deML data)
# ass_reads <- read.delim(file="/data/share/htp/prime-seq_NextGen/data/FC2024_05_02_PoP96_BA/deML_summary_to_read_in.txt") %>%
#   select(RG, assigned, total) %>%
#   dplyr::rename("project"="RG") %>%
#   filter(str_detect(project, "primeseq")) %>%
#   mutate(project = str_extract(project, "...$")) %>%
#   mutate(across(project, factor, levels=c("old", "new", "PTO"))) %>%
#   mutate(non_assigned = total-assigned, 
#          fract_non_assigned = paste(round(non_assigned/total, digits=3)*100, "%")) %>%
#   pivot_longer(cols=c(2,4), names_to = "category", values_to = "reads") %>%
#   mutate(fract_non_assigned = ifelse(category == "non_assigned", fract_non_assigned, NA))

## Trimming data
trim_df <- map_df(projects, function(n) {
  read_delim(paste0("/data/share/htp/prime-seq_NextGen/other_methods/BRBseq/02_trimming/", n, ".txt"),
             col_names = c("category", "reads")) %>%
    mutate(reads = str_remove_all(reads, "[ ,%\\)]")) %>%
    separate(reads, into = c("reads", "percentage"), sep = "\\(", convert = TRUE) %>%
    mutate(condition = n)
}) %>%
  mutate(reads = as.numeric(reads))

## Barcodes assignment
file_paths1 <- list.files(pattern = "kept_barcodes_binned\\.txt$", recursive = TRUE)

nreads_barcodes <- map_df(seq_along(file_paths1), function(i) {
  barcodes <- read.csv(file_paths1[i])
  barcodes %>%
    mutate(project = projects[i])
}) %>%
  bind_rows()

nreads <- nreads_barcodes %>%
  group_by(project) %>%
  summarise(barcode_assigned = sum(n))


## Read UMI data
read_umi_raw <- map_df(projects, function(n) {
  gene_counts <- read_rds(paste0("/data/share/htp/prime-seq_NextGen/other_methods/BRBseq/03_zUMIs/", n, "/zUMIs_output/stats/", n, ".bc.READcounts.rds")) %>%
    dplyr::rename(read_count = N, SampleID = RG)
  umi_counts <- read_delim(paste0("/data/share/htp/prime-seq_NextGen/other_methods/BRBseq/03_zUMIs/", n, "/zUMIs_output/stats/", n, ".UMIcounts.txt")) %>%
    dplyr::rename(umi_count = Count)
  
  full_join(gene_counts, umi_counts, by = c("SampleID", "type")) %>%
    mutate(project = n)
}) 
read_umi <- read_umi_raw %>%
  filter(type %in% c("Exon", "Intron") & SampleID != "bad") %>%
  mutate(umi_fraction = umi_count / read_count)

read_umi_summary <- read_umi %>%
  group_by(project) %>%
  summarize(UMI = sum(umi_count), .groups = 'drop')

# combine data
read_funnel_1  <-
  trim_df %>%
  select(-percentage) %>%
  pivot_wider(names_from = category, values_from = reads) %>%
  group_by(condition) %>%
  mutate(trimmed = `Total read pairs processed` - `Pairs that were too short`) %>%
  mutate(index_assigned = `Total read pairs processed`) %>%
  select(project = condition, index_assigned, trimmed) %>%
  full_join(
    nreads,
    by="project"
  ) %>%
  pivot_longer(cols=-1, names_to = "step", values_to = "reads") %>%
  bind_rows(rpc_all_smpl %>% 
              select(BC = RG, N, type, project) %>%
              ungroup() %>% 
              group_by(project, BC) %>% 
              filter(type %in% c("Intron", "Exon")) %>% 
              summarise(inex=sum(N), .groups = 'drop') %>%
              full_join(read_umi %>% 
                          select(BC = SampleID, UMI= umi_count, type, project) %>%
                          group_by(project, BC) %>%
                          filter(type %in% c("Intron", "Exon")) %>% 
                          summarise(UMI=sum(UMI), .groups = 'drop'),
                        by=c("project", "BC")) %>%
              pivot_longer(cols=-c(BC, project), names_to = "step", values_to = "reads")
  )  %>%
  ungroup() %>%
  mutate(step = factor(step, levels = c("total", "index_assigned", "trimmed", "barcode_assigned", "inex", "UMI"))) %>%
  mutate(max_rep=max(reads), 
         reads_norm = reads/max_rep) %>% # norm to total reads per  rep
  mutate(reads_norm = reads_norm*ifelse(!is.na(BC), 4, 1))

read_funnel_avg_1 <- read_funnel_1 %>%
  ungroup() %>%
  group_by(project, step) %>%
  summarise(reads=sum(reads)) %>% # summarise per step
  ungroup() %>%
  mutate(reads_rel = reads/max(reads))

#------------------------------------------------------------------------------------------------------------------------------------------------------------
# data2 ----
## Mapping Fractions
setwd("/data/share/htp/prime-seq_NextGen/other_methods/BRBseq_degraded_2mn/03_zUMIs")
file_paths <- list.files(pattern = "\\.readspercell\\.txt$", recursive=T)
file_paths

projects <- str_extract(file_paths, "([^/]*)(?=\\.readspercell\\.txt)")
names <- projects

rpc_all_smpl2 <- map_df(seq_along(file_paths), function(i) {
  read.csv(file_paths[i], sep= "\t") %>%
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
    select(c(RG, N, type, project, fraction, fraction_type_project))
}) %>%
  #filter(type != "Unmapped") %>%
  mutate(type = factor(type, levels=c("Exon", "Intron", "Intergenic", "Ambiguity", "Unmapped")))

## deML (no deML data)
# ass_reads <- read.delim(file="/data/share/htp/prime-seq_NextGen/data/FC2024_05_02_PoP96_BA/deML_summary_to_read_in.txt") %>%
#   select(RG, assigned, total) %>%
#   dplyr::rename("project"="RG") %>%
#   filter(str_detect(project, "primeseq")) %>%
#   mutate(project = str_extract(project, "...$")) %>%
#   mutate(across(project, factor, levels=c("old", "new", "PTO"))) %>%
#   mutate(non_assigned = total-assigned, 
#          fract_non_assigned = paste(round(non_assigned/total, digits=3)*100, "%")) %>%
#   pivot_longer(cols=c(2,4), names_to = "category", values_to = "reads") %>%
#   mutate(fract_non_assigned = ifelse(category == "non_assigned", fract_non_assigned, NA))

## Trimming data
trim_df2 <- map_df(projects, function(n) {
  read_delim(paste0("/data/share/htp/prime-seq_NextGen/other_methods/BRBseq_degraded_2mn/02_trimming/", n, ".txt"),
             col_names = c("category", "reads"), delim = "\t") %>%
    mutate(category = str_remove_all(category, "[:]")
      ,reads = str_remove_all(reads, "[ ,%\\)]")) %>%
    separate(reads, into = c("reads", "percentage"), sep = "\\(", convert = TRUE) %>%
    mutate(condition = n)
}) %>%
  mutate(reads = as.numeric(reads))

## Barcodes assignment
file_paths1 <- list.files(pattern = "kept_barcodes_binned\\.txt$", recursive = TRUE)

nreads_barcodes2 <- map_df(seq_along(file_paths1), function(i) {
  barcodes <- read.csv(file_paths1[i])
  barcodes %>%
    mutate(project = projects[i])
}) %>%
  bind_rows()

nreads2 <- nreads_barcodes2 %>%
  group_by(project) %>%
  summarise(barcode_assigned = sum(n))


## Read UMI data
read_umi_raw2 <- map_df(projects, function(n) {
  gene_counts <- read_rds(paste0("/data/share/htp/prime-seq_NextGen/other_methods/BRBseq_degraded_2mn/03_zUMIs/", n, "/zUMIs_output/stats/", n, ".bc.READcounts.rds")) %>%
    dplyr::rename(read_count = N, SampleID = RG)
  umi_counts <- read_delim(paste0("/data/share/htp/prime-seq_NextGen/other_methods/BRBseq_degraded_2mn/03_zUMIs/", n, "/zUMIs_output/stats/", n, ".UMIcounts.txt")) %>%
    dplyr::rename(umi_count = Count)
  
  full_join(gene_counts, umi_counts, by = c("SampleID", "type")) %>%
    mutate(project = n)
}) 
read_umi2 <- read_umi_raw2 %>%
  filter(type %in% c("Exon", "Intron") & SampleID != "bad") %>%
  mutate(umi_fraction = umi_count / read_count)

read_umi_summary2 <- read_umi2 %>%
  group_by(project) %>%
  summarize(UMI = sum(umi_count), .groups = 'drop')

# combine data
read_funnel_2  <-
    trim_df2 %>%
      select(-percentage) %>%
      pivot_wider(names_from = category, values_from = reads) %>%
      group_by(condition) %>%
      mutate(trimmed = `Total read pairs processed` - `Pairs that were too short`) %>%
      mutate(index_assigned = `Total read pairs processed`) %>%
      select(project = condition, index_assigned, trimmed) %>%
  full_join(
    nreads2,
    by="project"
  ) %>%
  pivot_longer(cols=-1, names_to = "step", values_to = "reads") %>%
  bind_rows(rpc_all_smpl2 %>% 
              select(BC = RG, N, type, project) %>%
              ungroup() %>% 
              group_by(project, BC) %>% 
              filter(type %in% c("Intron", "Exon")) %>% 
              summarise(inex=sum(N), .groups = 'drop') %>%
              full_join(read_umi2 %>% 
                          select(BC = SampleID, UMI= umi_count, type, project) %>%
                          group_by(project, BC) %>%
                          filter(type %in% c("Intron", "Exon")) %>% 
                          summarise(UMI=sum(UMI), .groups = 'drop'),
                        by=c("project", "BC")) %>%
              pivot_longer(cols=-c(BC, project), names_to = "step", values_to = "reads")
  )  %>%
  ungroup() %>%
  mutate(step = factor(step, levels = c("total", "index_assigned", "trimmed", "barcode_assigned", "inex", "UMI"))) %>%
  mutate(max_rep=max(reads), 
         reads_norm = reads/max_rep) %>% # norm to total reads per  rep
  mutate(reads_norm = reads_norm*ifelse(!is.na(BC), 4, 1))

read_funnel_avg_2 <- read_funnel_2 %>%
  ungroup() %>%
  group_by(project, step) %>%
  summarise(reads=sum(reads)) %>% # summarise per step
  ungroup() %>%
  mutate(reads_rel = reads/max(reads))

#------------------------------------------------------------------------------------------------------------------------------------------------------------
# data3 ----
## Mapping Fractions
setwd("/data/share/htp/prime-seq_NextGen/other_methods/BRBseq_degraded_1mn/03_zUMIs")
file_paths <- list.files(pattern = "\\.readspercell\\.txt$", recursive=T)
file_paths

projects <- str_extract(file_paths, "([^/]*)(?=\\.readspercell\\.txt)")
names <- projects

rpc_all_smpl3 <- map_df(seq_along(file_paths), function(i) {
  read.csv(file_paths[i], sep= "\t") %>%
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
    select(c(RG, N, type, project, fraction, fraction_type_project))
}) %>%
  #filter(type != "Unmapped") %>%
  mutate(type = factor(type, levels=c("Exon", "Intron", "Intergenic", "Ambiguity", "Unmapped")))

## deML (no deML data)
# ass_reads <- read.delim(file="/data/share/htp/prime-seq_NextGen/data/FC2024_05_02_PoP96_BA/deML_summary_to_read_in.txt") %>%
#   select(RG, assigned, total) %>%
#   dplyr::rename("project"="RG") %>%
#   filter(str_detect(project, "primeseq")) %>%
#   mutate(project = str_extract(project, "...$")) %>%
#   mutate(across(project, factor, levels=c("old", "new", "PTO"))) %>%
#   mutate(non_assigned = total-assigned, 
#          fract_non_assigned = paste(round(non_assigned/total, digits=3)*100, "%")) %>%
#   pivot_longer(cols=c(2,4), names_to = "category", values_to = "reads") %>%
#   mutate(fract_non_assigned = ifelse(category == "non_assigned", fract_non_assigned, NA))

## Trimming data
trim_df3 <- map_df(projects, function(n) {
  read_delim(paste0("/data/share/htp/prime-seq_NextGen/other_methods/BRBseq_degraded_1mn/02_trimming/", n, ".txt"),
             col_names = c("category", "reads")) %>%
    mutate(reads = str_remove_all(reads, "[ ,%\\)]")) %>%
    separate(reads, into = c("reads", "percentage"), sep = "\\(", convert = TRUE) %>%
    mutate(condition = n)
}) %>%
  mutate(reads = as.numeric(reads))

## Barcodes assignment
file_paths1 <- list.files(pattern = "kept_barcodes_binned\\.txt$", recursive = TRUE)

nreads_barcodes3 <- map_df(seq_along(file_paths1), function(i) {
  barcodes <- read.csv(file_paths1[i])
  barcodes %>%
    mutate(project = projects[i])
}) %>%
  bind_rows()

nreads3 <- nreads_barcodes3 %>%
  group_by(project) %>%
  summarise(barcode_assigned = sum(n))


## Read UMI data
read_umi_raw3 <- map_df(projects, function(n) {
  gene_counts <- read_rds(paste0("/data/share/htp/prime-seq_NextGen/other_methods/BRBseq_degraded_1mn/03_zUMIs/", n, "/zUMIs_output/stats/", n, ".bc.READcounts.rds")) %>%
    dplyr::rename(read_count = N, SampleID = RG)
  umi_counts <- read_delim(paste0("/data/share/htp/prime-seq_NextGen/other_methods/BRBseq_degraded_1mn/03_zUMIs/", n, "/zUMIs_output/stats/", n, ".UMIcounts.txt")) %>%
    dplyr::rename(umi_count = Count)
  
  full_join(gene_counts, umi_counts, by = c("SampleID", "type")) %>%
    mutate(project = n)
}) 
read_umi3 <- read_umi_raw3 %>%
  filter(type %in% c("Exon", "Intron") & SampleID != "bad") %>%
  mutate(umi_fraction = umi_count / read_count)

read_umi_summary3 <- read_umi3 %>%
  group_by(project) %>%
  summarize(UMI = sum(umi_count), .groups = 'drop')

# combine data
read_funnel_3  <-
  trim_df3 %>%
  select(-percentage) %>%
  pivot_wider(names_from = category, values_from = reads) %>%
  group_by(condition) %>%
  mutate(trimmed = `Total read pairs processed` - `Pairs that were too short`) %>%
  mutate(index_assigned = `Total read pairs processed`) %>%
  select(project = condition, index_assigned, trimmed) %>%
  full_join(
    nreads3,
    by="project"
  ) %>%
  pivot_longer(cols=-1, names_to = "step", values_to = "reads") %>%
  bind_rows(rpc_all_smpl3 %>% 
              select(BC = RG, N, type, project) %>%
              ungroup() %>% 
              group_by(project, BC) %>% 
              filter(type %in% c("Intron", "Exon")) %>% 
              summarise(inex=sum(N), .groups = 'drop') %>%
              full_join(read_umi3 %>% 
                          select(BC = SampleID, UMI= umi_count, type, project) %>%
                          group_by(project, BC) %>%
                          filter(type %in% c("Intron", "Exon")) %>% 
                          summarise(UMI=sum(UMI), .groups = 'drop'),
                        by=c("project", "BC")) %>%
              pivot_longer(cols=-c(BC, project), names_to = "step", values_to = "reads")
  )  %>%
  ungroup() %>%
  mutate(step = factor(step, levels = c("total", "index_assigned", "trimmed", "barcode_assigned", "inex", "UMI"))) %>%
  mutate(max_rep=max(reads), 
         reads_norm = reads/max_rep) %>% # norm to total reads per  rep
  mutate(reads_norm = reads_norm*ifelse(!is.na(BC), 4, 1))

read_funnel_avg_3 <- read_funnel_3 %>%
  ungroup() %>%
  group_by(project, step) %>%
  summarise(reads=sum(reads)) %>% # summarise per step
  ungroup() %>%
  mutate(reads_rel = reads/max(reads))



#------------------------------------------------------------------------------------------------------------------------------------------------------------
# data4 ----
## Mapping Fractions
setwd("/data/share/htp/prime-seq_NextGen/other_methods/BRBseq_degraded_nt/03_zUMIs")
file_paths <- list.files(pattern = "\\.readspercell\\.txt$", recursive=T)
file_paths

projects <- str_extract(file_paths, "([^/]*)(?=\\.readspercell\\.txt)")
names <- projects

rpc_all_smpl4 <- map_df(seq_along(file_paths), function(i) {
  read.csv(file_paths[i], sep= "\t") %>%
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
    select(c(RG, N, type, project, fraction, fraction_type_project))
}) %>%
  #filter(type != "Unmapped") %>%
  mutate(type = factor(type, levels=c("Exon", "Intron", "Intergenic", "Ambiguity", "Unmapped")))


## Trimming data
trim_df4 <- map_df(projects, function(n) {
  read_delim(paste0("/data/share/htp/prime-seq_NextGen/other_methods/BRBseq_degraded_nt/02_trimming/", n, ".txt"),
             col_names = c("category", "reads")) %>%
    mutate(category = str_remove_all(category, "[:]")
           ,reads = str_remove_all(reads, "[ ,%\\)]")) %>%
    separate(reads, into = c("reads", "percentage"), sep = "\\(", convert = TRUE) %>%
    mutate(condition = n)
}) %>%
  mutate(reads = as.numeric(reads))

## Barcodes assignment
file_paths1 <- list.files(pattern = "kept_barcodes_binned\\.txt$", recursive = TRUE)

nreads_barcodes4 <- map_df(seq_along(file_paths1), function(i) {
  barcodes <- read.csv(file_paths1[i])
  barcodes %>%
    mutate(project = projects[i])
}) %>%
  bind_rows()

nreads4 <- nreads_barcodes4 %>%
  group_by(project) %>%
  summarise(barcode_assigned = sum(n))


## Read UMI data
read_umi_raw4 <- map_df(projects, function(n) {
  gene_counts <- read_rds(paste0("/data/share/htp/prime-seq_NextGen/other_methods/BRBseq_degraded_nt/03_zUMIs/", n, "/zUMIs_output/stats/", n, ".bc.READcounts.rds")) %>%
    dplyr::rename(read_count = N, SampleID = RG)
  umi_counts <- read_delim(paste0("/data/share/htp/prime-seq_NextGen/other_methods/BRBseq_degraded_nt/03_zUMIs/", n, "/zUMIs_output/stats/", n, ".UMIcounts.txt")) %>%
    dplyr::rename(umi_count = Count)
  
  full_join(gene_counts, umi_counts, by = c("SampleID", "type")) %>%
    mutate(project = n)
}) 
read_umi4 <- read_umi_raw4 %>%
  filter(type %in% c("Exon", "Intron") & SampleID != "bad") %>%
  mutate(umi_fraction = umi_count / read_count)

read_umi_summary4 <- read_umi4 %>%
  group_by(project) %>%
  summarize(UMI = sum(umi_count), .groups = 'drop')

# combine data
read_funnel_4  <-
  trim_df4 %>%
  select(-percentage) %>%
  pivot_wider(names_from = category, values_from = reads) %>%
  group_by(condition) %>%
  mutate(trimmed = `Total read pairs processed` - `Pairs that were too short`) %>%
  mutate(index_assigned = `Total read pairs processed`) %>%
  select(project = condition, index_assigned, trimmed) %>%
  full_join(
    nreads4,
    by="project"
  ) %>%
  pivot_longer(cols=-1, names_to = "step", values_to = "reads") %>%
  bind_rows(rpc_all_smpl4 %>% 
              select(BC = RG, N, type, project) %>%
              ungroup() %>% 
              group_by(project, BC) %>% 
              filter(type %in% c("Intron", "Exon")) %>% 
              summarise(inex=sum(N), .groups = 'drop') %>%
              full_join(read_umi4 %>% 
                          select(BC = SampleID, UMI= umi_count, type, project) %>%
                          group_by(project, BC) %>%
                          filter(type %in% c("Intron", "Exon")) %>% 
                          summarise(UMI=sum(UMI), .groups = 'drop'),
                        by=c("project", "BC")) %>%
              pivot_longer(cols=-c(BC, project), names_to = "step", values_to = "reads")
  )  %>%
  ungroup() %>%
  mutate(step = factor(step, levels = c("total", "index_assigned", "trimmed", "barcode_assigned", "inex", "UMI"))) %>%
  mutate(max_rep=max(reads), 
         reads_norm = reads/max_rep) %>% # norm to total reads per  rep
  mutate(reads_norm = reads_norm*ifelse(!is.na(BC), 4, 1))

read_funnel_avg_4 <- read_funnel_4 %>%
  ungroup() %>%
  group_by(project, step) %>%
  summarise(reads=sum(reads)) %>% # summarise per step
  ungroup() %>%
  mutate(reads_rel = reads/max(reads))

#---------------------------------------------------------------------------------------------------------------------------------
read_funnel_BRB <- 
  rbind(read_funnel_1,
        read_funnel_2,
        read_funnel_3,
        read_funnel_4)

read_funnel_BRB_avg <- 
  rbind(read_funnel_avg_1,
        read_funnel_avg_2,
        read_funnel_avg_3,
        read_funnel_avg_4)



BRB_funnel <- ggplot()+
  # stat_summary(data=read_funnel_BRB, 
  #              aes(y=step, x=reads_norm, color=project), 
  #              fun.data = "mean_sdl",
  #              fun.args = list(mult = 1),
  #              geom="errorbar", 
  #              position=position_dodge(0.2), 
  #              width=.2,
  #              size=.5) + 
  geom_point(data=read_funnel_BRB, aes(y=step, x=reads_norm, color=project, group=project), alpha=0.7) + 
  geom_line(data=read_funnel_BRB_avg, aes(y=step, x=reads_rel, color=project, group=project)) + 
  scale_x_continuous(position = "top") + 
  ylab("Processing Steps \n \u27F5") +
  xlab("Relative Reads")+
  scale_y_discrete(limits = rev(levels(read_funnel_BRB$step)))+
  labs(color="Protocol")+
  theme_pubr(legend = "right")

BRB_funnel

ggsave("/data/share/htp/prime-seq_NextGen/figures/fig_SX_BRB_funnel.pdf",
       BRB_funnel,
       width=8,
       height=4,
       device = cairo_pdf)
