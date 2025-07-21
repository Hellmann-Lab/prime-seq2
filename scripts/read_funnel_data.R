# rep 1 ----
# read in data

## Mapping Fractions
file_paths <- list.files(path = here("data/FC2024_05_02_PoP96_BA/03_zUMIs"),
                         pattern = "\\.readspercell\\.txt$", recursive=T)
file_paths

projects <- str_extract(file_paths, "(...)(?=\\.readspercell\\.txt)")
names <- projects

rpc_all_smpl <- map_df(seq_along(file_paths), function(i) {
  read.csv(here("data/FC2024_05_02_PoP96_BA/03_zUMIs", file_paths[i]), sep= "\t") %>%
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

## deML
ass_reads <- read.delim(file=here("data/FC2024_05_02_PoP96_BA/deML_summary_to_read_in.txt")) %>%
  dplyr::select(RG, assigned, total) %>%
  dplyr::rename("project"="RG") %>%
  filter(str_detect(project, "primeseq")) %>%
  mutate(project = str_extract(project, "...$")) %>%
  mutate(across(project, factor, levels=c("old", "new", "PTO"))) %>%
  mutate(non_assigned = total-assigned, 
         fract_non_assigned = paste(round(non_assigned/total, digits=3)*100, "%")) %>%
  pivot_longer(cols=c(2,4), names_to = "category", values_to = "reads") %>%
  mutate(fract_non_assigned = ifelse(category == "non_assigned", fract_non_assigned, NA))

## Trimming data
trim_df <- map_df(projects, function(n) {
  read_delim(here(paste0("data/FC2024_05_02_PoP96_BA/02_trimming/", n, ".txt")),
             col_names = c("category", "reads")) %>%
    mutate(reads = str_remove_all(reads, "[ ,%\\)]")) %>%
    separate(reads, into = c("reads", "percentage"), sep = "\\(", convert = TRUE) %>%
    mutate(condition = n)
}) %>%
  mutate(reads = as.numeric(reads))

## Barcodes assignment
file_paths1 <- list.files(path = here("data/FC2024_05_02_PoP96_BA"),
                          pattern = "kept_barcodes_binned\\.txt$", recursive = TRUE)

nreads_barcodes <- map_df(seq_along(file_paths1), function(i) {
  barcodes <- read.csv(here("data/FC2024_05_02_PoP96_BA", file_paths1[i]))
  barcodes %>%
    mutate(project = projects[i])
}) %>%
  bind_rows()

nreads <- nreads_barcodes %>%
  group_by(project) %>%
  summarise('Barcode assigned' = sum(n))


## Read UMI data
read_umi_raw <- map_df(projects, function(n) {
  gene_counts <- read_rds(here(paste0("data/FC2024_05_02_PoP96_BA/03_zUMIs/", n, "/zUMIs_output/stats/OldNewPTO_", n, ".bc.READcounts.rds"))) %>%
    dplyr::rename(read_count = N, SampleID = RG)
  umi_counts <- read_delim(here(paste0("data/FC2024_05_02_PoP96_BA/03_zUMIs/", n, "/zUMIs_output/stats/OldNewPTO_", n, ".UMIcounts.txt"))) %>%
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
read_funnel1  <-
  ass_reads %>%
  filter(category == "assigned") %>%
  dplyr::select(project,  'Total Reads' = total, 'Index assigned' = reads) %>%
  full_join(
    trim_df %>%
      dplyr::select(-percentage) %>%
      tidyr::pivot_wider(names_from = category, values_from = reads) %>%
      group_by(condition) %>%
      mutate('Trimmed & Filtered' = `Total read pairs processed` - `Pairs that were too short`) %>%
      dplyr::select(project = condition, 'Trimmed & Filtered'), 
    by="project"
  ) %>%
  full_join(
    nreads,
    by="project"
  ) %>%
  pivot_longer(cols=-1, names_to = "step", values_to = "reads") %>%
  bind_rows(rpc_all_smpl %>% 
              dplyr::select(BC = RG, N, type, project) %>%
              ungroup() %>% 
              group_by(project, BC) %>% 
              filter(type %in% c("Intron", "Exon")) %>% 
              summarise('Exon & Intron selected'=sum(N), .groups = 'drop') %>%
              full_join(read_umi %>% 
                          dplyr::select(BC = SampleID, UMI= umi_count, type, project) %>%
                          group_by(project, BC) %>%
                          filter(type %in% c("Intron", "Exon")) %>% 
                          summarise('UMI collapsed'=sum(UMI), .groups = 'drop'),
                        by=c("project", "BC")) %>%
              pivot_longer(cols=-c(BC, project), names_to = "step", values_to = "reads")
  )  %>%
  ungroup()

# rep 2 ----
# read in data

## Mapping Fractions

file_paths <- list.files(path = here("data/FC2024_06_01_PoP96_FP_BA/03_zUMIs"),
                         pattern = "\\.readspercell\\.txt$", recursive=T)
file_paths

projects <- str_extract(file_paths, "(...)(?=\\.readspercell\\.txt)")
names <- projects

rpc_all_smpl2 <- map_df(seq_along(file_paths), function(i) {
  read.csv(here("data/FC2024_06_01_PoP96_FP_BA/03_zUMIs", file_paths[i]), sep= "\t") %>%
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

## deML
ass_reads2 <- read.delim(file=here("data/FC2024_06_01_PoP96_FP_BA/deML_summary_to_read_in.txt")) %>%
  dplyr::select(RG, assigned, total) %>%
  dplyr::rename("project"="RG") %>%
  filter(str_detect(project, "primeseq")) %>%
  mutate(project = str_extract(project, "...$")) %>%
  mutate(across(project, factor, levels=c("old", "new", "PTO"))) %>%
  mutate(non_assigned = total-assigned, 
         fract_non_assigned = paste(round(non_assigned/total, digits=3)*100, "%")) %>%
  pivot_longer(cols=c(2,4), names_to = "category", values_to = "reads") %>%
  mutate(fract_non_assigned = ifelse(category == "non_assigned", fract_non_assigned, NA))

## Trimming data
trim_df2 <- map_df(projects, function(n) {
  read_delim(here(paste0("data/FC2024_06_01_PoP96_FP_BA/02_trimming/", n, ".txt")),
             col_names = c("category", "reads")) %>%
    mutate(reads = str_remove_all(reads, "[ ,%\\)]")) %>%
    separate(reads, into = c("reads", "percentage"), sep = "\\(", convert = TRUE) %>%
    mutate(condition = n)
}) %>%
  mutate(reads = as.numeric(reads))

## Barcodes assignment
file_paths1 <- list.files(path = here("data/FC2024_06_01_PoP96_FP_BA"),
                          pattern = "kept_barcodes_binned\\.txt$", recursive = TRUE)

nreads_barcodes2 <- map_df(seq_along(file_paths1), function(i) {
  barcodes <- read.csv(here("data/FC2024_06_01_PoP96_FP_BA", file_paths1[i]))
  barcodes %>%
    mutate(project = projects[i])
}) %>%
  bind_rows()

nreads2 <- nreads_barcodes2 %>%
  group_by(project) %>%
  summarise('Barcode assigned' = sum(n))


## Read UMI data
read_umi_raw2 <- map_df(projects, function(n) {
  gene_counts <- read_rds(here(paste0("data/FC2024_06_01_PoP96_FP_BA/03_zUMIs/", n, "/zUMIs_output/stats/OldNewPTO_", n, ".bc.READcounts.rds"))) %>%
    dplyr::rename(read_count = N, SampleID = RG)
  umi_counts <- read_delim(here(paste0("data/FC2024_06_01_PoP96_FP_BA/03_zUMIs/", n, "/zUMIs_output/stats/OldNewPTO_", n, ".UMIcounts.txt"))) %>%
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
read_funnel2  <-
  ass_reads2 %>%
  filter(category == "assigned") %>%
  dplyr::select(project,  'Total Reads' = total, 'Index assigned' = reads) %>%
  full_join(
    trim_df2 %>%
      dplyr::select(-percentage) %>%
      pivot_wider(names_from = category, values_from = reads) %>%
      group_by(condition) %>%
      mutate('Trimmed & Filtered' = `Total read pairs processed` - `Pairs that were too short`) %>%
      dplyr::select(project = condition, 'Trimmed & Filtered'), 
    by="project"
  ) %>%
  full_join(
    nreads2,
    by="project"
  ) %>%
  pivot_longer(cols=-1, names_to = "step", values_to = "reads") %>%
  bind_rows(rpc_all_smpl2 %>% 
              dplyr::select(BC = RG, N, type, project) %>%
              ungroup() %>% 
              group_by(project, BC) %>% 
              filter(type %in% c("Intron", "Exon")) %>% 
              summarise('Exon & Intron selected'=sum(N), .groups = 'drop') %>%
              full_join(read_umi2 %>% 
                          dplyr::select(BC = SampleID, UMI= umi_count, type, project) %>%
                          group_by(project, BC) %>%
                          filter(type %in% c("Intron", "Exon")) %>% 
                          summarise('UMI collapsed'=sum(UMI), .groups = 'drop'),
                        by=c("project", "BC")) %>%
              pivot_longer(cols=-c(BC, project), names_to = "step", values_to = "reads")
  )  %>%
  ungroup()

# rep 3 ----
# read in data

## Mapping Fractions
file_paths <- list.files(path = here("data/FC2024_08_01_PoP96III/03_zUMIs"),
                         pattern = "\\.readspercell\\.txt$", recursive=T)
file_paths

projects <- str_extract(file_paths, "(...)(?=\\.readspercell\\.txt)")
names <- projects

rpc_all_smpl3 <- map_df(seq_along(file_paths), function(i) {
  read.csv(here("data/FC2024_08_01_PoP96III/03_zUMIs", file_paths[i]), sep= "\t") %>%
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

## deML
ass_reads3 <- read.delim(file=here("data/FC2024_08_01_PoP96III/deML_summary_to_read_in.txt")) %>%
  dplyr::select(RG, assigned, total) %>%
  dplyr::rename("project"="RG") %>%
  filter(str_detect(project, "primeseq_PoP")) %>%
  mutate(project = str_extract(project, "...$")) %>%
  mutate(across(project, factor, levels=c("old", "new", "PTO"))) %>%
  mutate(non_assigned = total-assigned, 
         fract_non_assigned = paste(round(non_assigned/total, digits=3)*100, "%")) %>%
  pivot_longer(cols=c(2,4), names_to = "category", values_to = "reads") %>%
  mutate(fract_non_assigned = ifelse(category == "non_assigned", fract_non_assigned, NA))

## Trimming data
trim_df3 <- map_df(projects, function(n) {
  read_delim(here(paste0("data/FC2024_08_01_PoP96III/02_trimming/", n, ".txt")),
             col_names = c("category", "reads")) %>%
    mutate(reads = str_remove_all(reads, "[ ,%\\)]")) %>%
    separate(reads, into = c("reads", "percentage"), sep = "\\(", convert = TRUE) %>%
    mutate(condition = n)
}) %>%
  mutate(reads = as.numeric(reads))

## Barcodes assignment
file_paths1 <- list.files(path = here("data/FC2024_08_01_PoP96III"),
                          pattern = "kept_barcodes_binned\\.txt$", recursive = TRUE)

nreads_barcodes3 <- map_df(seq_along(file_paths1), function(i) {
  barcodes <- read.csv(here("data/FC2024_08_01_PoP96III", file_paths1[i]))
  barcodes %>%
    mutate(project = projects[i])
}) %>%
  bind_rows()

nreads3 <- nreads_barcodes3 %>%
  group_by(project) %>%
  summarise('Barcode assigned' = sum(n))


## Read UMI data
read_umi_raw3 <- map_df(projects, function(n) {
  gene_counts <- read_rds(here(paste0("data/FC2024_08_01_PoP96III/03_zUMIs/", n, "/zUMIs_output/stats/OldNewPTO_", n, ".bc.READcounts.rds"))) %>%
    dplyr::rename(read_count = N, SampleID = RG)
  umi_counts <- read_delim(here(paste0("data/FC2024_08_01_PoP96III/03_zUMIs/", n, "/zUMIs_output/stats/OldNewPTO_", n, ".UMIcounts.txt"))) %>%
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
read_funnel3  <-
  ass_reads3 %>%
  filter(category == "assigned") %>%
  dplyr::select(project, 'Total Reads' = total, 'Index assigned' = reads) %>%
  full_join(
    trim_df3 %>%
      dplyr::select(-percentage) %>%
      pivot_wider(names_from = category, values_from = reads) %>%
      group_by(condition) %>%
      mutate('Trimmed & Filtered' = `Total read pairs processed` - `Pairs that were too short`) %>%
      dplyr::select(project = condition, 'Trimmed & Filtered'), 
    by="project"
  ) %>%
  full_join(
    nreads3,
    by="project"
  ) %>%
  pivot_longer(cols=-1, names_to = "step", values_to = "reads") %>%
  bind_rows(rpc_all_smpl3 %>% 
              dplyr::select(BC = RG, N, type, project) %>%
              ungroup() %>% 
              group_by(project, BC) %>% 
              filter(type %in% c("Intron", "Exon")) %>% 
              summarise('Exon & Intron selected'=sum(N), .groups = 'drop') %>%
              full_join(read_umi3 %>% 
                          dplyr::select(BC = SampleID, UMI= umi_count, type, project) %>%
                          group_by(project, BC) %>%
                          filter(type %in% c("Intron", "Exon")) %>% 
                          summarise('UMI collapsed'=sum(UMI), .groups = 'drop'),
                        by=c("project", "BC")) %>%
              pivot_longer(cols=-c(BC, project), names_to = "step", values_to = "reads")
  )  %>%
  ungroup()

# rep 4 ----
# read in data

## Mapping Fractions
file_paths <- list.files(path = here("data/FC2024_11_01_PoP64_1/03_zUMIs"),
                         pattern = "\\.readspercell\\.txt$", recursive=T)
file_paths

projects <- str_extract(file_paths, "(...)(?=\\.readspercell\\.txt)")
names <- projects

rpc_all_smpl4 <- map_df(seq_along(file_paths), function(i) {
  read.csv(here("data/FC2024_11_01_PoP64_1/03_zUMIs", file_paths[i]), sep= "\t") %>%
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

## deML
ass_reads4 <- read.delim(file=here("data/FC2024_11_01_PoP64_1/deML_summary_to_read_in.txt")) %>%
  dplyr::select(RG, assigned, total) %>%
  dplyr::rename("project"="RG") %>%
  filter(str_detect(project, "primeseq_PoP")) %>%
  mutate(project = str_replace(project, "_corrected", "")) %>%
  mutate(project = str_extract(project, "...$")) %>%
  mutate(across(project, factor, levels=c("old", "new", "PTO"))) %>%
  mutate(non_assigned = total-assigned, 
         fract_non_assigned = paste(round(non_assigned/total, digits=3)*100, "%")) %>%
  pivot_longer(cols=c(2,4), names_to = "category", values_to = "reads") %>%
  mutate(fract_non_assigned = ifelse(category == "non_assigned", fract_non_assigned, NA))

## Trimming data
trim_df4 <- map_df(projects, function(n) {
  read_delim(here(paste0("data/FC2024_11_01_PoP64_1/02_trimming/", n, ".txt")),
             col_names = c("category", "reads")) %>%
    mutate(reads = str_remove_all(reads, "[ ,%\\)]")) %>%
    separate(reads, into = c("reads", "percentage"), sep = "\\(", convert = TRUE) %>%
    mutate(condition = n)
}) %>%
  mutate(reads = as.numeric(reads))

## Barcodes assignment
file_paths1 <- list.files(path = here("data/FC2024_11_01_PoP64_1"),
                          pattern = "kept_barcodes_binned\\.txt$", recursive = TRUE)

nreads_barcodes4 <- map_df(seq_along(file_paths1), function(i) {
  barcodes <- read.csv(here("data/FC2024_11_01_PoP64_1", file_paths1[i]))
  barcodes %>%
    mutate(project = projects[i])
}) %>%
  bind_rows()

nreads4 <- nreads_barcodes4 %>%
  group_by(project) %>%
  summarise('Barcode assigned' = sum(n))


## Read UMI data
read_umi_raw4 <- map_df(projects, function(n) {
  gene_counts <- read_rds(here(paste0("data/FC2024_11_01_PoP64_1/03_zUMIs/", n, "/zUMIs_output/stats/PoP64_", n, ".bc.READcounts.rds"))) %>%
    dplyr::rename(read_count = N, SampleID = RG)
  umi_counts <- read_delim(here(paste0("data/FC2024_11_01_PoP64_1/03_zUMIs/", n, "/zUMIs_output/stats/PoP64_", n, ".UMIcounts.txt"))) %>%
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
read_funnel4  <-
  ass_reads4 %>%
  filter(category == "assigned") %>%
  dplyr::select(project, 'Total Reads' = total, 'Index assigned' = reads) %>%
  full_join(
    trim_df4 %>%
      dplyr::select(-percentage) %>%
      pivot_wider(names_from = category, values_from = reads) %>%
      group_by(condition) %>%
      mutate('Trimmed & Filtered' = `Total read pairs processed` - `Pairs that were too short`) %>%
      dplyr::select(project = condition, 'Trimmed & Filtered'), 
    by="project"
  ) %>%
  full_join(
    nreads4,
    by="project"
  ) %>%
  pivot_longer(cols=-1, names_to = "step", values_to = "reads") %>%
  bind_rows(rpc_all_smpl4 %>% 
              dplyr::select(BC = RG, N, type, project) %>%
              ungroup() %>% 
              group_by(project, BC) %>% 
              filter(type %in% c("Intron", "Exon")) %>% 
              summarise('Exon & Intron selected'=sum(N), .groups = 'drop') %>%
              full_join(read_umi4 %>% 
                          dplyr::select(BC = SampleID, UMI= umi_count, type, project) %>%
                          group_by(project, BC) %>%
                          filter(type %in% c("Intron", "Exon")) %>% 
                          summarise('UMI collapsed'=sum(UMI), .groups = 'drop'),
                        by=c("project", "BC")) %>%
              pivot_longer(cols=-c(BC, project), names_to = "step", values_to = "reads")
  )  %>%
  ungroup()

# combine ----
read_funnel_all <- rbind(read_funnel1 %>% mutate(rep = 1), 
                         read_funnel2 %>% mutate(rep = 2), 
                         read_funnel3 %>% mutate(rep = 3), 
                         read_funnel4 %>% mutate(rep = 4)) %>%
  filter(project != "new") %>%
  filter(!(BC %in% unwanted_BCs)) %>% # remove 3 BCs per rep that were only used once
  mutate(project = ifelse(project == "PTO", "improved", "original"),
         step = factor(step, levels = c("Total Reads", "Index assigned", "Trimmed & Filtered", "Barcode assigned", "Exon & Intron selected", "UMI collapsed"))) %>% # bind
  group_by(rep) %>%
  mutate(max_rep=max(reads), 
         reads_norm = reads/max_rep) %>% # norm to total reads per  rep
  mutate(reads_norm = reads_norm*ifelse(!is.na(BC), 8, 1)) %>% # mutliply by BC number
  mutate(project = factor(project, levels = c("original", "improved")))

# calculate average
read_funnel_avg_all <- read_funnel_all %>%
  ungroup() %>%
  group_by(project, step, rep) %>%
  summarise(reads=sum(reads)) %>%
  group_by(rep) %>%
  mutate(reads_rel_rep = reads/max(reads)) %>%
  group_by(project, step) %>%
  summarise(reads_rel = mean(reads_rel_rep))
