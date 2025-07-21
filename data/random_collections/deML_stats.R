library(dplyr)
library(stringr)
library(ggplot2)

stats <- read.delim("/data/share/htp/prime-seq_NextGen/data/random_collections/deML_stats_updatedEB.txt", skip = 1)

stats$fraction.exp.total.reads <- gsub(",", ".", stats$fraction.exp.total.reads) %>% as.numeric
stats$assigned. <- gsub(",", ".", stats$assigned.)

stats <- stats %>%
  mutate(protocol_simple = factor(ifelse(str_detect(protocol, "old"), "old", 
                                         ifelse(str_detect(protocol, "PTO"), "PTO", protocol)),
                                  levels = c("old", "new", "PTO", "NA", "hightemp")),
         assigned. = str_replace(assigned., "%", "") %>% as.numeric,
         unknown. = str_replace(unknown., "%", "") %>% as.numeric) %>% 
  filter(FC!="23_01_01") %>% 
  filter(!(FC=="24_05_02" & protocol=="old")) %>%
  filter(!(FC=="24_03_01" & protocol=="old")) %>%
  filter(!(FC=="24_05_01" & protocol=="old"))

filtered_stats <- stats %>% 
  filter(protocol_simple %in% c("old", "new", "PTO") & !is.na(protocol_simple))

plot_readfracs <- filtered_stats %>% mutate(FC=FC, reads=total.reads.FC)


# Read the data and perform all transformations in one go
filtered_stats <- read.delim("/data/share/htp/prime-seq_NextGen/data/random_collections/deML_stats_updatedEB.txt", skip = 1) %>%
  mutate(fraction.exp.total.reads = as.numeric(gsub(",", ".", fraction.exp.total.reads)),
         assigned. = as.numeric(str_replace(gsub(",", ".", assigned.), "%", "")),
         unknown. = as.numeric(str_replace(unknown., "%", "")),
         protocol_simple = factor(ifelse(str_detect(protocol, "old"), "old", 
                                         ifelse(str_detect(protocol, "PTO"), "PTO", protocol)),
                                  levels = c("old", "new", "PTO", "NA", "hightemp"))) %>%
  filter(FC != "23_01_01" & !(FC == "24_05_02" & protocol == "old") & 
           !(FC == "24_03_01" & protocol == "old") & !(FC == "24_05_01" & protocol == "old")) %>%
  filter(protocol_simple %in% c("old", "new", "PTO") & !is.na(protocol_simple))



#### Fraction of expected flow cell reads:
ggplot(filtered_stats[filtered_stats$protocol=='PTO'|filtered_stats$protocol=='old'|filtered_stats$protocol=='old -RNA',],
       aes(x=protocol_simple, y=fraction.exp.total.reads*100, fill=protocol_simple)) +
  #geom_boxplot() +
  geom_violin() +
  #coord_flip()+
  ggbeeswarm::geom_beeswarm(aes(color=FC.type)) +
  scale_color_manual(values = c("forestgreen", "black", "grey40"))

ggplot(filtered_stats, aes(x=protocol_simple ,y=assigned., fill=protocol_simple)) +
  geom_violin() +
  ggbeeswarm::geom_beeswarm() +
  scale_y_log10()

ggplot(filtered_stats, aes(x=protocol_simple ,y=unknown., fill=protocol_simple)) +
  geom_violin() +
  ggbeeswarm::geom_beeswarm() +
  scale_y_log10()

#### whats all in there

samples_table <- filtered_stats %>% select(c("species", "tissue_or_cell_type", "cell_type_simple")) %>% subset(tissue_or_cell_type != "") %>% distinct %>% arrange(species)

write.csv(samples_table, "/data/share/htp/prime-seq_NextGen/data/random_collections/samples_table.csv")

colnames(filtered_stats)


ggplot(filtered_stats[filtered_stats$protocol=='PTO'|filtered_stats$protocol=='old'|filtered_stats$protocol=='old -RNA',],
       aes(x=protocol_simple, y=fraction.exp.total.reads*100, fill=protocol_simple)) +
  #geom_boxplot() +
  geom_violin() +
  #coord_flip()+
  ggbeeswarm::geom_beeswarm(aes(color=cell_type_simple))

