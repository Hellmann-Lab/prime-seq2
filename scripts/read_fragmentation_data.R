file1 <- here("data/Fragmentation/2100 expert_High Sensitivity DNA Assay_DE72901629_2023-05-10_09-52-27.xml")
file2 <- here("data/Fragmentation/2100 expert_High Sensitivity DNA Assay_DE72901629_2023-05-11_12-10-32.xml")
file3 <- here("data/Fragmentation/2100 expert_High Sensitivity DNA Assay_DE72901629_2023-06-01_16-52-38.xml")

data1 <- read.bioanalyzer(file1)
data2 <- read.bioanalyzer(file2)
data3 <- read.bioanalyzer(file3)


data_lib <- subset(rbind(subset(data1, 
                                !startsWith(as.character(sample.name), "brain_5min_80ng")
), 
data2, 
data3
), 
startsWith(as.character(sample.name), "brain") | startsWith(as.character(sample.name), "Felix Frag")
)

df <- region.ratio(
  data_lib,
  sum.variable = "concentration",
  c(100, 8000),
  c(250, 650)
) %>% 
  as.data.frame() %>%
  dplyr::rename(fract = 'concentration ratio in length 250-650/100-8000') %>%
  mutate(sample.name = data_lib[["samples"]][,3])%>%
  mutate(
    time = as.numeric(sub(":30", ".5", sub("min", "", str_split(sample.name, "_", simplify = TRUE)[, 2]))),
    amount = as.numeric(sub("ng", "", str_split(sample.name, "_", simplify = TRUE)[, 3])),
    rep = str_split(sample.name, "_", simplify = TRUE)[, 4]) %>%
  arrange(amount, time, rep)

df_mean <- df %>%
  group_by(time, amount) %>%
  summarise(fract_mean = mean(fract))