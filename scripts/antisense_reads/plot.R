# Load required libraries
library(tibble)
library(here)
library(dplyr)
library(purrr)
library(readr)
library(ggplot2)
library(ggpubr)
library(tidyr)

# Set working directory
setwd(here("scripts/antisense_reads"))

# Read the antisense intergenic counts from the main script
antisense_data <- read.delim("counts_antisense_intergenics_all_projects.txt") %>%
  dplyr::rename(Anti_sense_reads = Count)

# Define project names and their corresponding reads per cell file paths
project_files <- list(
  "PoP64_PTO" = here("data/FC2024_11_01_PoP64_1/03_zUMIs/PTO/zUMIs_output/stats/PoP64_PTO.readspercell.txt"),
  "PoP64_old" = here("data/FC2024_11_01_PoP64_1/03_zUMIs/old/zUMIs_output/stats/PoP64_old.readspercell.txt"),
  "PoP96_3_PTO" = here("data/FC2024_08_01_PoP96III/03_zUMIs/PTO/zUMIs_output/stats/OldNewPTO_PTO.readspercell.txt"),
  "PoP96_3_old" = here("data/FC2024_08_01_PoP96III/03_zUMIs/old/zUMIs_output/stats/OldNewPTO_old.readspercell.txt"),
  "PoP96_2_PTO" = here("data/FC2024_06_01_PoP96_FP_BA/03_zUMIs/PTO/zUMIs_output/stats/OldNewPTO_PTO.readspercell.txt"),
  "PoP96_2_old" = here("data/FC2024_06_01_PoP96_FP_BA/03_zUMIs/old/zUMIs_output/stats/OldNewPTO_old.readspercell.txt"),
  "PoP96_1_PTO" = here("data/FC2024_05_02_PoP96_BA/03_zUMIs/PTO/zUMIs_output/stats/OldNewPTO_PTO.readspercell.txt"),
  "PoP96_1_old" = here("data/FC2024_05_02_PoP96_BA/03_zUMIs/old/zUMIs_output/stats/OldNewPTO_old.readspercell.txt")
)

# Define project names in the order they should appear in the plot
project_names <- names(project_files)

# Read intergenic data for each project
intergenics <- map_df(names(project_files), function(project) {
  file_path <- project_files[[project]]
  if (file.exists(file_path)) {
    read_delim(file_path) %>%
      filter(type == "Intergenic") %>%
      mutate(project = project) %>%
      summarise(Intergenic_reads = sum(N), .by = c(project))
  } else {
    warning(paste("File not found for project:", project, "-", file_path))
    tibble(project = project, Intergenic_reads = 0)
  }
})

# Read all reads data for each project
all <- map_df(names(project_files), function(project) {
  file_path <- project_files[[project]]
  if (file.exists(file_path)) {
    read_delim(file_path) %>%
      mutate(project = project) %>%
      summarise(All_reads = sum(N), .by = c(project))
  } else {
    warning(paste("File not found for project:", project, "-", file_path))
    tibble(project = project, All_reads = 0)
  }
})

# Define colors for the plot
colors_projects <- c("Intergenic: anti-sense" = "#FF6F61",
                     "Intergenic: other" = "#6B5B95")

# Create the plot
projects_intergenics_plot <- antisense_data %>%
  full_join(intergenics, by = c("Project" = "project")) %>%
  full_join(all, by = c("Project" = "project")) %>%
  mutate('Intergenic: anti-sense' = Anti_sense_reads / All_reads,
         'Intergenic: other' = (Intergenic_reads - Anti_sense_reads) / All_reads) %>%
  pivot_longer(cols = 5:6, names_to = "category", values_to = "fraction") %>%
  # mutate(Project = factor(Project, levels = project_names)) %>%
  ggplot(aes(y = fraction, x = Project, fill = category)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.2f", fraction)), position = position_stack(vjust = 0.5), 
            size = 2.5, color = "white", fontface = "bold") +
  theme_pubr(legend = "bottom") +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,    # right-justify so labels sit neatly under the ticks
    vjust = 1)) +
  labs(y = "Read fraction", x = "", fill = "") +
  scale_fill_manual(values = colors_projects) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# Display the plot
projects_intergenics_plot

# Save the plot
ggsave("antisense_intergenics_plot.pdf", projects_intergenics_plot, width = 10, height = 6)
# ggsave("antisense_intergenics_plot.png", projects_intergenics_plot, width = 10, height = 6, dpi = 300)