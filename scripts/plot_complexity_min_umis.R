# Plot detected-gene complexity across UMI thresholds for reviewer response.

library(tidyverse)
library(ggbeeswarm)
library(ggpubr)
library(patchwork)
library(here)

here::i_am("scripts/plot_complexity_min_umis.R")

thresholds <- c(1L, 5L, 10L)
threshold_levels <- c(">=1 UMI", ">=5 UMIs", ">=10 UMIs")

color_vector <- c(
  "42 \u00B0C" = "#DDC7E6",
  "50 \u00B0C" = "#3F1A4C",
  "original" = "#969994",
  "resuspended" = "#D28A4B"
)

validate_complexity <- function(data, required_cols, expected_ds_level, label) {
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(label, " is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  observed_thresholds <- sort(unique(data$min_umis))
  if (!identical(observed_thresholds, thresholds)) {
    stop(
      label, " has unexpected min_umis values: ",
      paste(observed_thresholds, collapse = ", ")
    )
  }

  observed_ds_levels <- sort(unique(data$ds_level))
  if (!identical(observed_ds_levels, expected_ds_level)) {
    stop(
      label, " has unexpected ds_level values: ",
      paste(observed_ds_levels, collapse = ", ")
    )
  }
}

add_threshold_label <- function(data) {
  data %>%
    mutate(
      threshold_label = factor(
        case_when(
          min_umis == 1L ~ ">=1 UMI",
          min_umis == 5L ~ ">=5 UMIs",
          min_umis == 10L ~ ">=10 UMIs"
        ),
        levels = threshold_levels
      )
    )
}

significance_label <- function(p_value) {
  case_when(
    p_value <= 0.0001 ~ "****",
    p_value <= 0.001 ~ "***",
    p_value <= 0.01 ~ "**",
    p_value <= 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

calculate_comparison_pvalues <- function(data, x_col, y_col, facet_col, comparisons) {
  comparison_tbl <- tibble(
    comparison_id = seq_along(comparisons),
    group1 = map_chr(comparisons, 1),
    group2 = map_chr(comparisons, 2)
  )

  data %>%
    group_by(.data[[facet_col]]) %>%
    group_modify(~ {
      map_dfr(seq_len(nrow(comparison_tbl)), function(comparison_idx) {
        group1 <- comparison_tbl$group1[comparison_idx]
        group2 <- comparison_tbl$group2[comparison_idx]
        group1_values <- .x[[y_col]][.x[[x_col]] == group1]
        group2_values <- .x[[y_col]][.x[[x_col]] == group2]

        tibble(
          comparison_id = comparison_idx,
          group1 = group1,
          group2 = group2,
          p = t.test(group1_values, group2_values, paired = FALSE)$p.value
        )
      })
    }) %>%
    ungroup() %>%
    mutate(p.signif = significance_label(p))
}

add_facet_y_positions <- function(
    pvalues,
    data,
    facet_col,
    y_col,
    step_multiplier = 0.16,
    start_multiplier = 0.08) {
  facet_ranges <- data %>%
    group_by(.data[[facet_col]]) %>%
    summarise(
      y_min = min(.data[[y_col]], na.rm = TRUE),
      y_max = max(.data[[y_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_range = y_max - y_min,
      y_range = if_else(y_range > 0, y_range, pmax(abs(y_max), 1)),
      y_start = y_max + y_range * start_multiplier,
      y_step = y_range * step_multiplier
    )

  pvalues %>%
    left_join(facet_ranges, by = facet_col) %>%
    mutate(y.position = y_start + (comparison_level - 1) * y_step)
}

rt_comparisons <- list(
  c("Maxima H-\n50 \u00B0C", "Maxima H-\n42 \u00B0C"),
  c("Superscript IV\n50 \u00B0C", "Superscript IV\n42 \u00B0C"),
  c("Maxima H-\n50 \u00B0C", "Superscript IV\n50 \u00B0C"),
  c("Maxima H-\n42 \u00B0C", "Superscript IV\n42 \u00B0C")
)

rt_comparison_levels <- tibble(
  comparison_id = seq_along(rt_comparisons),
  comparison_level = c(1, 1, 2, 3)
)

resusp_comparisons <- list(c("original", "resuspended"))

rt_complexity_raw <- readRDS(here("data/RTase/complexity_RTase_min_umis.rds"))
resusp_complexity_raw <- readRDS(here("data/BBB_Resusp/complexity_min_umis.rds"))

validate_complexity(
  rt_complexity_raw,
  required_cols = c("RG", "n", "project", "replicate", "ds_level", "min_umis", "ds_level2"),
  expected_ds_level = 7000000L,
  label = "RTase complexity"
)

validate_complexity(
  resusp_complexity_raw,
  required_cols = c("RG", "n", "project", "replicate", "ds_level", "min_umis"),
  expected_ds_level = 1250000L,
  label = "Resuspension complexity"
)

rt_complexity <- rt_complexity_raw %>%
  as_tibble() %>%
  filter(replicate == 1L, ds_level == 7000000L) %>%
  add_threshold_label() %>%
  mutate(
    temp = paste(str_extract(project, "..$"), "\u00B0C"),
    temp = factor(temp, levels = c("42 \u00B0C", "50 \u00B0C")),
    enzyme = if_else(grepl("Maxima_..", project), "Maxima H-", "Superscript IV")
  ) %>%
  mutate(
    temp_enz_cat = factor(
      paste0(enzyme, "\n", temp),
      levels = c(
        "Maxima H-\n42 \u00B0C",
        "Maxima H-\n50 \u00B0C",
        "Superscript IV\n42 \u00B0C",
        "Superscript IV\n50 \u00B0C"
      )
    )
  )

rt_complexity_pvalues <- calculate_comparison_pvalues(
  rt_complexity,
  x_col = "temp_enz_cat",
  y_col = "n",
  facet_col = "threshold_label",
  comparisons = rt_comparisons
) %>%
  left_join(rt_comparison_levels, by = "comparison_id") %>%
  add_facet_y_positions(
    data = rt_complexity,
    facet_col = "threshold_label",
    y_col = "n",
    step_multiplier = 0.18,
    start_multiplier = 0.08
  )

rt_complexity_plot <- ggplot(
  rt_complexity,
  aes(y = n, x = temp_enz_cat, color = temp)
) +
  geom_quasirandom(alpha = 0.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white") +
  stat_pvalue_manual(
    rt_complexity_pvalues,
    label = "p.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.3,
    size = 3,
    inherit.aes = FALSE
  ) +
  facet_grid(threshold_label ~ ., scales = "free_y") +
  labs(y = "Detected Genes", x = "") +
  theme_pubr(legend = "none") +
  scale_color_manual(values = color_vector) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.12))) +
  coord_cartesian(clip = "off") +
  theme(axis.text.x = element_text(size = 13, angle = 90, vjust = 0.5, hjust = 1),
    strip.text.x = element_text(size = 13),
    strip.text.y = element_text(size = 13)
  )

resusp_complexity <- resusp_complexity_raw %>%
  as_tibble() %>%
  filter(project != "BBB", ds_level == 1250000L) %>%
  group_by(project, RG, min_umis) %>%
  summarise(n = mean(n), .groups = "drop") %>%
  add_threshold_label() %>%
  mutate(
    project = factor(
      if_else(project == "Std", "original", "resuspended"),
      levels = c("original", "resuspended")
    )
  )

resusp_complexity_pvalues <- calculate_comparison_pvalues(
  resusp_complexity,
  x_col = "project",
  y_col = "n",
  facet_col = "threshold_label",
  comparisons = resusp_comparisons
) %>%
  mutate(comparison_level = 1) %>%
  add_facet_y_positions(
    data = resusp_complexity,
    facet_col = "threshold_label",
    y_col = "n",
    step_multiplier = 0.12,
    start_multiplier = 0.08
  )

resusp_complexity_plot <- ggplot(
  resusp_complexity,
  aes(y = n, x = project, color = project)
) +
  geom_violin() +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0, linewidth = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, stroke = 0.6, fill = "white") +
  stat_pvalue_manual(
    resusp_complexity_pvalues,
    label = "p.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.3,
    size = 3,
    inherit.aes = FALSE
  ) +
  facet_wrap(~threshold_label, nrow = 1, scales = "free_y") +
  labs(y = "Detected Genes", x = "") +
  theme_pubr(legend = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.12))) +
  coord_cartesian(clip = "off") +
  theme(axis.text.x = element_text(size = 13, angle = 90, vjust = 0.5, hjust = 1),
    strip.text.x = element_text(size = 13),
    strip.text.y = element_text(size = 13)
  )+
  scale_color_manual(values = color_vector)

complexity_min_umis_plot <- (rt_complexity_plot | resusp_complexity_plot) +
  plot_layout(widths = c(1, 1.5)) +
  plot_annotation(tag_levels = "A")

output_file <- here("figures/fig_S_complexity_min_umis.pdf")

ggsave(
  plot = complexity_min_umis_plot,
  filename = output_file,
  height = 5.5,
  width = 9,
  device = cairo_pdf
)

if (!file.exists(output_file) || file.info(output_file)$size == 0) {
  stop("Failed to create non-empty output file: ", output_file)
}
