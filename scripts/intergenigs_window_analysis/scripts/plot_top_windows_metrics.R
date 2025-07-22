#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(gridExtra)
})

# Number of top locations to plot (set to 0 to plot all)
num_top_locations <- 20  # Change this value to plot fewer or more locations

# Define colors for conditions
pool_colors <- c("80ng" = "#C4C082",
                 "320ng" = "#7D7A3B",
                 "920ng" = "#454321")

# Read top 100 windows
top_windows <- fread("/data/share/htp/prime-seq_NextGen/intergenics_analysis/plot/top_100_windows_mean_across_conditions.txt")
head(top_windows)

# Filter to only include the top X locations if specified
if (num_top_locations > 0) {
    top_windows <- top_windows[1:min(num_top_locations, nrow(top_windows)), ]
    cat("Plotting top", nrow(top_windows), "locations\n")
} else {
    cat("Plotting all", nrow(top_windows), "locations\n")
}

# Read metrics files
metrics_920ng <- fread("/data/share/htp/prime-seq_NextGen/intergenics_analysis/output/poolsize_920ng.filtered.Aligned.GeneTagged.sorted_window_metrics.txt")
metrics_320ng <- fread("/data/share/htp/prime-seq_NextGen/intergenics_analysis/output/poolsize_320ng.filtered.Aligned.GeneTagged.sorted_window_metrics.txt")
metrics_80ng <- fread("/data/share/htp/prime-seq_NextGen/intergenics_analysis/output/poolsize_80ng.filtered.Aligned.GeneTagged.sorted_window_metrics.txt")
head(metrics_920ng)

# Function to standardize chromosome names
standardize_chr <- function(chr) {
    # Remove 'chr' prefix if present
    chr <- gsub("^chr", "", chr)
    return(chr)
}

# Function to match windows and get metrics
get_window_metrics <- function(metrics_dt, windows_dt) {
    # Standardize chromosome names in both data tables
    metrics_dt[, chr := standardize_chr(chr)]
    windows_dt[, chr := standardize_chr(chr)]
    
    # Create window keys
    metrics_dt[, window_key := paste(chr, start, end, sep=":")]
    windows_dt[, window_key := paste(chr, start, end, sep=":")]
    
    # Merge the data
    matched_metrics <- merge(
        windows_dt[, .(window_key, location)],
        metrics_dt,
        by="window_key"
    )
    
    cat("Number of matched windows:", nrow(matched_metrics), "\n")
    
    matched_metrics[, window_key := NULL]
    return(matched_metrics)
}

# Get metrics for each condition
metrics_920ng <- get_window_metrics(metrics_920ng, top_windows)
head(metrics_920ng)
metrics_320ng <- get_window_metrics(metrics_320ng, top_windows)
metrics_80ng <- get_window_metrics(metrics_80ng, top_windows)

# Add condition labels with proper ordering
metrics_80ng[, condition := "80ng"]
metrics_320ng[, condition := "320ng"]
metrics_920ng[, condition := "920ng"]

# Combine all metrics
all_metrics <- rbindlist(list(metrics_80ng, metrics_320ng, metrics_920ng))

# Order locations by the ranking in top_windows
location_order <- top_windows$location
all_metrics[, location := factor(location, levels = location_order)]

# Set condition as a factor with specific order
all_metrics[, condition := factor(condition, levels = c("80ng", "320ng", "920ng"))]

# Create MAPQ plot
mapq_plot <- ggplot(all_metrics, aes(x = location, y = mean_mapq, color = condition)) +
    geom_point() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(x = "Window", y = "Mean MAPQ", title = "Mean MAPQ by Window and Condition") +
    scale_color_manual(values = pool_colors)

# Create strand ratio plot with explanation
strand_plot <- ggplot(all_metrics, aes(x = location, y = strand_ratio, color = condition)) +
    geom_point() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(
        x = "Window", 
        y = "Forward fraction", 
        title = "Strand Ratio by Window and Condition"
    ) +
    scale_color_manual(values = pool_colors) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
    annotate("text", x = Inf, y = 0.5, hjust = 1.1, vjust = 1.1, 
             label = "Equal forward/reverse", color = "gray50", size = 3)

# Create output filename based on number of locations
output_filename <- paste0("/data/share/htp/prime-seq_NextGen/intergenics_analysis/output/top_windows_metrics_plots")
if (num_top_locations > 0) {
    output_filename <- paste0(output_filename, "_top", num_top_locations)
}
output_filename <- paste0(output_filename, ".pdf")

# Save plots
pdf(output_filename, width = 12, height = 10)  # Increased height to accommodate subtitle
grid.arrange(mapq_plot, strand_plot, ncol = 1)
dev.off()