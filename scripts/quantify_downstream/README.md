# Quantify Intergenic Reads Downstream of Genes

## What it does
This script quantifies intergenic reads that map within 10kb downstream of genes for multiple BAM files. It identifies reads that are not assigned to any gene (intergenic) and determines how many of these fall within a specified distance downstream of gene 3' ends.

## How it works

### Major Steps:
1. **Gene Range Processing**: Extracts gene coordinates from GTF annotation file and converts to BED format
2. **Intergenic Read Extraction**: 
   - Reads BAM files and filters for reads with ES="Unassigned_NoFeatures" and IS="Unassigned_NoFeatures"
   - Optionally filters by barcode whitelist
   - Converts intergenic reads to BED format
3. **Downstream Analysis**: Uses bedtools to find intergenic reads within 10kb downstream of gene 3' ends
4. **Statistical Analysis**: Calculates distance distributions, mean/median distances, and gene-level summaries
5. **Visualization**: Creates individual and combined histograms with mean/median lines and comparison plots

## Inputs Required

### Files:
- **GTF annotation file**: Gene annotation in GTF format (e.g., GENCODE)
- **BAM files**: Aligned reads with ES, IS, and BC tags (from zUMIs pipeline)
- **Barcode whitelist files**: Text files with valid barcodes (one per line)
- **Sample configuration file**: `bam_files.txt` with columns:
  - `project`: Sample name
  - `bam_path`: Path to BAM file
  - `BC_WL_path`: Path to barcode whitelist file

### Example `bam_files.txt`:
```
project	bam_path	BC_WL_path
PoP64_PTO	/data/path/sample.bam	/data/path/barcodes.txt
```

### Dependencies:
- R with Bioconductor packages: `GenomicAlignments`, `GenomicRanges`, `rtracklayer`, `stringr`, `dplyr`, `ggplot2`
- Command line tools: `bedtools`, `gtf2bed`, `grep`, `awk`, `sortBed`

## Outputs Created

### Files in `results/` directory:
- **`gene_ranges.bed`**: Gene coordinates in BED format
- **`{sample}_intergenic_reads.bed`**: Intergenic reads for each sample
- **`{sample}_distance_histogram.pdf`**: Individual histogram with mean/median lines
- **`{sample}_gene_summary.txt`**: Gene-level read counts and distances
- **`combined_distance_histograms.pdf`**: All samples in faceted plot
- **`comparison_plot.pdf`**: Bar plot comparing downstream fractions across samples
- **`summary_report.txt`**: Summary statistics for all samples

### Summary Statistics:
- Total intergenic reads per sample
- Number of reads within 10kb downstream
- Fraction of intergenic reads that are downstream
- Mean and median distances from gene 3' ends
- Number of unique genes with downstream reads

## Usage
```bash
Rscript quantify_intergenic_reads.R
```

## Configuration
Modify these variables in the script:
- `GTF_FILE`: Path to GTF annotation file
- `DOWNSTREAM_DISTANCE`: Distance threshold (default: 10000 bp)
- `OUTPUT_DIR`: Output directory path 