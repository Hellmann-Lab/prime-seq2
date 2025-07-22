#!/bin/bash

#########################################################
#### Run Intergenic Read Downstream Analysis ####
#########################################################

# This script runs the intergenic read analysis for all BAM files
# It checks dependencies and provides proper error handling

set -e  # Exit on any error

echo "=== Intergenic Read Downstream Analysis ==="
echo "Started at: $(date)"
echo

#### Check dependencies ####
###########################

echo "Checking dependencies..."

# Check if R is available
if ! command -v R &> /dev/null; then
    echo "ERROR: R is not installed or not in PATH"
    exit 1
fi

# Check if bedtools is available
if ! command -v bedtools &> /dev/null; then
    echo "ERROR: bedtools is not installed or not in PATH"
    echo "Please install bedtools: https://bedtools.readthedocs.io/en/latest/content/installation.html"
    exit 1
fi

# Check if gtf2bed is available (part of bedtools)
if ! command -v gtf2bed &> /dev/null; then
    echo "ERROR: gtf2bed is not available. This is part of bedtools."
    echo "Please ensure bedtools is properly installed."
    exit 1
fi

# Check if sortBed is available (part of bedtools)
if ! command -v sortBed &> /dev/null; then
    echo "ERROR: sortBed is not available. This is part of bedtools."
    echo "Please ensure bedtools is properly installed."
    exit 1
fi

echo "✓ All dependencies found"

#### Check input files ####
##########################

echo "Checking input files..."

# Check if BAM files list exists
if [ ! -f "/data/share/htp/prime-seq_NextGen/scripts/quantify_downstream/bam_files.txt" ]; then
    echo "ERROR: bam_files.txt not found at /data/share/htp/prime-seq_NextGen/scripts/quantify_downstream/bam_files.txt"
    exit 1
fi

# Check if GTF file exists
GTF_FILE="/data/share/htp/Felix_genotyping/zUMIs_tests/own_genomes/mus_musculus/gencode.vM34.primary_assembly.annotation.gtf"
if [ ! -f "$GTF_FILE" ]; then
    echo "ERROR: GTF file not found at $GTF_FILE"
    exit 1
fi

# Check if all BAM files and barcode whitelists exist
echo "Checking BAM files and barcode whitelists..."
while IFS=$'\t' read -r project bam_path bc_wl_path; do
    if [ "$project" != "project" ]; then  # Skip header
        if [ ! -f "$bam_path" ]; then
            echo "ERROR: BAM file not found: $bam_path"
            exit 1
        fi
        
        if [ ! -f "$bc_wl_path" ]; then
            echo "ERROR: Barcode whitelist not found: $bc_wl_path"
            exit 1
        fi
        
        # Check if BAM index exists
        if [ ! -f "${bam_path}.bai" ] && [ ! -f "${bam_path%.bam}.bai" ]; then
            echo "WARNING: BAM index not found for $bam_path"
            echo "Attempting to create index..."
            samtools index "$bam_path"
        fi
        
        echo "✓ $project: BAM and whitelist found"
    fi
done < /data/share/htp/prime-seq_NextGen/scripts/quantify_downstream/bam_files.txt

echo "✓ All input files found"

#### Create output directory ####
################################

OUTPUT_DIR="/data/share/htp/prime-seq_NextGen/scripts/quantify_downstream/results"
mkdir -p "$OUTPUT_DIR"

echo "Output directory: $OUTPUT_DIR"

#### Run the analysis ####
#########################

echo "Starting R analysis..."
echo

# Run the R script
/opt/R/4.4.1/lib/R/bin/Rscript /data/share/htp/prime-seq_NextGen/scripts/quantify_downstream/quantify_intergenic_reads.R

#### Check results ####
######################

echo
echo "=== Checking Results ==="

if [ -f "$OUTPUT_DIR/summary_report.txt" ]; then
    echo "✓ Summary report created"
    echo
    echo "Summary of results:"
    cat "$OUTPUT_DIR/summary_report.txt"
else
    echo "WARNING: Summary report not found"
fi

# List all output files
echo
echo "Generated files:"
ls -la "$OUTPUT_DIR/"

echo
echo "=== Analysis Complete ==="
echo "Finished at: $(date)"
echo "Results saved in: $OUTPUT_DIR" 