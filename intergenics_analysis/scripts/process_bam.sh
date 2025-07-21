#!/bin/bash

# Check if input BAM file is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <input.bam>"
    exit 1
fi

# Get absolute path of input BAM
INPUT_BAM=$(readlink -f "$1")

# Create output directory in the same location as the script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUTPUT_DIR="$SCRIPT_DIR/../output"
mkdir -p "$OUTPUT_DIR"

# Set output prefix
OUTPUT_PREFIX="$OUTPUT_DIR/$(basename ${INPUT_BAM%.bam})"

# Extract intergenic reads and create BED file in one step
echo "Processing BAM file..."
samtools view -h "$INPUT_BAM" | \
    grep -E '^@|ES:Z:Unassigned_NoFeatures.*IS:Z:Unassigned_NoFeatures' | \
    samtools view - | \
    awk '{umi="NA"; for(i=1;i<=NF;i++) if($i ~ /^UB:Z:/) {split($i,a,":"); umi=a[3]}; print $3 "\t" $4 "\t" $4+length($10) "\t" $1 "\t" umi}' > "${OUTPUT_PREFIX}_intergenic.bed"

# Run window counting
echo "Counting reads and UMIs per window..."
/opt/R/4.4.1/lib/R/bin/Rscript "$SCRIPT_DIR/count_windows.R" "${OUTPUT_PREFIX}_intergenic.bed"

# Run validation
echo "Validating counts..."
/opt/R/4.4.1/lib/R/bin/Rscript "$SCRIPT_DIR/validate_counts.R" "${OUTPUT_PREFIX}_intergenic_window_counts.txt" "${OUTPUT_PREFIX}_intergenic.bed"

# Run window metrics analysis
echo "Analyzing MAPQ and strand distribution per window..."
/opt/R/4.4.1/lib/R/bin/Rscript "$SCRIPT_DIR/analyze_window_metrics.R" "$INPUT_BAM"

if [ $? -eq 0 ]; then
    echo "Processing complete. Output files:"
    echo "- ${OUTPUT_PREFIX}_intergenic_window_counts.txt"
    echo "- ${OUTPUT_PREFIX}_window_metrics.txt"
else
    echo "Processing failed. Please check the output above."
    exit 1
fi 