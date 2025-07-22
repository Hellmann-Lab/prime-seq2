#!/bin/bash

# List of BAM files to process
BAMS=(
    "/data/share/htp/prime-seq_NextGen/data/FC2024_08_01_poolsize/03_zUMIs/80ng/poolsize_80ng.filtered.Aligned.GeneTagged.sorted.bam"
    "/data/share/htp/prime-seq_NextGen/data/FC2024_08_01_poolsize/03_zUMIs/320ng/poolsize_320ng.filtered.Aligned.GeneTagged.sorted.bam"
    "/data/share/htp/prime-seq_NextGen/data/FC2024_08_01_poolsize/03_zUMIs/920ng/poolsize_920ng.filtered.Aligned.GeneTagged.sorted.bam"
)

# Process each BAM file
for bam in "${BAMS[@]}"; do
    echo "Processing $bam..."
    ./scripts/process_bam.sh "$bam"
done 