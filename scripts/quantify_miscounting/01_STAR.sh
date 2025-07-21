#!/usr/bin/env bash
set -euo pipefail

#### 1) Settings from your zUMIs.yaml ####
PROJECT="${PROJECT:-PoP64_PTO}" # <-------------------------------------------------- adapt
THREADS=10
MEM_GB=80

# Single‐end cDNA file (prime‐seq SE layout)
FASTQ_CDNA="${FASTQ_CDNA:-/data/share/htp/prime-seq_NextGen/data/FC2024_11_01_PoP64_1/02_trimming/lane1_prime-seq_PoP64_PTO_corrected_r2_trimmed.fq.gz}" # <--adapt

# Reference
STAR_INDEX="/data/share/htp/Felix_genotyping/zUMIs_tests/own_genomes/mus_musculus/STAR_idx_mouse_GRCm39/"
GTF_FILE="/data/share/htp/Felix_genotyping/zUMIs_tests/own_genomes/mus_musculus/gencode.vM34.primary_assembly.annotation.gtf"

# Executables
STAR_EXEC="/data/share/htp/prime-seq_Data_Processing/00_tools/zUMIs_2.9.7e/zUMIs-env/bin/STAR"
SAMTOOLS="samtools"
PIGZ="pigz"

# Output dirs
OUTDIR="/home/felix/prime-seq_NextGen/scripts/quantify_miscounting/STAR_out"
mkdir -p "$OUTDIR"

#### 2) “ALL‐BEST” STAR run (capture every equal‐best hit) ####
echo "### Running ALL‐BEST STAR (up to 50 best alignments per read) ###"
"$STAR_EXEC" \
  --runThreadN "$THREADS" \
  --genomeDir "$STAR_INDEX" \
  --sjdbGTFfile "$GTF_FILE" \
  --readFilesIn "$FASTQ_CDNA" \
  --readFilesCommand "$PIGZ" -dc \
  --twopassMode Basic \
  --outFilterMultimapNmax 50 \
  --outSAMmultNmax 50 \
  --outSAMprimaryFlag AllBestScore \
  --limitBAMsortRAM $(( MEM_GB * 1000000000 )) \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix "${OUTDIR}/${PROJECT}_allbest_"
