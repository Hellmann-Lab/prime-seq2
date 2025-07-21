#!/bin/bash
#SBATCH --job-name=STAR
#SBATCH --output=/home/felix/prime-seq_NextGen/scripts/quantify_miscounting/STAR_out/STAR_%j.out
#SBATCH --error=/home/felix/prime-seq_NextGen/scripts/quantify_miscounting/STAR_out/STAR_%j.err
#SBATCH --cpus-per-task=10
#SBATCH --mem=80G

PROJECT="${PROJECT:-PoP64_PTO}"

bash 01_STAR.sh && bash 02_quantify_multimappers.sh && /opt/R/4.4.1/lib/R/bin/Rscript 03_summarize_multimappers.R

# Remove intermediate files
OUTDIR="/home/felix/prime-seq_NextGen/scripts/quantify_miscounting/STAR_out"
PREFIX="$OUTDIR/${PROJECT}_allbest"
rm -f ${PREFIX}_reads_all.txt ${PREFIX}_reads_exon.txt ${PREFIX}_reads_intron.txt ${PREFIX}_reads_any_intergenic.txt ${PREFIX}_reads_intergenic.txt ${PREFIX}_SJ.out.tab  ${PREFIX}_Log* ${PREFIX}_Aligned.sortedByCoord.out.bam ${PREFIX}_Aligned.sortedByCoord.out.bam.bai ${PREFIX}.bam ${PREFIX}.bam.bai ${PREFIX}.bed
rm -rf ${PREFIX}__STAR*
rm -rf "$OUTDIR/annotation_${PROJECT}"
