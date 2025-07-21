#!/usr/bin/env bash
set -euo pipefail

#### Config (adjust if paths differ) ####
OUTDIR="/home/felix/prime-seq_NextGen/scripts/quantify_miscounting/STAR_out"
PROJECT="${PROJECT:-PoP64_PTO}"
GTF="/data/share/htp/Felix_genotyping/zUMIs_tests/own_genomes/mus_musculus/gencode.vM34.primary_assembly.annotation.gtf"
BAM_ALLBEST="${OUTDIR}/${PROJECT}_allbest_Aligned.sortedByCoord.out.bam"
PREFIX="${OUTDIR}/${PROJECT}_allbest"

# tools
SAMTOOLS="samtools"
BEDTOOLS="bedtools"
AWK="awk"

#### 1) Index (if not already) ####
$SAMTOOLS index "$BAM_ALLBEST"

#### 2) Build exon & intron BEDs (once per annotation) ####
anno_dir="${OUTDIR}/annotation_${PROJECT}"
mkdir -p "$anno_dir"

# Exons
$AWK '$3=="exon" {
  split($10,f,"\"");
  print $1, $4-1, $5, f[2]
}' OFS="\t" "$GTF" \
  > "${anno_dir}/exons.bed"

# Introns = gene body minus exons
$AWK '$3=="gene" {
  split($10,f,"\"");
  print $1, $4-1, $5, f[2]
}' OFS="\t" "$GTF" \
  > "${anno_dir}/genes.bed"

$BEDTOOLS subtract \
  -a "${anno_dir}/genes.bed" \
  -b "${anno_dir}/exons.bed" \
  > "${anno_dir}/introns.bed"


#### 3) Extract multi‐mapping read IDs ####
# -F 0x100 -> not secondary alignment
# multiple primary alignments (uniq -d)
$SAMTOOLS view -F 0x100 "$BAM_ALLBEST" | awk '{print $1}' | sort | uniq -d > "${PREFIX}_reads_all.txt"


#### 4) Subset BAM to multi‐mappers only ####
$SAMTOOLS view -h "$BAM_ALLBEST" | \
  awk 'NR==FNR{r[$1]; next} /^@/ {print; next} ($1 in r)' "${PREFIX}_reads_all.txt" - | \
  $SAMTOOLS view -bS - > "${PREFIX}.bam"
$SAMTOOLS index "${PREFIX}.bam"


#### 5) Convert to BED (one record per alignment) ####
$BEDTOOLS bamtobed \
  -i "${PREFIX}.bam" \
  > "${PREFIX}.bed"

# 5b) reads with *any* intergenic alignment (at least one alignment that misses both exons&introns)
# merge exon+intron BEDs
cat "${anno_dir}/exons.bed" "${anno_dir}/introns.bed" > "${anno_dir}/genic.bed"

# find alignments not hitting genic regions
bedtools intersect -v -a "${PREFIX}.bed" -b "${anno_dir}/genic.bed" \
  | cut -f4 \
  | sort -u \
  > "${PREFIX}_reads_any_intergenic.txt"


#### 6) Find which reads hit exons, introns, or neither ####
# exon
$BEDTOOLS intersect -wa -u \
  -a "${PREFIX}.bed" \
  -b "${anno_dir}/exons.bed" \
  | cut -f4 \
  | sort -u \
  > "${PREFIX}_reads_exon.txt"

# intron
$BEDTOOLS intersect -wa -u \
  -a "${PREFIX}.bed" \
  -b "${anno_dir}/introns.bed" \
| $BEDTOOLS intersect -v \
  -a - \
  -b "${anno_dir}/exons.bed" \
| cut -f4 \
| sort -u \
> "${PREFIX}_reads_intron.txt"

# # intergenic = all minus (exon ∪ intron)
# comm -23 \
#   <(sort "${PREFIX}_reads_all.txt") \
#   <(sort -m "${PREFIX}_reads_exon.txt" "${PREFIX}_reads_intron.txt") \
#   > "${PREFIX}_reads_intergenic.txt"

echo "Done extracting multi‐mapper categories in ${PREFIX}_*"
