# Aim: Quantify what fraction of the genome is covered by Exon / Intron annotations
## ================================================================
## 1.  Load packages
## ------------------------------------------------
library(GenomicFeatures)
library(GenomeInfoDb)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm39)

## ================================================================
## 2.  Bring in the GENCODE mouse annotation
## ------------------------------------------------
gtf <- "/data/share/htp/Felix_genotyping/zUMIs_tests/own_genomes/mus_musculus/gencode.vM34.primary_assembly.annotation.gtf"
txdb <- makeTxDbFromGFF(gtf, format = "gtf")

## Harmonise to UCSC “chr*” style used by the BSgenome object
seqlevelsStyle(txdb) <- "UCSC"

## Optional: keep only main chromosomes (chr1-19,X,Y,M)
# keep <- standardChromosomes(txdb)      # same call later for the genome
# txdb <- keepSeqlevels(txdb, keep, pruning.mode = "coarse")

## ================================================================
## 3.  Build disjoint exonic & intronic parts (strand-aware)
## ------------------------------------------------
exon.parts   <- exonicParts  (txdb, linked.to.single.gene.only = FALSE)
intron.parts <- intronicParts(txdb, linked.to.single.gene.only = FALSE)

exon.u   <- reduce(exon.parts,   ignore.strand = FALSE)  # unique exons
intron.u <- reduce(intron.parts, ignore.strand = FALSE)  # unique introns
annot.u  <- reduce(c(exon.u, intron.u), ignore.strand = FALSE)

## ================================================================
## 4.  Intergenic space (genome − all annotated bases, but keep strand)
##    Intergenic intervals come back with strand “*”.
## ------------------------------------------------
# genome.gr <- as(seqinfo(BSgenome.Mmusculus.UCSC.mm39)[keep], "GRanges")
genome.gr <- as(seqinfo(BSgenome.Mmusculus.UCSC.mm39), "GRanges")
intergenic.raw <- setdiff(genome.gr, annot.u, ignore.strand = FALSE)

## Copy intergenic length to both strands so every class has “+” and “–”
intergenic <- unstrand(intergenic.raw)          # drop the “*”
strand(intergenic) <- "+"                       # first copy
intergenic.plus   <- intergenic
strand(intergenic) <- "-"                       # second copy
intergenic.minus  <- intergenic
intergenic <- c(intergenic.plus, intergenic.minus)

## ================================================================
## 5.  Helper to sum base-pairs by strand (always returns + and –)
## ------------------------------------------------
sum_by_strand <- function(gr, strands = c("+", "-")) {
  out <- tapply(width(gr),
                factor(strand(gr), levels = strands),
                sum, default = 0)
  as.numeric(out)                     # return unnamed numeric vector
}

## ================================================================
## 6.  Build the matrix of absolute bp
## ------------------------------------------------
strands <- c("+", "-")

bp <- rbind(
  exon       = sum_by_strand(exon.u,   strands),
  intron     = sum_by_strand(intron.u, strands),
  intergenic = sum_by_strand(intergenic, strands)
)

colnames(bp) <- strands

## Add a denominator row (genome length per strand)
bp <- rbind(bp, genome = colSums(bp))

## Fractions
fractions <- sweep(bp[1:3, ], 2, bp["genome", ], "/")

## ================================================================
## 7.  Show results
## ------------------------------------------------
round(fractions, 4)

# intergenic
cat("Fraction that is not annotated as Exons or Introns:", round(mean(fractions["intergenic",]),4)*100, "%")
cat("Fraction that is annotated as Exons or Introns:", round(mean(rowSums(fractions[1:2,])),4)*100, "%")
