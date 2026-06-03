# load packages
library(data.table)
library(parallel)
library(Rsamtools)
library(here)

here::i_am("scripts/Complexity_Resusp_min_umis.R")

reads2genes_new_ds <- function(featfile, bccount, inex, chunk, cores, downsampling) {
  chunk_bcs <- bccount[chunkID == chunk]$XC
  idxstats <- Rsamtools::idxstatsBam(featfile)
  if ("*" %in% idxstats$seqnames) {
    idxstats <- idxstats[idxstats$seqnames != "*", ]
    idxstats$seqnames <- as.character(idxstats$seqnames)
  }

  taglist <- c("BC", "UB", "GE")
  if (inex) {
    taglist <- c(taglist, "GI")
  }

  rsamtools_reads <- parallel::mclapply(seq_len(nrow(idxstats)), function(x) {
    parms <- Rsamtools::ScanBamParam(
      tag = taglist,
      what = "pos",
      tagFilter = list(BC = chunk_bcs),
      which = GenomicRanges::GRanges(
        seqnames = idxstats[x, "seqnames"],
        ranges = IRanges::IRanges(1, idxstats[x, "seqlength"])
      )
    )

    dat <- Rsamtools::scanBam(file = featfile, param = parms)
    if (inex) {
      data.table(
        RG = dat[[1]]$tag$BC,
        UB = dat[[1]]$tag$UB,
        GE = dat[[1]]$tag$GE,
        GEin = dat[[1]]$tag$GI
      )
    } else {
      data.table(
        RG = dat[[1]]$tag$BC,
        UB = dat[[1]]$tag$UB,
        GE = dat[[1]]$tag$GE
      )
    }
  }, mc.cores = cores)

  rsamtools_reads <- rbindlist(rsamtools_reads, fill = TRUE, use.names = TRUE)

  if (downsampling != FALSE && nrow(rsamtools_reads) > downsampling) {
    rsamtools_reads <- rsamtools_reads[sample(.N, downsampling)]
  }

  if (inex) {
    rsamtools_reads[, ftype := "NA"]
    rsamtools_reads[!is.na(GEin), ftype := "intron"]
    rsamtools_reads[!is.na(GE), ftype := "exon"]
    rsamtools_reads[is.na(GE), GE := GEin]
    rsamtools_reads[, GEin := NULL]
  } else {
    rsamtools_reads[, ftype := "NA"]
    rsamtools_reads[!is.na(GE), ftype := "exon"]
  }
  setkey(rsamtools_reads, RG)

  rsamtools_reads
}

splitRG <- function(bccount, mem, hamdist) {
  if (is.null(mem) || mem == 0) {
    maxR <- Inf
  } else {
    maxR <- floor(mem * 1000 * 4500)
  }
  if (maxR > 2e+09) {
    maxR <- 2e+09
  }
  if (hamdist > 0) {
    maxR <- floor(maxR / 3)
  }

  print(paste(maxR, "Reads per chunk"))
  nc <- nrow(bccount)
  cs <- 0
  chunkID <- 1
  bccount[, chunkID := 1]
  if (sum(bccount$n) > maxR) {
    for (i in seq_len(nc)) {
      cs <- cs + bccount[i]$n
      if (bccount[i]$n > maxR) {
        print(paste(
          "Warning: Barcode", bccount[i]$XC,
          "has more reads than allowed for the memory limit. Proceeding anyway..."
        ))
      }
      if (cs >= maxR) {
        if (i > 1) {
          chunkID <- chunkID + 1
        }
        cs <- bccount[i, n]
      }
      bccount[i, chunkID := chunkID]
    }
  }
  bccount
}

count_thresholded_complexity <- function(reads, barcodes, thresholds = c(1L, 5L, 10L)) {
  gene_umis <- as.data.table(reads)[!is.na(GE) & !is.na(UB), .(RG, GE, UB)]
  gene_umis <- unique(gene_umis)[, .(umi_n = uniqueN(UB)), by = .(RG, GE)]

  rbindlist(lapply(thresholds, function(min_umis) {
    counts <- gene_umis[umi_n >= min_umis, .(n = .N), by = RG]
    counts <- merge(data.table(RG = barcodes), counts, by = "RG", all.x = TRUE)
    counts[is.na(n), n := 0L]
    counts[, min_umis := min_umis]
    counts
  }), use.names = TRUE)
}

validate_complexity <- function(complexity, expected_ds_level, thresholds = c(1L, 5L, 10L)) {
  stopifnot(setequal(sort(unique(complexity$min_umis)), thresholds))
  stopifnot(identical(sort(unique(complexity$ds_level)), expected_ds_level))

  meta_cols <- c("RG", "project", "replicate", "ds_level", "min_umis")
  stopifnot(!any(vapply(complexity[, ..meta_cols], anyNA, logical(1))))

  wide <- dcast(
    complexity,
    RG + project + replicate + ds_level ~ min_umis,
    value.var = "n"
  )
  stopifnot(all(wide[["1"]] >= wide[["5"]]))
  stopifnot(all(wide[["5"]] >= wide[["10"]]))
}

samples <- c("Std", "BBB", "Resusp")
ds_level <- 1250000L
replicates <- seq_len(10L)

setDTthreads(threads = 30)

results <- list()
result_idx <- 1L

for (sample_idx in seq_along(samples)) {
  sample_name <- samples[[sample_idx]]

  for (replicate_idx in replicates) {
    print(paste0("sample ", sample_name, ", rep ", replicate_idx))

    inputBAM <- paste0(
      "/data/share/htp/Felix_genotyping/primeSeq_tests/BBB_Resusp/trimmed/",
      sample_name, "/", sample_name, ".filtered.Aligned.GeneTagged.sorted.bam"
    )
    bccount_file <- paste0(
      "/data/share/htp/Felix_genotyping/primeSeq_tests/BBB_Resusp/trimmed/",
      sample_name, "/zUMIs_output/", sample_name, "kept_barcodes_binned.txt"
    )

    bccount <- fread(bccount_file)
    bccount <- splitRG(bccount = bccount, mem = 50, hamdist = 0)

    set.seed(20260522 + sample_idx * 100 + replicate_idx)
    reads <- reads2genes_new_ds(
      featfile = inputBAM,
      bccount = bccount,
      inex = TRUE,
      chunk = 1,
      cores = 20,
      downsampling = ds_level
    )
    print(nrow(reads))

    counts <- count_thresholded_complexity(reads, barcodes = bccount$XC)
    counts[, `:=`(
      project = sample_name,
      replicate = replicate_idx,
      ds_level = ds_level
    )]
    setcolorder(counts, c("RG", "n", "project", "replicate", "ds_level", "min_umis"))

    results[[result_idx]] <- counts
    result_idx <- result_idx + 1L
    rm(reads, counts)
  }
}

complexity_min_umis <- rbindlist(results, use.names = TRUE)
validate_complexity(complexity_min_umis, expected_ds_level = ds_level)

saveRDS(
  complexity_min_umis,
  here("data/BBB_Resusp/complexity_min_umis.rds")
)
