# load packages
library(tidyverse)
library(ShortRead)
library(data.table)
library(parallel)
library(Rsamtools)
library(here)

here::i_am("scripts/Complexity_RTase.R")

# define functions
reads2genes_new_ds <- function(featfile, bccount, inex, chunk, cores, downsampling, keepUnassigned = FALSE){
  chunk_bcs <- bccount[chunkID==chunk]$XC
  idxstats <- Rsamtools::idxstatsBam(featfile)
  if("*" %in% idxstats$seqnames){
    idxstats <- idxstats[idxstats$seqnames != "*", ]
    idxstats$seqnames <- as.character(idxstats$seqnames)
  }
  taglist <- c("BC", "UB","GE")
  if(inex){
    taglist <- c(taglist, "GI")
  }
  
  rsamtools_reads <- parallel::mclapply(1:nrow(idxstats), function(x) {
    if(FALSE){
      parms <- ScanBamParam(tag=taglist,
                            what="pos",
                            flag = scanBamFlag(isFirstMateRead = TRUE),
                            tagFilter = list(BC = chunk_bcs),
                            which = GRanges(seqnames = idxstats[x,"seqnames"], ranges = IRanges(1,idxstats[x,"seqlength"])))
    }else{
      parms <- ScanBamParam(tag=taglist,
                            what="pos",
                            tagFilter = list(BC = chunk_bcs),
                            which = GRanges(seqnames = idxstats[x,"seqnames"], ranges = IRanges(1,idxstats[x,"seqlength"])))
    }
    
    dat <- scanBam(file = featfile, param = parms)
    if(inex){
      dt <- data.table(RG = dat[[1]]$tag$BC, UB = dat[[1]]$tag$UB, GE = dat[[1]]$tag$GE, GEin = dat[[1]]$tag$GI)
    }else{
      dt <- data.table(RG = dat[[1]]$tag$BC, UB = dat[[1]]$tag$UB, GE = dat[[1]]$tag$GE)
    }
    return(dt)
  }, mc.cores = cores)
  rsamtools_reads <- rbindlist(rsamtools_reads, fill = TRUE, use.names = TRUE)
  
  #downsampling
  if (downsampling!=FALSE){
    if (nrow(rsamtools_reads) > downsampling){
      rsamtools_reads <- rsamtools_reads[sample(.N, downsampling)]
    }
  }
  # return( rsamtools_reads)
  if(inex){
    rsamtools_reads[ , ftype :="NA"][
      is.na(GEin)==F, ftype :="intron"][
        is.na(GE)==F  , ftype:="exon"][
          is.na(GE)     , GE:=GEin][
            ,GEin:=NULL ]
  }else{
    rsamtools_reads[, ftype :="NA"][
      is.na(GE)==F, ftype :="exon"]
  }
  setkey(rsamtools_reads,RG)
  
  # if(keepUnassigned){
  #   return( rsamtools_reads )
  # }else{
  #   return( rsamtools_reads[GE!="NA"] )
  # }
}

splitRG<-function(bccount,mem,hamdist){
  if(is.null(mem) || mem==0){
    maxR<- Inf
  }else{
    maxR<- floor( mem*1000 * 4500 )
  }
  #if( (maxR > 2e+09 & opt$read_layout == "SE") | (maxR > 1e+09 & opt$read_layout == "PE") ){
  #  maxR <- ifelse(opt$read_layout == "SE",2e+09,1e+09)
  #}
  if(maxR > 2e+09){
    maxR <- 2e+09
  }
  if(hamdist>0){ #multicore hamming distance takes a lot of memory
    #ram_factor <- ifelse(opt$num_threads>10, 5, 2)
    ram_factor <- 3
    maxR <- floor( maxR/ram_factor )
  }
  
  print(paste(maxR,"Reads per chunk"))
  nc<-nrow(bccount)
  cs=0
  chunkID=1
  bccount[,chunkID:=1]
  if(sum(bccount$n) > maxR) {
    for(i in 1:nc){
      cs=cs+bccount[i]$n
      if(bccount[i]$n>maxR){
        print(paste("Warning: Barcode",bccount[i]$XC,"has more reads than allowed for the memory limit!
                  Proceeding anyway..."))
      }
      if(cs>=maxR){
        if(i > 1){ #if the first BC exceeds the limit, keep chunkID 1
          chunkID=chunkID+1
        }
        cs=bccount[i][,"n"]
      }
      bccount[i][,"chunkID"]=chunkID
    }
  }
  return(bccount)
}

# set samples
samples <- c("Maxima_42", "Maxima_50", "SS_42", "SS_50")

# set threads
setDTthreads(threads=30)

complexity <- data.frame(RG=character(),
                         n=integer(),
                         project=character(),
                         replicate=integer(),
                         ds_level=integer())

# set ds-level & files, then run
for (k in c(14000000, 5000000, 1000000)){
  print(k)
  for (i in samples){
    print(paste0("sample ", i))
    for (j in 1:3){
      inputBAM = paste0("/data/share/htp/Felix_genotyping/primeSeq_tests/RTase_test/trimmed/",i,"/",i,".filtered.Aligned.GeneTagged.sorted.bam")
      bccount = paste0("/data/share/htp/Felix_genotyping/primeSeq_tests/RTase_test/trimmed/",i,"/zUMIs_output/",i,"kept_barcodes_binned.txt")
      bccount <- fread(bccount)
      bccount<-splitRG(bccount=bccount, mem= 50, hamdist = 0)
      reads <- reads2genes_new_ds(featfile = inputBAM,
                                  bccount  = bccount,
                                  inex     = TRUE,
                                  chunk    = 1,
                                  cores    = 20,
                                  downsampling = k)
      print(nrow(reads))
      reads <- reads %>% filter(!is.na(GE)) %>% dplyr::select(-UB, -ftype) %>% distinct() %>% dplyr::count(RG) %>% 
        mutate(project=i) %>%
        mutate(replicate=j) %>%
        mutate(ds_level=k)
      complexity <- rbind(complexity, reads)
      rm(reads)
    }
  }
}

saveRDS(complexity,
        here("data/RTase/complexity_RTase.rds"))