!/usr/bin/env Rscript

require(optparse)

############# optparse options ############# 
option_list = list(
  make_option(c("--annotation"), type="character", default=NULL, 
              help="Sample-Annotation file name. Needs to contain sample barcode (if >192samples: i7+i5+BC merged!) as first column and sample-names as second column!", metavar="character"),
  make_option(c("--outfolder"), type="character", default=NULL, 
              help="output folder in which count matrix csv-files are to be written", metavar="character"),
  make_option(c("--use_names"), type="character", default="no", 
              help="Whether to use sample names as column names for count matrices; valida values: 'yes' or 'no'", metavar="character"), #rather use logical/boolean?
  make_option(c("--project_folder"), type="character", default=NULL, 
              help="project_folder (i.e. folder containing the subfolders 00_fastq, 01_trimming, etc.", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
############# 



############# import sample-annotation ############# 

sample_anno_full<-read.table(opt$annotation,sep="\t",header=T,fill = T) #read sample-annotation file; use "fill" to avoid errors due to missing values

#check if sample-annotation contains NA values, if so produce output to let the user know:
if (any(is.na(sample_anno_full))){write(paste(Sys.time(), " Your sample-annotation file contains NA-values! You might want to recheck this!\nNote: Only rows containing non-NA values in SampleBC-column will be further processed!", sep=""), stderr())}
sample_anno_full=sample_anno_full[sample_anno_full[,1]!="",]  #reduce sample-annotation to rows for which the sampleBC is not empty (sanity check...)

############# 

#############  Check if sample names are unique ############# 
# !! if NOT unique, throw a warning if "use_names==no" or ERROR&EXIT if use_names==T!!!
if (opt$use_names=="yes" && length(sample_anno_full[,2])!=length(unique(sample_anno_full[,2]))){
  write(paste(Sys.time(), " ERROR: Your sample-names inside the sample-annotation file are NOT unique and hence cannot be used to substitue the sample-barcodes! Count matrix processing will be aborted now...", sep=""), stderr())
  quit()
}
if (opt$use_names=="no" && length(sample_anno_full[,2])!=length(unique(sample_anno_full[,2]))){
  write(paste(Sys.time(), " WARNING: Your sample-names inside the sample-annotation file are NOT unique! Re-check your sample-annotation before submission!", sep=""), stderr())
}
############# 


############# import dge object ############# 

# to find dge object hand over project folder in pipeline and search within expected subfolder
dge_file=list.files(path = paste(opt$project_folder,"/03_zUMIs/zUMIs_output/expression/",sep=""),pattern = ".*.dgecounts.rds$",full.names = T)

dgecounts=readRDS(dge_file)
dgecounts_UMI_ex=as.matrix(dgecounts$umicount$exon$all)
dgecounts_UMI_inex=as.matrix(dgecounts$umicount$inex$all)

############# 

############# filter count-table for expected barcodes ############# 

# separately for ex & inex table, in case they differ (which shouldnt be the case actually, but anyways - copy&pasting is just so cool... :D )

### exon UMI counts
dgecounts_UMI_ex<-dgecounts_UMI_ex[,colnames(dgecounts_UMI_ex)%in%sample_anno_full[,1]] # filter for cell barcodes present in sample-anno (BCs MUST be in first column!)
sample_anno_ex<-sample_anno_full[sample_anno_full[,1]%in%colnames(dgecounts_UMI_ex),]           # in case some BCs were not detected, we also filter the other way around
sample_anno_ex<-sample_anno_ex[match(colnames(dgecounts_UMI_ex),sample_anno_ex[,1]),]             # now expected&observed BCs should be the same, so we can re-arrange the sample-annotation to match the order of the countmatrix-columns

### if #columns of count-matrix does NOT match number of samples in (unfiltered) sample-annotation, throw a WARNING message!!
if (ncol(dgecounts_UMI_ex) != nrow(sample_anno_full)){
  write(paste(Sys.time(), " WARNING: Number of samples in count-matrix does not match number of samples in sample-annotation! Please re-check before submission!!", sep=""), stderr())
}


### inex UMI counts
dgecounts_UMI_inex<-dgecounts_UMI_inex[,colnames(dgecounts_UMI_inex)%in%sample_anno_full[,1]] # filter for cell barcodes present in sample-anno (BCs MUST be in first column!)
sample_anno_inex<-sample_anno_full[sample_anno_full[,1]%in%colnames(dgecounts_UMI_inex),]           # in case some BCs were not detected, we also filter the other way around
sample_anno_inex<-sample_anno_inex[match(colnames(dgecounts_UMI_inex),sample_anno_inex[,1]),]             # now expected&observed BCs should be the same, so we can re-arrange the sample-annotation to match the order of the countmatrix-columns

### if #columns of count-matrix does NOT match number of samples in (unfiltered) sample-annotation, throw a WARNING message!!
if (ncol(dgecounts_UMI_inex) != nrow(sample_anno_full)){
  write(paste(Sys.time(), " WARNING: Number of samples in count-matrix does not match number of samples in sample-annotation! Please re-check before submission!!", sep=""), stderr())
}

############# 

############# Optional: Change columns-names ############# 

# if wanted, column-names are converted to sample-names
if (opt$use_names == "yes"){
  base::colnames(dgecounts_UMI_ex)=sample_anno_inex[,1]
  base::colnames(dgecounts_UMI_inex)=sample_anno_inex[,1]
}
  
############# 

############# rownames to column in order to create a nice table --?include optional conversion of gene_ids to gene_names??############# 
# do using base R to not make further dependencies?
dgecounts_UMI_ex <- cbind(rownames(dgecounts_UMI_ex), data.frame(dgecounts_UMI_ex, row.names=NULL))
colnames(dgecounts_UMI_ex)[1]="gene_id"

dgecounts_UMI_inex <- cbind(rownames(dgecounts_UMI_inex), data.frame(dgecounts_UMI_inex, row.names=NULL))
colnames(dgecounts_UMI_inex)[1]="gene_id"
############# 



############# Output count matrices as csv (also include excel file-types??) ############# 

write.table(x = dgecounts_UMI_ex,file = paste0(opt$outfolder,"/CountMatrix_exonic_UMIs.txt"),sep=",",col.names=T,row.names=F,quote = F) # output csv
write.table(x = dgecounts_UMI_inex,file = paste0(opt$outfolder,"/CountMatrix_InEx_UMIs.txt"),sep=",",col.names=T,row.names=F,quote = F) # output csv

############# 
