#!/usr/bin/env Rscript
require(optparse)

########### READ options ###########
option_list = list(
  make_option(c("--bcfile"), type="character", default=NULL, 
              help="Sample-Annotation file name. Needs to contain sample barcode (if >192samples: i7+i5+BC merged!) as first column and sample-names as second column!", metavar="character"),
  make_option(c("--outfolder"), type="character", default=NULL, 
              help="output folder", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
###########



sample_anno_full<-read.table(opt$bcfile,sep="\t",header=T,fill = T) #read sample-annotation file; use "fill" to avoid errors due to missing values
#check if sample-annotation contains NA values, if so produce output to let the user know:
if (any(is.na(sample_anno_full))){write(paste(Sys.time(), " Your sample-annotation file contains NA-values! You might want to recheck this!\nNote: Only rows containing non-NA values in SampleBC-column will be further processed!", sep=""), stderr())}

sample_anno_full=sample_anno_full[sample_anno_full[,1]!="",]#reduce sample-annotation to rows for which the sampleBC is not empty (sanity check...)

sample_anno<-data.frame(Index=sample_anno_full[,1],Name=sample_anno_full[,2]) #create short table for demultiplexing

colnames(sample_anno)<-c("#Index","Name") #re-define colnames - 1st column needs a hashtag in front of its name!



write.table(sample_anno,paste0(opt$outfolder,"/index.txt"),sep="\t",col.names=T,row.names=F,quote = F) # output demultiplexing list -> input for deML

write.table(sample_anno_full,paste0(opt$outfolder,"/sample_annotation.txt"),sep="\t",col.names=T,row.names=F,quote = F) # put copy of full sample-annotation file into sample-demultiplexing folder
