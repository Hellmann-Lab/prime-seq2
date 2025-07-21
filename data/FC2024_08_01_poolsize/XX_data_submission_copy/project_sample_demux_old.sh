#!/bin/bash
# usage: 
# This script assumes that it's stsarted within a standardized prime-seq processing folder!
# Hence, the "raw" fastq files should be present within the "00_fastq" folder and file-names should either end with "read1.fastq.gz" or "r1.fq.gz". (= gene-center demux-style or deML demux-style)
#
# 1st parameter: name of the project. will be used as prefix for the demultiplexed fastq-files.
# 2nd parameter: absolute path to sample annotation table (tabbed text file!). Needs to contain unique sample barcodes within first column and sample-names in second column (as in template-sample-annotation)
# 3rd parameter: use barcodes or samplenames (2nd column in annotation) for creation of file names? valid values: "BC" or "sample" (Note: sample names need to be unique if utilized!!)
# 4th parameter: utilize only sample-barcodes or i7+i5+BC for demultiplexing? valid values: "normal" or "extended" (Note: if using "extended" barcodes in sample annotation .....)
# 5th parameter: number of threads to utilize for cutadapt (trimming of sample-barcodes from r1 fastq)
# 6th parameter: additionally output raw count matrix from zUMIs? -- valid values: "yes", "no"
# 
# Example: bash project_sample_demux.sh project_name /path/to/project_folder/sample_annotation.txt BC extended 10 yes


name=PoP96_2_old
bcfile=/data/share/htp/prime-seq_NextGen/data/FC2024_06_01_PoP96_FP_BA/XX_data_submission/SampleAnnotation_old.txt
extendedBC=normal # extendedBC must be set to "normal" to use only sample-barcodes! otherwise, i7+i5+BC concatenation is assumed!!
BcOrName=BC
threads=10
OutputMatrix=no
# if using sample-names for file-names, we also use them within the count-matrix! To make that eas



#NOTE: we assume that both utilized scripts (catfq.pl and make_deML_index.R) are located within the correct subfolders within the project-folder!
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" #absolute path to where THIS shell script is located
project_folder=`dirname ${SCRIPT_DIR}` #project folder


# define input and output by standardized project-folder structure:
infolder=$project_folder/00_fastq
outfolder=$project_folder/XX_data_submission

#check which names the raw fastqs have:
if test -f `ls ${project_folder}/00_fastq/*r1.fq.gz`; #CASE1: demultiplexed using deML
then 
	bcumi=`ls ${project_folder}/00_fastq/*r1.fq.gz | xargs -n1 basename`
	cDNA=`ls ${project_folder}/00_fastq/*r2.fq.gz | xargs -n1 basename`
elif test -f `ls ${project_folder}/00_fastq/*read1.fastq.gz`; #CASE2: demultiplexed using deML but NOT including index reads!
then 
	bcumi=`ls ${project_folder}/00_fastq/read1*.fastq.gz | xargs -n1 basename`
	cDNA=`ls ${project_folder}/00_fastq/read2*.fastq.gz | xargs -n1 basename`
else
	echo "ERROR: No fastq files matching the known naming-patterns were found within the 00_fastq subfolder! Please re-check your fastq files!"
	exit
fi

#create new subfolder which will contain all demultiplexed fastqs
mkdir $outfolder/SubmissionData_$name 
outfolder=$outfolder/SubmissionData_$name

#touch $outfolder/message.out
#chmod 775 $outfolder/message.out

cd $outfolder

# if the last command line parameter (5th parameter; "extendedBCs") is missing, only the prime-seq sample-barcode will be utilized for demultiplexing ( and NOT the index reads)
if [ $extendedBC = "normal" ]
then
	#if no concatenation is necessary, we create a "fake job" -- so we dont need to change the script downstream regarding SLURM-JobIDs :D
	jid1=$(sbatch --wrap 'sleep 3s; echo "No extended barcodes usage specified, continuing with R1 as sole BC file"') # output info that only sample-barcode will be utilized (and shortly wait in order to have the job still in queue when submitting the next one...)
	j1=`echo $jid1 | cut -f4 -d " "` ## extract jobid from "Submitted batch job jobid" - we can use this to queue further jobs and have them started after this one finishes
elif [ $extendedBC = "extended" ]
then
	if test -f `ls ${project_folder}/00_fastq/*i1.fq.gz`; then
		i7=`ls ${project_folder}/00_fastq/*i1.fq.gz | xargs -n1 basename`
		i5=`ls ${project_folder}/00_fastq/*i2.fq.gz | xargs -n1 basename`
	elif test -f `ls ${project_folder}/00_fastq/*index1.fastq.gz`; then
		i7=`ls ${project_folder}/00_fastq/index1*.fastq.gz | xargs -n1 basename`
		i5=`ls ${project_folder}/00_fastq/index2*.fastq.gz | xargs -n1 basename`
	else
		echo "You wanted to demux based on i7+i5+BC, but no index-fastq's matching the expected names are present within the 00_fastq folder! Please re-check your fastq files!"
		exit
	fi
	jid1=$(sbatch -J cat.$name --wrap "perl $project_folder/XX_data_submission/scripts/catfq.pl $infolder/$i7 $infolder/$i5 $infolder/$bcumi $infolder/$name.i7i5BC_R1.fq $threads") #call the PERL script for merging the fastq files
	j1=`echo $jid1 | cut -f4 -d " "` ## extract jobid from "Submitted batch job jobid"

	bcumi=$name.i7i5BC_R1.fq.gz # replace file-name within bcumi-variable to match the merged (i7+i5+bcumi) file for further processing
else  # if last command-line parameter is provided i7, i5 and BC will be concatenated for demultiplexing
	echo 'ERROR: Please provide a valid value for the extendedBC-choice (parameter 5??). Valid values are "normal" or "extended".'
	exit
fi


# trim down read to only sample-BC in order to use it as "index read" for demultiplexing
# if using only BC we cut down to the first 16 nucleotides and if using i7+i5+BC we cut down to 8+8+16=32 nucleotides:
if [ $extendedBC = "normal" ]
then
	jid2=$(sbatch -J cut.$name --dependency afterany:$j1 --wrap "cutadapt -j $threads -u -16 -o $infolder/$name.BC.fq.gz $infolder/$bcumi") #get rid of the last 16 bases within original R1 using cutadapt (this is the UMI, which will not be utilized for demultiplexing)
	j2=`echo $jid2 | cut -f4 -d " "` ## extract jobid from "Submitted batch job jobid"
else
	jid2=$(sbatch -J cut.$name --dependency afterany:$j1 --wrap "cutadapt -j $threads -u -32 -o $infolder/$name.BC.fq.gz $infolder/$bcumi")
	j2=`echo $jid2 | cut -f4 -d " "` ## extract jobid from "Submitted batch job jobid"
fi


# create sample-index table for demultiplexing - copy either 2x 1st column or 1st and 2nd column of sample-annotation into new text file:
if [ $BcOrName = "BC" ];
then 
	#awk: check if first column is not empty (--broken annotation file) and if so, copy first column to new file
	awk -F'\t' 'BEGIN { OFS = FS } { if ($1 != "") print $1,$1}'  $bcfile > $outfolder/index.txt #no need to create a SLURM job for this tiny bit of work...
elif [ $BcOrName = "sample" ];
then
	awk -F'\t' 'BEGIN { OFS = FS } { if ($1 != "") print $1,$2}'  $bcfile > $outfolder/index.txt #no need to create a SLURM job for this tiny bit of work...
else
	echo 'Value for 3rd parameter (file-names) is invalid! Please provide either "BC" or "sample" as value!'
	exit
fi

# replace column names by hash-tagged header for deML compatibility:
sed -i "1s/.*/#index\tname/" $outfolder/index.txt 


# start the actual demultiplexing using deML:
jid3=$(sbatch -J deML.$name  --dependency afterok:$j2 --wrap "deML --summary $outfolder/${name}_deML_summary.txt --error $outfolder/${name}_deML_error.txt -i $outfolder/index.txt -if1 $infolder/${name}.BC.fq.gz -f $infolder/$bcumi -r $infolder/$cDNA -o $outfolder/$name") #start demultiplexing using deML
j3=`echo $jid3 | cut -f4 -d " "` ## extract jobid from "Submitted batch job jobid"


# delete any unwanted file output:
jid4=$(sbatch -J clean.$name --dependency afterok:$j3 --wrap "rm $infolder/${name}.i7i5BC_R1.fq.gz $infolder/${name}.BC.fq.gz $outfolder/index.txt $outfolder/*_unknown_*.fq.gz $outfolder/*fail.fq.gz") # cleaning step: remove unassigned read and produce a text-file containing the md5 hashes of all produced fastq files
j4=`echo $jid4 | cut -f4 -d " "` ## extract jobid from "Submitted batch job jobid"


# create a text file containing md5 hashes for the final fastq files (only after unwanted fastqs are deleted!) --this doesnt work yet; dependency never satisfied!; set to "afterany" as some files to be deleted may not have existed...
sbatch -J md5.$name --dependency afterany:$j4 --wrap "md5sum *.fq.gz >md5sums.txt" # produce a text-file containing the md5 hashes of all produced fastq files

# create count-matrices in csv format
if [ $OutputMatrix =  "yes" ]
then
sbatch -J CM.$name --wrap "Rscript ${project_folder}/XX_data_submission/scripts/create_count_table.R --annotation $bcfile --outfolder $outfolder --use_names=$BcOrName --project_folder $project_folder"
fi
