#!/bin/bash

#conda activate cutadapt #SFB server uses conda environment to run cutadapt (avoiding PYTHON versioning issues...)

threads=16
SLURM_RAM=16

# trimming of polyA's at the end of R2 using cutadapt
# included minimal read-length for trimming (40 nt for now - 28nt for R1)!


SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" #get absolute location of this script:
project_folder=`dirname ${SCRIPT_DIR}` #determine the project folder - this can then be used to define the input and output paths accordingly!



#try to adapt input based on common naming patterns for fq files:
if test -f `ls ${project_folder}/00_fastq/*r1.fq.gz`; #CASE1: demultiplexed using deML
then 
	r1_file=`ls ${project_folder}/00_fastq/*r1.fq.gz | xargs -n1 basename`
	r2_file=`ls ${project_folder}/00_fastq/*r2.fq.gz | xargs -n1 basename`
	fastq_name=${r1_file%"r1.fq.gz"} #remove anything after r1.fq.gz and use this as prefix for output files
	
elif test -f `ls ${project_folder}/00_fastq/*read1.fastq.gz`; #CASE2: demultiplexed using deML but NOT including index reads!
then 
	r1_file=`ls ${project_folder}/00_fastq/read1*.fastq.gz | xargs -n1 basename`
	r2_file=`ls ${project_folder}/00_fastq/read2*.fastq.gz | xargs -n1 basename`
	fastq_name=${r1_file%"read1.fastq.gz"} #remove anything after read1.fastq.gz and use this as prefix for output files
	
else
	echo "ERROR: No fastq files matching the known naming-patterns were found within the 00_fastq subfolder! Please re-check your fastq files!"
	exit
fi

#we now use these variables to construct the SLURM job processing the data included within the current project folder:
sbatch \
--error="${SCRIPT_DIR}/cutadapt.%J.err" \
--output="${SCRIPT_DIR}/cutadapt.%J.out" \
--cpus-per-task=$threads \
--mem=$SLURM_RAM \
--wrap "cutadapt -A A{30} -o ${project_folder}/02_trimming/${fastq_name}r1_trimmed.fq.gz -p ${project_folder}/02_trimming/${fastq_name}r2_trimmed.fq.gz ${project_folder}/00_fastq/${r1_file} ${project_folder}/00_fastq/${r2_file} --discard-casava --minimum-length 28:40 -j ${threads}"

