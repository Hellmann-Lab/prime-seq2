#!/bin/bash
#conda activate cutadapt #SFB server uses conda environment to run cutadapt (avoiding PYTHON versioning issues...)

threads=4
SLURM_RAM=4

# trimming of polyA's at the end of R2 using cutadapt
# included minimal read-length for trimming (40 nt for now - 28nt for R1)!


SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" #get absolute location of this script:
project_folder=`dirname ${SCRIPT_DIR}` #determine the project folder - this can then be used to define the input and output paths accordingly!

#project_folder="/data/share/htp/prime-seq_NextGen/data/FC2024_03_01_old_new_PTO/"

#determine the ending of the file-names

#if test -f `ls ${project_folder}/00_fastq/*r1.fastq.gz`; #CASE1: demultiplexed using deML and including index reads!
if test -n "$(find ${project_folder}/00_fastq/ -maxdepth 1 -name '*r1.fastq.gz' -print -quit)"
then 
file_suffix1="r1.fastq.gz"
file_suffix2="r2.fastq.gz"
#elif test -f `ls ${project_folder}/00_fastq/*read1.fastq.gz`
elif test -n "$(find ${project_folder}/00_fastq/ -maxdepth 1 -name '*read1.fastq.gz' -print -quit)"
then
file_suffix1="read1.fastq.gz"
file_suffix2="read2.fastq.gz"
else
	echo "ERROR: No fastq files matching the known naming-patterns were found within the 00_fastq subfolder! Please re-check your fastq files!"
	exit
fi

# now start trimming for each sample:
for i in `ls ${project_folder}/00_fastq/*${file_suffix1} | xargs -n1 basename`; 
do
	#sample=`cut -d _ -f 1-4 <<< $i` #trim file-name down to sample-name
	#sample_name=${i%$file_suffix} #remove suffix and use this as prefix for output files
	sample_name=${i%${file_suffix1}} 
	r1_file=${sample_name}${file_suffix1}
	r2_file=${sample_name}${file_suffix2}
	#echo ${sample_name}
	sbatch \
--error="${project_folder}/02_trimming/${sample_name}cutadapt.%J.err" \
--output="${project_folder}/02_trimming/${sample_name}cutadapt.%J.out" \
--cpus-per-task=$threads \
--mem=$SLURM_RAM \
--wrap "cutadapt -A A{30} -o ${project_folder}/02_trimming/${sample_name}r1_trimmed.fastq.gz -p ${project_folder}/02_trimming/${sample_name}r2_trimmed.fastq.gz ${project_folder}/00_fastq/${r1_file} ${project_folder}/00_fastq/${r2_file} --discard-casava --minimum-length 28:40 -j ${threads}"
done
