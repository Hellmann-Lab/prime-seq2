To demultiplex your prime-seq to individual samples, you need to call the bash script located within this directory.
It needs the following parameters:

1st parameter: name of the project. will be used as prefix for the demultiplexed fastq-files.
2nd parameter: absolute path to sample annotation table (tabbed text file!). Needs to contain unique sample barcodes within first column and sample-names in second column (as in template-sample-annotation)
3rd parameter: use barcodes or samplenames (2nd column in annotation) for creation of file names? valid values: "BC" or "sample" (Note: sample names need to be unique if utilized!!)
4th parameter: utilize only sample-barcodes or i7+i5+BC for demultiplexing? valid values: "normal" or "extended" (Note: if using "extended" barcodes in sample annotation .....)
5th parameter: number of threads to utilize for cutadapt (trimming of sample-barcodes from r1 fastq)
6th parameter: additionally output raw count matrix from zUMIs? -- valid values: "yes", "no"

 
Note:
- the script needs to be located within your individual project folder, which also needs to contain the standardized folder-structure
- hence, the demultiplexed fastq files also need to be present inside the "00_fastq" subfolder within the project folder
- the sample annotation:
	- is best created using the sample annotation template :)
	- needs to be a tab-separated text file
	- needs to include the unique sample-barcode within the first column and unique sample-names within the second column (names can also be the same as the unique BC in column 1)
	- sample-barcodes need to be unique! if barcodes themselves are not unique, you need to provide the concatenated i7+i5+BC sequence in this column!
- if sample-names (second column within sample annotation file!) are used for file-names, they will also be used for export of the count matrices. Make sure they are unique!
- to ensure that everything worked as expected, check your output (also SLURM out-files) BEFORE submitting your data to GEO/SRA/ENA!

Examples:

Demultiplex samples within project folder only using the "normal" prime-seq sample-BC in the R1 fastq :
bash project_sample_demux.sh project_name /path/to/project_folder/sample_annotation.txt BC normal 10 no

Demultiplex samples within project folder using the extended BC (prime-seq sample-BC AND the i7+i5 index reads), replace the long BCs by sample-names & output count matrices:
bash project_sample_demux.sh project_name /path/to/project_folder/sample_annotation.txt sample extended 10 yes

