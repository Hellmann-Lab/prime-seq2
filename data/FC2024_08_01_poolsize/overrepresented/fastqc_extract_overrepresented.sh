#!/bin/bash

for names in 80 320 920
do 
unzip /data/share/htp/prime-seq_NextGen/data/FC2024_08_01_poolsize/01_fastqc/lane1_primeseq_poolsize_${names}_r2_fastqc.zip
sed -n '/#Sequence/,$p' lane1_primeseq_poolsize_${names}_r2_fastqc/fastqc_data.txt | sed '/>>END_MODULE/q' | sed '/>>END_MODULE/d'  | sed 's/^#//' > ${names}.out.txt
rm -r lane1_primeseq_poolsize_${names}_r2_fastqc
done