#!/bin/bash

for names in old new PTO
do 
unzip /data/share/htp/prime-seq_NextGen/data/FC2024_05_02_PoP96_BA/01_fastqc/lane1_primeseq_PoP96_BA_${names}_r2_fastqc.zip
sed -n '/#Sequence/,$p' lane1_primeseq_PoP96_BA_${names}_r2_fastqc/fastqc_data.txt | sed '/>>END_MODULE/q' | sed '/>>END_MODULE/d'  | sed 's/^#//' > ${names}.out.txt
rm -r lane1_primeseq_PoP96_BA_${names}_r2_fastqc
done