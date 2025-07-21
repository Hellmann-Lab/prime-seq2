#!/bin/bash
set -euo pipefail

while IFS=$'\t' read -r PROJECT FASTQ_CDNA; do
  echo "Submitting job for $PROJECT with $FASTQ_CDNA"
  export PROJECT FASTQ_CDNA
  sbatch 00_quantify_miscounting.sh
  sleep 1 # avoid overloading scheduler
done < to_run.txt 