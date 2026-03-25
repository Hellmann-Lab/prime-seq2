#!/bin/bash
#SBATCH --job-name=count_fqs

output="read_counts_2.tsv"
echo -e "File\tReads" > "$output"

for file in /data/share/htp/prime-seq_NextGen/data/data_collection_for_submission/second_submission/demultiplexed_fastqs/*.fq.gz; do
    [[ -f "$file" ]] || continue

    lines=$(gzip -cd "$file" | wc -l)
    reads=$((lines / 4))
    name=$(basename "$file")

    echo -e "$name\t$reads" >> "$output"
done

