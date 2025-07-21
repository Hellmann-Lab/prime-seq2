#!/bin/bash

#get absolute location of this script:
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

#determine the project folder
project_folder=`dirname ${SCRIPT_DIR}`

#get the flow cell number from a fastq file:
flowcell_name=`zcat ${project_folder}/00_fastq/*.gz| head -n 1`

echo "project name: $project_folder" > Flowcell_name.txt
echo "flow cell name: $flowcell_name" >> Flowcell_name.txt