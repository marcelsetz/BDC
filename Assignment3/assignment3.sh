#!/bin/bash

# Define the list of Fastq files to process
FASTQ_FILES=("/commons/Themas/Thema12/HPC/rnaseq.fastq")

# Function to process each Fastq file
process_fastq() {
    fastq_file=$1
    output_file="output.csv"
    python3 assignment3.py --chunk "$fastq_file" > "$output_file"
}

export -f process_fastq

# Run GNU Parallel to process the files
parallel process_fastq ::: "${FASTQ_FILES[@]}"