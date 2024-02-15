#!/bin/bash

# Example command:
# bash script_name.sh /path/to/dir1 /path/to/dir2 /path/to/dir3 ...

# Directories containing FASTQ files as command-line arguments
directories=("$@")

# Temporary associative array to hold unique sample identifiers
declare -A samples

# Find all FASTQ files and populate the samples associative array
# with unique sample identifiers
for dir in "${directories[@]}"; do
  fastq_files=$(find "$dir" -name '*.fastq.gz')
  for fq in $fastq_files; do
    base_name=$(basename "$fq")
    sample_id=$(echo "$base_name" | cut -d'_' -f1-1) # Here is where we split the file name to only include everything to the left of the first "_"
    samples["$sample_id"]=1
  done
done

# Directory to store concatenated files
mkdir -p cat_reads

# Loop over each unique sample identifier and concatenate files
for sample in "${!samples[@]}"; do
  # Concatenate the R1 files
  find "${directories[@]}" -name "${sample}_*_R1_*.fastq.gz" -exec cat {} + > cat_reads/${sample}_R1_cat.fastq.gz
  # Concatenate the R2 files
  find "${directories[@]}" -name "${sample}_*_R2_*.fastq.gz" -exec cat {} + > cat_reads/${sample}_R2_cat.fastq.gz
done
