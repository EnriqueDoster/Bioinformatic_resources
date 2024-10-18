#!/bin/bash
#SBATCH -J cat_reads -o cat_log.out -t 12:00:00 -p knl --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=10G

# Bash script to concate sample reads sequenced on multiple lanes
# Must be run in the working directory with the fastq files.

# You'll likely have to change this next command to match your reads.
# Here, we take everything to the left of the first "_" to identify
# unique sample prefixes. 

# Get a list of unique sample identifiers
samples=$(ls *.fastq.gz | cut -d'_' -f1-1 | sort | uniq)

# Loop over each unique sample identifier
mkdir cat_reads

for sample in $samples
do
  # Concatenate the R1 files
  cat ${sample}_*_R1_*.fastq.gz > cat_reads/${sample}_R1_cat.fastq.gz
  # Concatenate the R2 files
  cat ${sample}_*_R2_*.fastq.gz > cat_reads/${sample}_R2_cat.fastq.gz
done