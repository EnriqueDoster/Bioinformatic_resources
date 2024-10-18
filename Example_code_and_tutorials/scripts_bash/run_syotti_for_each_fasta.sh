#!/bin/bash

# Uses syotti and cd-hit

# Function is to run syotti on each fasta file in the input folder, place the output in the output_folder, concatenates all baits, runs cd-hit to remove redundant

# User-defined prefix for naming output folders and files
prefix="Mycoplasmopsis"  # Replace with your desired prefix

# Directory containing the files
input_folder="${prefix}_genbank_genomes/"
# Output directory for syotti outputs
output_folder="${prefix}_output_folder/"

# Create the output directories if they don't exist
mkdir -p "$output_folder"

# Concatenated output file
concatenated_output="${prefix}_concatenated_output.fna"

# Empty the concatenated output file if it already exists
> "$output_folder$concatenated_output"

# Loop through each file in the input folder
for file in "$input_folder"*a
do
    # Extract the base name of the file
    base_name=$(basename "$file")
    
    # Construct the output file name
    output_file="${output_folder}${prefix}_${base_name}_syotti_baits_l120_h40.fna"
    
    # Run syotti command on each file
    syotti design --bait-len 120 --hamming-distance 40 -s "$file" -r -o "$output_file"

    # Concatenate the output to the final file
    cat "$output_file" >> "$output_folder$concatenated_output"
done

# Run CD-HIT to remove sequence redundancies
cdhit -i "$output_folder$concatenated_output" -o "$output_folder/baits_${prefix}_cleaned_output.fna" -c 1 -n 5 -d 0

echo "Processing complete. All files concatenated and cleaned output stored in $output_folder"
