import os
import glob
import random
import numpy as np
from math import ceil

## Python script that takes as input a directory of genomes, subsets 20% of them and makes commands to subset reads to run against
# kraken databases made with the remaining 80% of genomes. 

# Configuration
input_directory = '/scratch/user/enriquedoster/test_PSV_db_subset/Labeled_PSV_Mh_genomes'
template_db = '/scratch/group/vero_research/databases/kraken2/Pasteurellaceae_Mh_PSV_db_benchmarking_template'
runs = 1  # Number of runs
start_confidence = 0.0
end_confidence = 1.0
confidence_step = 0.1
output_sbatch_filename = 'run_kraken_analysis.sbatch'
temp_root_dir = '/scratch/user/enriquedoster/test_PSV_db_subset/Temp_dir'  # Root for temporary directories

# Prepare sbatch file
with open(output_sbatch_filename, 'w') as sbatch_file:
    sbatch_file.write("#!/bin/bash\n")
    sbatch_file.write("#SBATCH -J build_db -o build_db.out -t 24:00:00 --nodes=1 --ntasks-per-node=6 --mem=40G\n\n")
    
    # Copy template DB and set up runs
    for run_number in range(1, runs + 1):
        run_db_name = f"Run_{run_number}_Pasteurellaceae_Mh_PSV_db_benchmarking"
        run_db_path = f"/scratch/group/vero_research/databases/kraken2/{run_db_name}"
        temp_dir = os.path.join(temp_root_dir, f"Run_{run_number}")
        sbatch_file.write(f"mkdir -p {temp_dir}\n")
        sbatch_file.write(f"cp -r {template_db}/* {run_db_path}\n")
        
        # File selection and command generation
        # Adjusted file pattern to match '.fasta', '.fna', or '.fa' files
        file_patterns = [os.path.join(input_directory, pattern) for pattern in ('*.fasta', '*.fna', '*.fa')]
        all_files = [file for pattern in file_patterns for file in glob.glob(pattern)]
        
        artificial_read_file_paths = []
        
        selected_files = random.sample(all_files, ceil(len(all_files) * 0.2))
        for file in selected_files:
            # Adjust output filename to work with different extensions and add "_reads"
            base_filename = os.path.basename(file)
            output_filename = base_filename.rsplit('.', 1)[0] + '.reads.' + base_filename.rsplit('.', 1)[1]
            output_file_path = os.path.join(temp_dir, output_filename)
            artificial_read_file_paths.append(output_file_path)
            sbatch_file.write(f"python make_genome_reads.py {file} 150 {output_file_path}\n")
        
        # Kraken build for non-selected files
        for file in set(all_files) - set(selected_files):
            sbatch_file.write(f"kraken2-build --add-to-library {file} --db {run_db_path}\n")
        
        # Add genomes to Kraken DB and build
        sbatch_file.write("for file in /scratch/user/enriquedoster/clean_Mh_database/genome_folders/Labeled_Pasteur_genbank_genomes/*a\ndo\n")
        sbatch_file.write(f"    kraken2-build --add-to-library $file --db {run_db_path}\ndone\n")
        sbatch_file.write(f"kraken2-build --build --db {run_db_path}\n")

        # Classification with varying confidence levels
        # Adjust the generation of confidence levels to avoid floating-point issues
        confidences = [start_confidence + x * confidence_step for x in range(int((end_confidence - start_confidence) / confidence_step) + 1)]
        
        # Loop through each file in artificial_read_file_paths for Kraken classification at each confidence level
        for confidence in confidences:
            for file in artificial_read_file_paths:
                # Adjust to use the output file path with confidence level in the filename
                output_basename = os.path.basename(file).rsplit('.', 1)[0]
                output_file = f"{output_basename}.output_conf_{confidence}_{run_number}.txt"
                report_file = f"{output_basename}.report_conf_{confidence}_{run_number}.txt"
                output_file_path = os.path.join(temp_dir, output_file)
                report_file_path = os.path.join(temp_dir, report_file)
                sbatch_file.write(f"kraken2 --db {run_db_path} --confidence {confidence} --report {report_file_path} --output {output_file_path} {file}\n")

    # Placeholder for final steps to merge results

print(f"Generated sbatch script: {output_sbatch_filename}")
