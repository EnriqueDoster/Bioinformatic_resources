import os
import glob
import random
import numpy as np
from math import ceil

# Make genome subsets from two directories, splits into 80% that go into kraken db and 20% which are subset into 150nt length reads.
# Creates a bash script will all commands, you need to run the resulting script by calling it with an sbatch script.

# Configuration
input_directory = '/scratch/user/enriquedoster/test_PSV_db_subset/Labeled_PSV_Mh_genomes/'
pasteur_genomes_directory = '/scratch/user/enriquedoster/clean_Mh_database/genome_folders/Labeled_Pasteur_genbank_genomes'
template_db = '/scratch/group/vero_research/databases/kraken2/Pasteurellaceae_Mh_PSV_db_benchmarking_template'
runs = 2  # Number of runs
start_confidence = 0.0
end_confidence = 1.0
confidence_step = 0.1
output_sbatch_filename = 'run_kraken_commands.sh'
temp_root_dir = '/scratch/user/enriquedoster/test_PSV_db_subset/Temp_dir'  # Root for temporary directories

# Function to select files and return selected and remaining files
def select_files(directory, selection_ratio=0.2):
    all_files = glob.glob(os.path.join(directory, '*.*a'))
    selected = random.sample(all_files, ceil(len(all_files) * selection_ratio))
    remaining = list(set(all_files) - set(selected))
    return selected, remaining

# Function to create sbatch commands for subsetted reads
def create_sbatch_commands_for_subset_reads(selected_files, temp_dir, run_db_path, sbatch_file):
    subsetted_read_paths = []  # Store paths of subsetted reads
    for file in selected_files:
        base_filename = os.path.basename(file)
        output_filename = f"{base_filename.rsplit('.', 1)[0]}.reads.{base_filename.rsplit('.', 1)[1]}"
        output_file_path = os.path.join(temp_dir, output_filename)
        subsetted_read_paths.append(output_file_path)
        sbatch_file.write(f"python make_genome_reads.py {file} 150 {output_file_path}\n")
    return subsetted_read_paths

# Prepare sbatch file
with open(output_sbatch_filename, 'w') as sbatch_file:
    #sbatch_file.write("#!/bin/bash\n")
    #sbatch_file.write("#SBATCH -J build_db -o build_db.out -t 24:00:00 --nodes=1 --ntasks-per-node=6 --mem=40G\n\n")
    
    for run_number in range(1, runs + 1):
        run_db_name = f"Run_{run_number}_Pasteurellaceae_Mh_PSV_db_benchmarking"
        run_db_path = f"/scratch/group/vero_research/databases/kraken2/{run_db_name}"
        temp_dir = os.path.join(temp_root_dir, f"Run_{run_number}")
        sbatch_file.write(f"mkdir -p {temp_dir}\n")
        sbatch_file.write(f"cp -r {template_db}/* {run_db_path}\n")
        
        # File selection
        main_selected_files, main_remaining_files = select_files(input_directory)
        pasteur_selected_files, pasteur_remaining_files = select_files(pasteur_genomes_directory)

        # Create subset reads and get their paths
        main_subsetted_read_paths = create_sbatch_commands_for_subset_reads(main_selected_files , temp_dir, run_db_path, sbatch_file)
        
        pasteur_subsetted_read_paths = create_sbatch_commands_for_subset_reads( pasteur_selected_files, temp_dir, run_db_path, sbatch_file)
        
        print(f"Adding {len(main_remaining_files)} Mh genomes and {len(pasteur_remaining_files)} Pasteurellaceae genomes to the DB")
        # Kraken build for remaining files
        for file in main_remaining_files + pasteur_remaining_files:
            sbatch_file.write(f"kraken2-build --add-to-library {file} --db {run_db_path}\n")
        
        sbatch_file.write(f"kraken2-build --build --db {run_db_path}\n")

        # Classification with varying confidence levels
        confidences = [round(start_confidence + x * confidence_step, 1) for x in range(int((end_confidence - start_confidence) / confidence_step) + 1)]
        
        print(f"Running classification on {len(main_subsetted_read_paths)} Mh sample reads and {len(pasteur_subsetted_read_paths)} Pasteurellaceae sample reads")
        # Run classification on subsetted reads
        for confidence in confidences:
            for read_path in main_subsetted_read_paths + pasteur_subsetted_read_paths:
                output_basename = os.path.basename(read_path).rsplit('.', 1)[0]
                output_file = f"{temp_dir}/{output_basename}.output_conf_{confidence}_run_{run_number}.txt"
                report_file = f"{temp_dir}/{output_basename}.report_conf_{confidence}_run_{run_number}.txt"
                sbatch_file.write(f"kraken2 --db {run_db_path} --confidence {confidence} --report {report_file} --output {output_file}\n")