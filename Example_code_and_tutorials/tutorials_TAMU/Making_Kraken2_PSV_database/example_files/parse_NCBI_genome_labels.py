import os
import pandas as pd
import gzip
import argparse

def read_fasta(file_path):
    open_func = gzip.open if file_path.endswith('.gz') else open
    with open_func(file_path, 'rt') as f:
        header, sequence = '', ''
        for line in f:
            if line.startswith(">"):
                if header:
                    yield header, sequence
                header, sequence = line.strip(), ''
            else:
                sequence += line.strip()
        if header:
            yield header, sequence

def write_fasta(file_path, sequences):
    open_func = gzip.open if file_path.endswith('.gz') else open
    with open_func(file_path, 'wt') as f:
        for header, sequence in sequences:
            f.write(f"{header}\n")
            f.write(f"{sequence}\n")

def main(metadata_file, genome_folder, output_folder):
    df = pd.read_csv(metadata_file, delimiter='\t')

    name_taxid_map = {}
    sequence_number = 1

    for _, row in df.iterrows():
        local_filename = os.path.basename(row['local_filename'])
        if local_filename.endswith('.fna.gz'):
            local_filename = local_filename.rsplit('.fna.gz', 1)[0]
        taxid = row['taxid']
        organism_name = row['organism_name']
        name_taxid_map[local_filename] = (taxid, organism_name)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for file_name in os.listdir(genome_folder):
        if file_name.endswith('.fna.gz') or file_name.endswith('.fna'):
            original_name = file_name.rsplit('.fna', 1)[0]
        elif file_name.endswith('.fasta.gz') or file_name.endswith('.fasta'):
            original_name = file_name.rsplit('.fasta', 1)[0]
        else:
            continue

        if original_name in name_taxid_map:
            taxid, organism_name = name_taxid_map[original_name]
            original_file_path = os.path.join(genome_folder, file_name)
            modified_sequences = []

            for _, sequence in read_fasta(original_file_path):
                new_header = f">sequence{sequence_number}|kraken:taxid|{taxid}|{organism_name}"
                sequence_number += 1
                modified_sequences.append((new_header, sequence))

            modified_file_path = os.path.join(output_folder, f"modified_{file_name}")
            write_fasta(modified_file_path, modified_sequences)
        else:
            print(f"Warning: {original_name} not found in metadata.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Modify genome file headers based on metadata and save to output folder.')
    parser.add_argument('metadata_file', type=str, help='Path to the metadata file.')
    parser.add_argument('genome_folder', type=str, help='Path to the folder containing genome files.')
    parser.add_argument('output_folder', type=str, help='Path to the folder to save modified genome files.')
    args = parser.parse_args()

    main(args.metadata_file, args.genome_folder, args.output_folder)
