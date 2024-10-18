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

def write_dmp_files(output_folder, label_taxid_map,parent_taxid, label_prefix):
    with open(os.path.join(output_folder, 'names.dmp'), 'w') as names_f, \
         open(os.path.join(output_folder, 'nodes.dmp'), 'w') as nodes_f, \
         open(os.path.join(output_folder, 'accession2taxid'), 'w') as acc_f:
        for label, taxid in label_taxid_map.items():
            # Include the label prefix in the sanitized label
            sanitized_label = f"{label_prefix}{label.replace(' ', '_')}"
            names_f.write(f"{taxid}\t|\t{sanitized_label}\t|\t\t|\tscientific name\t|\n")
            nodes_f.write(f"{taxid}\t|\t{parent_taxid}\t|\tstrain\t|\t\t|\t0\t|\t1\t|\t1\t|\t1\t|\t0\t|\t1\t|\t1\t|\t0\t|\t\t|\n")
            acc_f.write(f"{sanitized_label}\t{sanitized_label}.1\t{taxid}\t0\n")

def main(metadata_file, genome_folder, output_folder, label_prefix, parent_taxid):
    df = pd.read_csv(metadata_file, delimiter='\t')

    name_taxid_map = {}
    label_taxid_map = {}
    current_taxid = 1000000000
    sequence_number = 1

    for label in df['ProposedLabel'].unique():
        if not pd.isna(label):
            label_taxid_map[label] = current_taxid
            current_taxid += 1

    for _, row in df.iterrows():
        local_filename = os.path.basename(row['SequenceName'])
        if local_filename.endswith('.fna.gz'):
            local_filename = local_filename.rsplit('.fna.gz', 1)[0]
        
        if row['ProposedLabel'] not in label_taxid_map:
            raise ValueError(f"ProposedLabel '{row['ProposedLabel']}' not found in label_taxid_map.")
        
        taxid = label_taxid_map[row['ProposedLabel']]
        organism_name = row['SequenceName']
        name_taxid_map[local_filename] = (taxid, row['ProposedLabel'], organism_name)

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
            taxid, proposed_label, organism_name = name_taxid_map[original_name]
            original_file_path = os.path.join(genome_folder, file_name)
            modified_sequences = []

            for _, sequence in read_fasta(original_file_path):
                new_header = f">sequence{sequence_number}|kraken:taxid|{taxid} {label_prefix}{proposed_label}|{organism_name}"
                sequence_number += 1
                modified_sequences.append((new_header, sequence))

            modified_file_path = os.path.join(output_folder, f"modified_{file_name}")
            write_fasta(modified_file_path, modified_sequences)
        else:
            print(f"Warning: {original_name} not found in metadata.")

    write_dmp_files(output_folder, label_taxid_map, parent_taxid, label_prefix )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script to modify genome file headers based on metadata for use with kraken2. '
                                                 'The metadata file must be a TSV file containing at least the columns '
                                                 '`SequenceName` and `ProposedLabel`. The genome folder should contain '
                                                 'FASTA files (.fna or .fasta, with optional .gz compression).')

    parser.add_argument('--metadata_file', type=str, required=True, help='Path to the metadata file (TSV format).')
    parser.add_argument('--genome_folder', type=str, required=True, help='Path to the folder containing genome FASTA files.')
    parser.add_argument('--output_folder', type=str, required=True, help='Path to the folder to save modified genome files.')
    parser.add_argument('--label_prefix', type=str, default='', help='Prefix to add before each label (default: none).')
    parser.add_argument('--parent_taxid', type=str, default='75985', help='TaxID of the parent species (default: 75985 for Mannheimia haemolytica).')

    args = parser.parse_args()

    main(args.metadata_file, args.genome_folder, args.output_folder, args.label_prefix, args.parent_taxid)
