import argparse
import re
import os

# Usage: take as input fasta files and remove sequences shorter than "-l"
# Example command:
# python script_name.py -i file1.fasta file2.fasta -l 100 -o /path/to/output/directory


def process_fasta(file_path, min_length, output_dir):
    with open(file_path, 'r') as file:
        content = file.read()

    # Regular expression to match fasta sequences
    fasta_pattern = r'(>.*?\n)([^>]+)'
    
    # Find all sequences
    sequences = re.findall(fasta_pattern, content, re.DOTALL)
    
    # Filter sequences based on length
    filtered_sequences = [seq for seq in sequences if len(seq[1].replace('\n', '')) >= min_length]
    
    # Combine back into fasta format
    edited_content = ''.join([f'{header}{sequence}' for header, sequence in filtered_sequences])
    
    # Determine new file path with 'trimmed_' prefix in the specified output directory
    base_name = os.path.basename(file_path)
    new_file_path = os.path.join(output_dir, f'trimmed_{base_name}')

    # Write to the new file
    with open(new_file_path, 'w') as file:
        file.write(edited_content)

def main():
    parser = argparse.ArgumentParser(description='Process fasta files to remove short sequences and save in a specified directory.')
    parser.add_argument('-i', '--input_files', nargs='+', required=True, help='List of fasta files to process')
    parser.add_argument('-l', '--min_length', type=int, required=True, help='Minimum sequence length to keep')
    parser.add_argument('-o', '--output_dir', required=True, help='Directory to save output files')
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    for file_path in args.input_files:
        process_fasta(file_path, args.min_length, args.output_dir)

if __name__ == '__main__':
    main()
