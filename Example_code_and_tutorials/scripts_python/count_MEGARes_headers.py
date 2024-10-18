import argparse
from collections import defaultdict

# Example command:
# python count_MEGARes_headers.py 1


def read_fasta(file):
    sequences = {}
    with open(file, 'r') as f:
        sequence = ""
        header = ""
        for line in f:
            if line.startswith('>'):
                if header:
                    sequences[header] = sequence
                header = line.strip()
                sequence = ""
            else:
                sequence += line.strip()
        sequences[header] = sequence
    return sequences

def count_unique_sections(headers, index):
    unique_sections = set()
    for header in headers:
        section = header.split("|")[index]  
        unique_sections.add(section)
    return len(unique_sections)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Count unique sections in FASTA headers.')
    parser.add_argument('fasta_file', help='Path to the input FASTA file')
    parser.add_argument('index', type=int, help='Index of the section to consider (0-based)')

    args = parser.parse_args()

    sequences = read_fasta(args.fasta_file)
    total_unique_sections = count_unique_sections(sequences.keys(), args.index)

    print(f"Total unique sections for index {args.index}: {total_unique_sections}")
