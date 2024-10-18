import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random
# picks a nucleotide for each "N" in a sequence instead of generating all possible combinations. 
# usage: python script.py input.fasta output.fasta

def replace_n_randomly(sequence, replacements):
    sequence_list = list(sequence)
    for i, nucleotide in enumerate(sequence_list):
        if nucleotide == "N":
            sequence_list[i] = random.choice(replacements)
    return "".join(sequence_list)

def main(input_file, output_file):
    with open(input_file, "r") as input_handle, open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            if "N" in record.seq:
                new_seq = replace_n_randomly(str(record.seq), ["A", "T", "C", "G"])
                new_seq = Seq(new_seq)  # Convert the string to a Seq object
                new_record = SeqRecord(new_seq, id=f"{record.id}", description="")
                # Write the new record to the output file as a single-line FASTA
                output_handle.write(f">{new_record.id}\n{new_record.seq}\n")
            else:
                # Write the original record to the output file as a single-line FASTA
                output_handle.write(f">{record.id}\n{record.seq}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Randomly replace N in fasta sequences.')
    parser.add_argument('input_file', type=str, help='Input fasta file')
    parser.add_argument('output_file', type=str, help='Output fasta file')
    args = parser.parse_args()

    main(args.input_file, args.output_file)
