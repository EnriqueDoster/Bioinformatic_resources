import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import product

# usage: python script.py input.fasta output.fasta

def replace_n(sequence, replacements):
    positions = [i for i, x in enumerate(sequence) if x == "N"]
    for combination in product(replacements, repeat=len(positions)):
        sequence_list = list(sequence)
        for position, replacement in zip(positions, combination):
            sequence_list[position] = replacement
        yield "".join(sequence_list)

def main(input_file, output_file):
    with open(input_file, "r") as input_handle, open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            if "N" in record.seq:
                for i, new_seq in enumerate(replace_n(str(record.seq), ["A", "T", "C", "G"]), 1):
                    new_seq = Seq(new_seq)  # Convert the string to a Seq object
                    new_record = SeqRecord(new_seq, id=f"{record.id}_{i}", description="")
                    # Write the new record to the output file as a single-line FASTA
                    output_handle.write(f">{new_record.id}\n{new_record.seq}\n")
            else:
                # Write the original record to the output file as a single-line FASTA
                output_handle.write(f">{record.id}\n{record.seq}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Replace N in fasta sequences.')
    parser.add_argument('input_file', type=str, help='Input fasta file')
    parser.add_argument('output_file', type=str, help='Output fasta file')
    args = parser.parse_args()

    main(args.input_file, args.output_file)
