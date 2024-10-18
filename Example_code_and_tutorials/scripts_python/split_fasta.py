from Bio import SeqIO
import os

input_fasta = "../megares_database_v3.00.fasta"
output_dir = "individual_fastas"

os.makedirs(output_dir, exist_ok=True)

with open(input_fasta, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        # Split the ID by "|" and use the first part
        filename = record.id.split("|")[0]
        output_file = os.path.join(output_dir, f"{filename}.fasta")
        with open(output_file, "w") as output_handle:
            SeqIO.write(record, output_handle, "fasta")

