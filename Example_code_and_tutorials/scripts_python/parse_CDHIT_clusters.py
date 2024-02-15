from Bio import SeqIO

cluster_file = "meg_90id.fasta.clstr"
fasta_file = "megares_database_v3.00.fasta"
output_file = "meg_90id_rep_seqs.fasta"

# Parse the cluster file
clusters = []
with open(cluster_file, "r") as f:
    for line in f:
        if line.startswith(">"):
            cluster = []
            clusters.append(cluster)
        else:
            accession = line.split(">")[1].split("...")[0]
            cluster.append(accession)

# Create a dictionary for the input fasta sequences
seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# Extract the representative sequences
representative_sequences = []
for cluster in clusters:
    representative_accession = cluster[0]
    representative_seq = seq_dict[representative_accession]
    representative_sequences.append(representative_seq)

# Write the representative sequences to a fasta file
with open(output_file, "w") as f:
    SeqIO.write(representative_sequences, f, "fasta")

