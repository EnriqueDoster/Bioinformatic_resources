import os
import pandas as pd
from Bio import SeqIO
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

# Load FastANI output
fastani_output = "fastani_output.txt"
df = pd.read_csv(fastani_output, sep='\t', header=None, names=["query", "ref", "ani", "matches", "total"])

# Filter out comparisons below a certain ANI threshold (e.g., 95%)
threshold = 98
filtered_df = df[df["ani"] >= threshold]

# Create a distance matrix
distance_matrix = filtered_df.pivot(index="query", columns="ref", values="ani").fillna(0)
distance_matrix += distance_matrix.T
distance_matrix = 100 - distance_matrix

# Set the diagonal of the distance matrix to 0
for seq_id in distance_matrix.index:
    distance_matrix.at[seq_id, seq_id] = 0

# Convert distance matrix to condensed format
condensed_distance_matrix = squareform(distance_matrix)

# Set any negative distance to 0
condensed_distance_matrix[condensed_distance_matrix < 0] = 0

# Perform hierarchical clustering
linked = linkage(condensed_distance_matrix, method="average")

# Assign clusters based on a distance threshold
distance_threshold = 100 - threshold
cluster_labels = fcluster(linked, distance_threshold, criterion="distance")

# Create a dictionary of clusters
clusters = {}
for label, query in zip(cluster_labels, distance_matrix.index):
    if label not in clusters:
        clusters[label] = set()
    clusters[label].add(query)

# Define the output directory where individual FASTA files are stored
output_dir = "individual_fastas"

# Extract representative sequences from each cluster
representative_sequences = []
for cluster in clusters.values():
    representative_seq_id = min(cluster)
    # Remove the output_dir and ".fasta" extension from the representative_seq_id before reading the file
    representative_seq_id = representative_seq_id.replace(f"{output_dir}/", "").replace(".fasta", "")
    representative_seq_file = os.path.join(output_dir, f"{representative_seq_id}.fasta")
    representative_seq = SeqIO.read(representative_seq_file, "fasta")
    representative_sequences.append(representative_seq)

# Save representative sequences to a multi-FASTA file
output_fasta = "meg_ANI98_rep_sequences.fasta"
with open(output_fasta, "w") as output_handle:
    SeqIO.write(representative_sequences, output_handle, "fasta")

