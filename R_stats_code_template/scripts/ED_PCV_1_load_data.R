#### Mh++ ####

# Read in metadata
mhpp_map <- read.csv("data/metadata_mhpp.csv", header = T, row.names = 1)

# Read in count matrix
mhpp_otu.df <- read.csv("data/PSV_Results/strain_classification_kraken_matrix.csv", header = T, row.names = 1)
mhpp_otu <- otu_table(mhpp_otu.df, taxa_are_rows = T)

# Read in tax table
mhpp_taxa.df <- read.csv("data/PSV_Results/taxonomy_strain_classification_kraken_matrix.csv", header = T, row.names = 1)

# Make PS object
mhpp.ps <- merge_phyloseq(sample_data(mhpp_map), tax_table(as.matrix(mhpp_taxa.df)), mhpp_otu)
mhpp.ps

taxa_names <- rownames(tax_table(mhpp.ps))

# Filter to just find PSVs
selected_taxa <- taxa_names[grepl("Mh_PSV", taxa_names)]

# Filter taxa names based on the presence of "Mh" or "Mannheimia haemolytica"
#selected_taxa <- taxa_names[grepl("Mh", taxa_names) | grepl("Mannheimia haemolytica", taxa_names)]

# Use prune_taxa to keep only the selected taxa in the phyloseq object
mhpp_PSV.ps <- prune_taxa(selected_taxa, mhpp.ps)
