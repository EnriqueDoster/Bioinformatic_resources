

# 1) subset all samples such that you only keep samples with A, B, C time points
Habitating_SNV_data_noSNP.css <- subset_samples(SNV_data_noSNP.css, pool_timepoint %in% c("A","B","C"))

# 2) prune the subset samples to erase any 0s
Habitating_SNV_data_noSNP.css <- prune_taxa(taxa_sums(Habitating_SNV_data_noSNP.css) > 0, Habitating_SNV_data_noSNP.css)

# 3) Calculate beta diversity distance matrix
Habitating_SNV_data_noSNP.css.df <- as(sample_data(Habitating_SNV_data_noSNP.css),"data.frame")
Habitating_SNV_data_noSNP.css.dist <- vegdist(t(otu_table(Habitating_SNV_data_noSNP.css)), method = "bray")
Habitating_SNV_data_noSNP.css.ord <- vegan::metaMDS(comm = t(otu_table(Habitating_SNV_data_noSNP.css)), distance = "bray", try = 10, trymax = 100, autotransform = F)

#View Data frame to ensure that only timepoints A,B,C are included (looks good)

# 4) subset the distances for each pair of participant Id's

## Indices of samples from shared households

#103-105 Shared_2
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "103" | Habitating_SNV_data_noSNP.css.df$participant_id == "105") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5415667


#117-118 Shared_1
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "117" | Habitating_SNV_data_noSNP.css.df$participant_id == "118") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.4830927


#130-126 Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "130" | Habitating_SNV_data_noSNP.css.df$participant_id == "126") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.3596901



#113-121 Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "113" | Habitating_SNV_data_noSNP.css.df$participant_id == "121") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.600346


#108-129 Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "108" | Habitating_SNV_data_noSNP.css.df$participant_id == "129") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5424946



#115-135 Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "115" | Habitating_SNV_data_noSNP.css.df$participant_id == "135") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5286628



#101-120 Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "101" | Habitating_SNV_data_noSNP.css.df$participant_id == "120") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.29208



#122-132 Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "122" | Habitating_SNV_data_noSNP.css.df$participant_id == "132") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5695767


#104-102 Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "104" | Habitating_SNV_data_noSNP.css.df$participant_id == "102") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5477212


#111-110 Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "111" | Habitating_SNV_data_noSNP.css.df$participant_id == "110") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.3810635


#106-123 Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "106" | Habitating_SNV_data_noSNP.css.df$participant_id == "123") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.4002308


#128-124 Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "128" | Habitating_SNV_data_noSNP.css.df$participant_id == "124") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5957735



#114-133 Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "114" | Habitating_SNV_data_noSNP.css.df$participant_id == "133") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.525046



#Now compare Individuals from shared households to controls



#103-128 Shared_Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "103" | Habitating_SNV_data_noSNP.css.df$participant_id == "128") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5396353


#103-114 Shared_Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "103" | Habitating_SNV_data_noSNP.css.df$participant_id == "114") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5905242


#103-126 Shared_Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "103" | Habitating_SNV_data_noSNP.css.df$participant_id == "126") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5470155


#103-121 Shared_Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "103" | Habitating_SNV_data_noSNP.css.df$participant_id == "121") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.6372923


#105-128 Shared_Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "105" | Habitating_SNV_data_noSNP.css.df$participant_id == "128") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.4474943



#105-114 Shared_Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "105" | Habitating_SNV_data_noSNP.css.df$participant_id == "114") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.6289834


#105-126 Shared_Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "105" | Habitating_SNV_data_noSNP.css.df$participant_id == "126") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5415869


#105-121 Shared_Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "105" | Habitating_SNV_data_noSNP.css.df$participant_id == "121") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5852793


#117-128 Shared_Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "117" | Habitating_SNV_data_noSNP.css.df$participant_id == "128") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.6043026


#117-114 Shared_Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "117" | Habitating_SNV_data_noSNP.css.df$participant_id == "114") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.4160128


#117-126 Shared_Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "117" | Habitating_SNV_data_noSNP.css.df$participant_id == "126") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.4047874


#117-121 Shared_Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "117" | Habitating_SNV_data_noSNP.css.df$participant_id == "121") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5980149


#118-128 Shared_Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "118" | Habitating_SNV_data_noSNP.css.df$participant_id == "128") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.6357109


#118-114 Shared_Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "118" | Habitating_SNV_data_noSNP.css.df$participant_id == "114") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5042197


#118-126 Shared_Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "118" | Habitating_SNV_data_noSNP.css.df$participant_id == "126") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5156977

#118-121 Shared_Control
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$participant_id == "118" | Habitating_SNV_data_noSNP.css.df$participant_id == "121") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.6546043
