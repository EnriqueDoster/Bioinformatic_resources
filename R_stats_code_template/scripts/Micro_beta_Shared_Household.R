# Micro beta diversity by shared household status
#############################################################################################
####################################   CSS TRANSFORM    ###################################
data_micro.ps <- prune_taxa(taxa_sums(data_micro.ps) > 0, data_micro.ps)
any(taxa_sums(data_micro.ps)==0) # QUADRUPLE CHECKING - nope good.

data_micro.ps.css <- phyloseq_transform_css(data_micro.ps, log = F)
data_micro.ps.css.df <- as(sample_data(data_micro.ps.css), "data.frame")


###
####
# Only subset "Habitating" ######
####
### 

#Subset samples by timepoint- A,B,C
Habitating_data_micro.ps.css <- subset_samples(data_micro.ps.css, pool_timepoint %in% c("A","B","C"))

Habitating_data_micro.ps.css <- prune_taxa(taxa_sums(Habitating_data_micro.ps.css) > 0, Habitating_data_micro.ps.css)



# 3) Calculate beta diversity distance matrix
Habitating_data_micro.ps.css.df <- as(sample_data(Habitating_data_micro.ps.css),"data.frame")
Habitating_data_micro.ps.css.dist <- vegdist(t(otu_table(Habitating_data_micro.ps.css)), method = "bray")
Habitating_data_micro.ps.css.ord <- vegan::metaMDS(comm = t(otu_table(Habitating_data_micro.ps.css)), distance = "bray", try = 10, trymax = 100, autotransform = F)

#View Data frame to ensure that only timepoints A,B,C are included (looks good)

# 4) subset the distances for each pair of participant Id's

## Indices of samples from shared households

#103-105 Shared_2
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "103" | Habitating_data_micro.ps.css.df$participant_id == "105") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5298269


#117-118 Shared_1
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "117" | Habitating_data_micro.ps.css.df$participant_id == "118") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.522391


#130-126 Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "130" | Habitating_data_micro.ps.css.df$participant_id == "126") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5529644



#113-121 Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "113" | Habitating_data_micro.ps.css.df$participant_id == "121") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5917375


#108-129 Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "108" | Habitating_data_micro.ps.css.df$participant_id == "129") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.4503922



#115-135 Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "115" | Habitating_data_micro.ps.css.df$participant_id == "135") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5253418



#101-120 Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "101" | Habitating_data_micro.ps.css.df$participant_id == "120") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.494616



#122-132 Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "122" | Habitating_data_micro.ps.css.df$participant_id == "132") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.4883694


#104-102 Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "104" | Habitating_data_micro.ps.css.df$participant_id == "102") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5431302


#111-110 Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "111" | Habitating_data_micro.ps.css.df$participant_id == "110") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5476034


#106-123 Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "106" | Habitating_data_micro.ps.css.df$participant_id == "123") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.493665


#128-124 Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "128" | Habitating_data_micro.ps.css.df$participant_id == "124") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5246865



#114-133 Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "114" | Habitating_data_micro.ps.css.df$participant_id == "133") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5419941



#Now compare Individuals from shared households to controls



#103-128 Shared_Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "103" | Habitating_data_micro.ps.css.df$participant_id == "128") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5662211


#103-114 Shared_Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "103" | Habitating_data_micro.ps.css.df$participant_id == "114") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.6125238


#103-126 Shared_Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "103" | Habitating_data_micro.ps.css.df$participant_id == "126") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5266252


#103-121 Shared_Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "103" | Habitating_data_micro.ps.css.df$participant_id == "121") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.6394345


#105-128 Shared_Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "105" | Habitating_data_micro.ps.css.df$participant_id == "128") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5390553



#105-114 Shared_Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "105" | Habitating_data_micro.ps.css.df$participant_id == "114") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.6011602


#105-126 Shared_Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "105" | Habitating_data_micro.ps.css.df$participant_id == "126") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5608525


#105-121 Shared_Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "105" | Habitating_data_micro.ps.css.df$participant_id == "121") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.6047697


#117-128 Shared_Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "117" | Habitating_data_micro.ps.css.df$participant_id == "128") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.4774625


#117-114 Shared_Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "117" | Habitating_data_micro.ps.css.df$participant_id == "114") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.521449



#117-126 Shared_Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "117" | Habitating_data_micro.ps.css.df$participant_id == "126") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5331161


#117-121 Shared_Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "117" | Habitating_data_micro.ps.css.df$participant_id == "121") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5060688


#118-128 Shared_Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "118" | Habitating_data_micro.ps.css.df$participant_id == "128") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5039478


#118-114 Shared_Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "118" | Habitating_data_micro.ps.css.df$participant_id == "114") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5802849


#118-126 Shared_Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "118" | Habitating_data_micro.ps.css.df$participant_id == "126") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.5306929

#118-121 Shared_Control
shared_indices <- which(Habitating_data_micro.ps.css.df$participant_id == "118" | Habitating_data_micro.ps.css.df$participant_id == "121") ## modify this function to pick up
## both participant Ids that you want to compare. use the "|" symbol
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_data_micro.ps.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])
#0.6012312
