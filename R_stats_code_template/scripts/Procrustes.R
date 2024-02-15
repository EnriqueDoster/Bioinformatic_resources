# Procrustes analysis

# Load required packages
library(phyloseq)
library(vegan)

sample_names(data_micro.ps) <- sample_data(data_micro.ps)$ShortID
sample_names(data_noSNP) <- sample_data(data_noSNP)$ShortID

# Find the shared sample names between both phyloseq objects
shared_samples <- intersect(sample_names(data_micro.ps), sample_names(data_noSNP))

# Subset both phyloseq objects to include only the shared samples
data_micro.ps_subset <- prune_samples(sample_names(data_micro.ps) %in% shared_samples, data_micro.ps)
data_noSNP_subset <- prune_samples(sample_names(data_noSNP) %in% shared_samples, data_noSNP)


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
data_micro_dist <- distance(data_micro.ps_subset, method = "bray")
data_noSNP_dist <- distance(data_noSNP_subset, method = "bray")

# Perform PCoA ordination using the vegan::capscale function
data_micro_pcoa <- capscale(data_micro_dist ~ 1)
data_noSNP_pcoa <- capscale(data_noSNP_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
data_micro_scores <- vegan::scores(data_micro_pcoa, display = "sites")
data_noSNP_scores <- vegan::scores(data_noSNP_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result <- protest(data_micro_scores, data_noSNP_scores)

# Print Procrustes result summary
print(procrustes_result)

# Optional: Plot the Procrustes analysis
plot(procrustes_result)


###
##### Now by subset samples #####
###



#### Start with Start_CONV ####
# Subset the phyloseq objects to keep only samples with "Start_CONV" in the "Combined_order" column
data_micro.ps_Start_CONV <- subset_samples(data_micro.ps_subset, Combined_order == "Start_CONV")
data_noSNP_Start_CONV <- subset_samples(data_noSNP_subset, Combined_order == "Start_CONV")

# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
data_micro_dist_Start_CONV <- distance(data_micro.ps_Start_CONV, method = "bray")
data_noSNP_dist_Start_CONV <- distance(data_noSNP_Start_CONV, method = "bray")

# Perform PCoA ordination using the vegan::capscale function
data_micro_pcoa_Start_CONV <- capscale(data_micro_dist_Start_CONV ~ 1)
data_noSNP_pcoa_Start_CONV <- capscale(data_noSNP_dist_Start_CONV ~ 1)

# Extract ordination scores (site scores) from both ordination objects
data_micro_scores_Start_CONV <- vegan::scores(data_micro_pcoa_Start_CONV, display = "sites")
data_noSNP_scores_Start_CONV <- vegan::scores(data_noSNP_pcoa_Start_CONV, display = "sites")

# Run Procrustes analysis
procrustes_result_Start_CONV <- protest(data_micro_scores_Start_CONV, data_noSNP_scores_Start_CONV)

# Print Procrustes result summary
print(procrustes_result_Start_CONV)

#Procrustes Sum of Squares (m12 squared):        0.9365 
#Correlation in a symmetric Procrustes rotation: 0.2519 
#Significance:  0.442 

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Start_CONV)

#### End_CONV ####


# Subset the phyloseq objects to keep only samples with "End_CONV" in the "Combined_order" column
data_micro.ps_End_CONV <- subset_samples(data_micro.ps_subset, Combined_order == "End_CONV")
data_noSNP_End_CONV <- subset_samples(data_noSNP_subset, Combined_order == "End_CONV")

# Find the shared sample names between both phyloseq objects
shared_samples_End_CONV <- intersect(sample_names(data_micro.ps_End_CONV), sample_names(data_noSNP_End_CONV))

# Subset both phyloseq objects to include only the shared samples
data_micro.ps_End_CONV <- prune_samples(sample_names(data_micro.ps_End_CONV) %in% shared_samples_End_CONV, data_micro.ps_End_CONV)
data_noSNP_End_CONV <- prune_samples(sample_names(data_noSNP_End_CONV) %in% shared_samples_End_CONV, data_noSNP_End_CONV)

# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
data_micro_dist_End_CONV <- distance(data_micro.ps_End_CONV, method = "bray")
data_noSNP_dist_End_CONV <- distance(data_noSNP_End_CONV, method = "bray")

# Perform PCoA ordination using the vegan::capscale function
data_micro_pcoa_End_CONV <- capscale(data_micro_dist_End_CONV ~ 1)
data_noSNP_pcoa_End_CONV <- capscale(data_noSNP_dist_End_CONV ~ 1)

# Extract ordination scores (site scores) from both ordination objects
data_micro_scores_End_CONV <- vegan::scores(data_micro_pcoa_End_CONV, display = "sites")
data_noSNP_scores_End_CONV <- vegan::scores(data_noSNP_pcoa_End_CONV, display = "sites")

# Run Procrustes analysis
procrustes_result_End_CONV <- protest(data_micro_scores_End_CONV, data_noSNP_scores_End_CONV)

# Print Procrustes result summary
print(procrustes_result_End_CONV)

#Procrustes Sum of Squares (m12 squared):        0.8265 
#Correlation in a symmetric Procrustes rotation: 0.4166 
#Significance:  0.079 

# Optional: Plot the Procrustes analysis
plot(procrustes_result_End_CONV)


#### End_RWA ####


# Subset the phyloseq objects to keep only samples with "End_RWA" in the "Combined_order" column
data_micro.ps_End_RWA <- subset_samples(data_micro.ps_subset, Combined_order == "End_RWA")
data_noSNP_End_RWA <- subset_samples(data_noSNP_subset, Combined_order == "End_RWA")

# Find the shared sample names between both phyloseq objects
shared_samples_End_RWA <- intersect(sample_names(data_micro.ps_End_RWA), sample_names(data_noSNP_End_RWA))

# Subset both phyloseq objects to include only the shared samples
data_micro.ps_End_RWA <- prune_samples(sample_names(data_micro.ps_End_RWA) %in% shared_samples_End_RWA, data_micro.ps_End_RWA)
data_noSNP_End_RWA <- prune_samples(sample_names(data_noSNP_End_RWA) %in% shared_samples_End_RWA, data_noSNP_End_RWA)

# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
data_micro_dist_End_RWA <- distance(data_micro.ps_End_RWA, method = "bray")
data_noSNP_dist_End_RWA <- distance(data_noSNP_End_RWA, method = "bray")

# Perform PCoA ordination using the vegan::capscale function
data_micro_pcoa_End_RWA <- capscale(data_micro_dist_End_RWA ~ 1)
data_noSNP_pcoa_End_RWA <- capscale(data_noSNP_dist_End_RWA ~ 1)

# Extract ordination scores (site scores) from both ordination objects
data_micro_scores_End_RWA <- vegan::scores(data_micro_pcoa_End_RWA, display = "sites")
data_noSNP_scores_End_RWA <- vegan::scores(data_noSNP_pcoa_End_RWA, display = "sites")

# Run Procrustes analysis
procrustes_result_End_RWA <- protest(data_micro_scores_End_RWA, data_noSNP_scores_End_RWA)

# Print Procrustes result summary
print(procrustes_result_End_RWA)

# Procrustes Sum of Squares (m12 squared):        0.7928 
# Correlation in a symmetric Procrustes rotation: 0.4552 
# Significance:  0.029 

# Optional: Plot the Procrustes analysis
plot(procrustes_result_End_RWA)
