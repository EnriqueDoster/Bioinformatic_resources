## Procrustes analysis
##  Households

# Load required packages
#Load R packages
library(phyloseq)
library(vegan)



#### Participant 103_104 ####
# Subset the phyloseq objects to keep only samples within Shared_1 household
Shared_103 <- subset_samples(AMR_data.ps, participant_id %in% c("103"))
Shared_104 <- subset_samples(AMR_data.ps, participant_id %in% c("104"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Shared_103_dist <- distance(Shared_103, method = "bray")
Shared_104_dist <- distance(Shared_104, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Shared_103_pcoa <- capscale(Shared_103_dist ~ 1)
Shared_104_pcoa <- capscale(Shared_104_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Shared_103_scores <- vegan::scores(Shared_103_pcoa, display = "sites")
Shared_104_scores <- vegan::scores(Shared_104_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Shared_103_104 <- protest(Shared_103_scores, Shared_104_scores)

# Print Procrustes result summary
print(procrustes_result_Shared_103_104)

#Procrustes Sum of Squares (m12 squared):       0.6605 
#Correlation in a symmetric Procrustes rotation: 0.5826 
#Significance:  0.296
#Number of permutations: 5039

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Shared_103_104)








#### Participant 103_108 ####
# Subset the phyloseq objects to keep only samples within Shared_1 household
Shared_103 <- subset_samples(AMR_data.ps, participant_id %in% c("103"))
Shared_108 <- subset_samples(AMR_data.ps, participant_id %in% c("108"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Shared_103_dist <- distance(Shared_103, method = "bray")
Shared_108_dist <- distance(Shared_108, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Shared_103_pcoa <- capscale(Shared_103_dist ~ 1)
Shared_108_pcoa <- capscale(Shared_108_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Shared_103_scores <- vegan::scores(Shared_103_pcoa, display = "sites")
Shared_108_scores <- vegan::scores(Shared_108_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Shared_103_108 <- protest(Shared_103_scores, Shared_108_scores)

# Print Procrustes result summary
print(procrustes_result_Shared_103_108)

#Procrustes Sum of Squares (m12 squared):       0.5555 
#Correlation in a symmetric Procrustes rotation: 0.6667 
#Significance:  0.116
#Number of permutations: 5039

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Shared_103_108)





#### Participant 103_111 ####
# Subset the phyloseq objects to keep only samples within Shared_1 household
Shared_103 <- subset_samples(AMR_data.ps, participant_id %in% c("103"))
Shared_111 <- subset_samples(AMR_data.ps, participant_id %in% c("111"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Shared_103_dist <- distance(Shared_103, method = "bray")
Shared_111_dist <- distance(Shared_111, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Shared_103_pcoa <- capscale(Shared_103_dist ~ 1)
Shared_111_pcoa <- capscale(Shared_111_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Shared_103_scores <- vegan::scores(Shared_103_pcoa, display = "sites")
Shared_111_scores <- vegan::scores(Shared_111_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Shared_103_111 <- protest(Shared_103_scores, Shared_111_scores)

# Print Procrustes result summary
print(procrustes_result_Shared_103_111)

#Procrustes Sum of Squares (m12 squared):       0.7933  
#Correlation in a symmetric Procrustes rotation: 0.4546  
#Significance:  0.624
#Number of permutations: 5039

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Shared_103_111)





#### Participant 103_130 ####
# Subset the phyloseq objects to keep only samples within Shared_1 household
Shared_103 <- subset_samples(AMR_data.ps, participant_id %in% c("103"))
Shared_130 <- subset_samples(AMR_data.ps, participant_id %in% c("130"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Shared_103_dist <- distance(Shared_103, method = "bray")
Shared_130_dist <- distance(Shared_130, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Shared_103_pcoa <- capscale(Shared_103_dist ~ 1)
Shared_130_pcoa <- capscale(Shared_130_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Shared_103_scores <- vegan::scores(Shared_103_pcoa, display = "sites")
Shared_130_scores <- vegan::scores(Shared_130_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Shared_103_130 <- protest(Shared_103_scores, Shared_130_scores)

# Print Procrustes result summary
print(procrustes_result_Shared_103_130)

#Procrustes Sum of Squares (m12 squared):       0.6275  
#Correlation in a symmetric Procrustes rotation: 0.6103  
#Significance:  0.224
#Number of permutations: 5039

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Shared_103_130)






#### Participant 105_104 ####
# Subset the phyloseq objects to keep only samples within Shared_1 household
Shared_105 <- subset_samples(AMR_data.ps, participant_id %in% c("105"))
Shared_104 <- subset_samples(AMR_data.ps, participant_id %in% c("104"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Shared_105_dist <- distance(Shared_105, method = "bray")
Shared_104_dist <- distance(Shared_104, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Shared_105_pcoa <- capscale(Shared_105_dist ~ 1)
Shared_104_pcoa <- capscale(Shared_104_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Shared_105_scores <- vegan::scores(Shared_105_pcoa, display = "sites")
Shared_104_scores <- vegan::scores(Shared_104_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Shared_105_104 <- protest(Shared_105_scores, Shared_104_scores)

# Print Procrustes result summary
print(procrustes_result_Shared_105_104)

#Procrustes Sum of Squares (m12 squared):       0.3012  
#Correlation in a symmetric Procrustes rotation: 0.8359  
#Significance:  0.006 
#Number of permutations: 5039

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Shared_105_104)





#### Participant 105_108 ####
# Subset the phyloseq objects to keep only samples within Shared_1 household
Shared_105 <- subset_samples(AMR_data.ps, participant_id %in% c("105"))
Shared_108 <- subset_samples(AMR_data.ps, participant_id %in% c("108"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Shared_105_dist <- distance(Shared_105, method = "bray")
Shared_108_dist <- distance(Shared_108, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Shared_105_pcoa <- capscale(Shared_105_dist ~ 1)
Shared_108_pcoa <- capscale(Shared_108_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Shared_105_scores <- vegan::scores(Shared_105_pcoa, display = "sites")
Shared_108_scores <- vegan::scores(Shared_108_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Shared_105_108 <- protest(Shared_105_scores, Shared_108_scores)

# Print Procrustes result summary
print(procrustes_result_Shared_105_108)

#Procrustes Sum of Squares (m12 squared):       0.8579 
#Correlation in a symmetric Procrustes rotation: 0.377  
#Significance:  0.794 
#Number of permutations: 5039

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Shared_105_108)




#### Participant 105_111 ####
# Subset the phyloseq objects to keep only samples within Shared_1 household
Shared_105 <- subset_samples(AMR_data.ps, participant_id %in% c("105"))
Shared_111 <- subset_samples(AMR_data.ps, participant_id %in% c("111"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Shared_105_dist <- distance(Shared_105, method = "bray")
Shared_111_dist <- distance(Shared_111, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Shared_105_pcoa <- capscale(Shared_105_dist ~ 1)
Shared_111_pcoa <- capscale(Shared_111_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Shared_105_scores <- vegan::scores(Shared_105_pcoa, display = "sites")
Shared_111_scores <- vegan::scores(Shared_111_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Shared_105_111 <- protest(Shared_105_scores, Shared_111_scores)

# Print Procrustes result summary
print(procrustes_result_Shared_105_111)

#Procrustes Sum of Squares (m12 squared):       0.7326 
#Correlation in a symmetric Procrustes rotation: 0.5171  
#Significance:  0.415 
#Number of permutations: 5039

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Shared_105_111)





#### Participant 105_130 ####
# Subset the phyloseq objects to keep only samples within Shared_1 household
Shared_105 <- subset_samples(AMR_data.ps, participant_id %in% c("105"))
Shared_130 <- subset_samples(AMR_data.ps, participant_id %in% c("130"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Shared_105_dist <- distance(Shared_105, method = "bray")
Shared_130_dist <- distance(Shared_130, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Shared_105_pcoa <- capscale(Shared_105_dist ~ 1)
Shared_130_pcoa <- capscale(Shared_130_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Shared_105_scores <- vegan::scores(Shared_105_pcoa, display = "sites")
Shared_130_scores <- vegan::scores(Shared_130_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Shared_105_130 <- protest(Shared_105_scores, Shared_130_scores)

# Print Procrustes result summary
print(procrustes_result_Shared_105_130)

#Procrustes Sum of Squares (m12 squared):       0.9212 
#Correlation in a symmetric Procrustes rotation: 0.2807  
#Significance:  0.941 
#Number of permutations: 5039

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Shared_105_130)

