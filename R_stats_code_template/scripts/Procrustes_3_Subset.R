## Procrustes analysis
##  Households

# Load required packages
#Load R packages
library(phyloseq)
library(vegan)



#### Participant 117_113 ####
# Subset the phyloseq objects to keep only samples within Shared_1 household

First_subset_117 <- subset_samples(AMR_data.ps, participant_id %in% c("117"))
Shared_117 <- subset_samples(First_subset_117, pool_timepoint %in% c("A","B","C"))
Shared_113 <- subset_samples(AMR_data.ps, participant_id %in% c("113") & pool_timepoint %in% c("A", "B", "C"))
#Shared_113 <- subset_samples(AMR_data.ps, participant_id %in% c("113"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Shared_117_dist <- distance(Shared_117, method = "bray")
Shared_113_dist <- distance(Shared_113, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Shared_117_pcoa <- capscale(Shared_117_dist ~ 1)
Shared_113_pcoa <- capscale(Shared_113_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Shared_117_scores <- vegan::scores(Shared_117_pcoa, display = "sites")
Shared_113_scores <- vegan::scores(Shared_113_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Shared_117_113 <- protest(Shared_117_scores, Shared_113_scores)

# Print Procrustes result summary
print(procrustes_result_Shared_117_113)

#Procrustes Sum of Squares (m12 squared):       2.22e-16 
#Correlation in a symmetric Procrustes rotation: 1
#Significance:  1
#Number of permutations: 5

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Shared_117_113)




#### Participant 117_115 ####
# Subset the phyloseq objects to keep only samples within Shared_1 household

#First_subset_117 <- subset_samples(AMR_data.ps, participant_id %in% c("117"))
Shared_117 <- subset_samples(AMR_data.ps, participant_id %in% c("117"))
Shared_115 <- subset_samples(AMR_data.ps, participant_id %in% c("115") & pool_timepoint %in% c("A", "B", "C"))
#Shared_115 <- subset_samples(AMR_data.ps, participant_id %in% c("115"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Shared_117_dist <- distance(Shared_117, method = "bray")
Shared_115_dist <- distance(Shared_115, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Shared_117_pcoa <- capscale(Shared_117_dist ~ 1)
Shared_115_pcoa <- capscale(Shared_115_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Shared_117_scores <- vegan::scores(Shared_117_pcoa, display = "sites")
Shared_115_scores <- vegan::scores(Shared_115_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Shared_117_115 <- protest(Shared_117_scores, Shared_115_scores)

# Print Procrustes result summary
print(procrustes_result_Shared_117_115)

#Procrustes Sum of Squares (m12 squared):       0 
#Correlation in a symmetric Procrustes rotation: 1
#Significance:  1
#Number of permutations: 5

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Shared_117_115)



#### Participant 117_108 ####
# Subset the phyloseq objects to keep only samples within Shared_1 household

First_subset_117 <- subset_samples(AMR_data.ps, participant_id %in% c("117"))
Shared_117 <- subset_samples(First_subset_117, pool_timepoint %in% c("A","B","C"))
Shared_108 <- subset_samples(AMR_data.ps, participant_id %in% c("108") & pool_timepoint %in% c("A", "B", "C"))
#Shared_108 <- subset_samples(AMR_data.ps, participant_id %in% c("108"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Shared_117_dist <- distance(Shared_117, method = "bray")
Shared_108_dist <- distance(Shared_108, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Shared_117_pcoa <- capscale(Shared_117_dist ~ 1)
Shared_108_pcoa <- capscale(Shared_108_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Shared_117_scores <- vegan::scores(Shared_117_pcoa, display = "sites")
Shared_108_scores <- vegan::scores(Shared_108_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Shared_117_108 <- protest(Shared_117_scores, Shared_108_scores)

# Print Procrustes result summary
print(procrustes_result_Shared_117_108)

#Procrustes Sum of Squares (m12 squared):       0 
#Correlation in a symmetric Procrustes rotation: 1
#Significance:  1
#Number of permutations: 5

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Shared_117_108)





#### Participant 117_130 ####
# Subset the phyloseq objects to keep only samples within Shared_1 household

First_subset_117 <- subset_samples(AMR_data.ps, participant_id %in% c("117"))
Shared_117 <- subset_samples(First_subset_117, pool_timepoint %in% c("A","B","C"))
Shared_130 <- subset_samples(AMR_data.ps, participant_id %in% c("130") & pool_timepoint %in% c("A", "B", "C"))
#Shared_130 <- subset_samples(AMR_data.ps, participant_id %in% c("130"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Shared_117_dist <- distance(Shared_117, method = "bray")
Shared_130_dist <- distance(Shared_130, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Shared_117_pcoa <- capscale(Shared_117_dist ~ 1)
Shared_130_pcoa <- capscale(Shared_130_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Shared_117_scores <- vegan::scores(Shared_117_pcoa, display = "sites")
Shared_130_scores <- vegan::scores(Shared_130_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Shared_117_130 <- protest(Shared_117_scores, Shared_130_scores)

# Print Procrustes result summary
print(procrustes_result_Shared_117_130)

#Procrustes Sum of Squares (m12 squared):       0
#Correlation in a symmetric Procrustes rotation: 1
#Significance:  1
#Number of permutations: 5

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Shared_117_130)



#### Participant 113_130 ####
# Subset the phyloseq objects to keep only samples within Shared_1 household

First_subset_113 <- subset_samples(AMR_data.ps, participant_id %in% c("113"))
Shared_113 <- subset_samples(First_subset_113, pool_timepoint %in% c("A","B","C"))
Shared_130 <- subset_samples(AMR_data.ps, participant_id %in% c("130") & pool_timepoint %in% c("A", "B", "C"))
#Shared_130 <- subset_samples(AMR_data.ps, participant_id %in% c("130"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Shared_113_dist <- distance(Shared_113, method = "bray")
Shared_130_dist <- distance(Shared_130, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Shared_113_pcoa <- capscale(Shared_113_dist ~ 1)
Shared_130_pcoa <- capscale(Shared_130_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Shared_113_scores <- vegan::scores(Shared_113_pcoa, display = "sites")
Shared_130_scores <- vegan::scores(Shared_130_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Shared_113_130 <- protest(Shared_113_scores, Shared_130_scores)

# Print Procrustes result summary
print(procrustes_result_Shared_113_130)

#Procrustes Sum of Squares (m12 squared):       2.22e-16 
#Correlation in a symmetric Procrustes rotation: 1
#Significance:  1
#Number of permutations: 5

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Shared_113_130)





#### Participant 118_113 ####
# Subset the phyloseq objects to keep only samples within Shared_1 household

First_subset_118 <- subset_samples(AMR_data.ps, participant_id %in% c("118"))
Shared_118 <- subset_samples(First_subset_118, pool_timepoint %in% c("A","B","C"))
Shared_113 <- subset_samples(AMR_data.ps, participant_id %in% c("113") & pool_timepoint %in% c("A", "B", "C"))
#Shared_113 <- subset_samples(AMR_data.ps, participant_id %in% c("113"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Shared_118_dist <- distance(Shared_118, method = "bray")
Shared_113_dist <- distance(Shared_113, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Shared_118_pcoa <- capscale(Shared_118_dist ~ 1)
Shared_113_pcoa <- capscale(Shared_113_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Shared_118_scores <- vegan::scores(Shared_118_pcoa, display = "sites")
Shared_113_scores <- vegan::scores(Shared_113_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Shared_118_113 <- protest(Shared_118_scores, Shared_113_scores)

# Print Procrustes result summary
print(procrustes_result_Shared_118_113)

#Procrustes Sum of Squares (m12 squared):        -4.441e-16  
#Correlation in a symmetric Procrustes rotation: 1
#Significance:  1
#Number of permutations: 5

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Shared_118_113)




#### Participant 118_115 ####
# Subset the phyloseq objects to keep only samples within Shared_1 household

First_subset_118 <- subset_samples(AMR_data.ps, participant_id %in% c("118"))
Shared_118 <- subset_samples(First_subset_118, pool_timepoint %in% c("A","B","C"))
Shared_115 <- subset_samples(AMR_data.ps, participant_id %in% c("115") & pool_timepoint %in% c("A", "B", "C"))
#Shared_115 <- subset_samples(AMR_data.ps, participant_id %in% c("115"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Shared_118_dist <- distance(Shared_118, method = "bray")
Shared_115_dist <- distance(Shared_115, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Shared_118_pcoa <- capscale(Shared_118_dist ~ 1)
Shared_115_pcoa <- capscale(Shared_115_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Shared_118_scores <- vegan::scores(Shared_118_pcoa, display = "sites")
Shared_115_scores <- vegan::scores(Shared_115_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Shared_118_115 <- protest(Shared_118_scores, Shared_115_scores)

# Print Procrustes result summary
print(procrustes_result_Shared_118_115)

#Procrustes Sum of Squares (m12 squared):       0
#Correlation in a symmetric Procrustes rotation: 1
#Significance:  1
#Number of permutations: 5

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Shared_118_115)





#### Participant 118_108 ####
# Subset the phyloseq objects to keep only samples within Shared_1 household

First_subset_118 <- subset_samples(AMR_data.ps, participant_id %in% c("118"))
Shared_118 <- subset_samples(First_subset_118, pool_timepoint %in% c("A","B","C"))
Shared_108 <- subset_samples(AMR_data.ps, participant_id %in% c("108") & pool_timepoint %in% c("A", "B", "C"))
#Shared_108 <- subset_samples(AMR_data.ps, participant_id %in% c("108"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Shared_118_dist <- distance(Shared_118, method = "bray")
Shared_108_dist <- distance(Shared_108, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Shared_118_pcoa <- capscale(Shared_118_dist ~ 1)
Shared_108_pcoa <- capscale(Shared_108_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Shared_118_scores <- vegan::scores(Shared_118_pcoa, display = "sites")
Shared_108_scores <- vegan::scores(Shared_108_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Shared_118_108 <- protest(Shared_118_scores, Shared_108_scores)

# Print Procrustes result summary
print(procrustes_result_Shared_118_108)

#Procrustes Sum of Squares (m12 squared):       0
#Correlation in a symmetric Procrustes rotation: 1
#Significance:  1
#Number of permutations: 5

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Shared_118_108)





#### Participant 118_111 ####
# Subset the phyloseq objects to keep only samples within Shared_1 household

First_subset_118 <- subset_samples(AMR_data.ps, participant_id %in% c("118"))
Shared_118 <- subset_samples(First_subset_118, pool_timepoint %in% c("A","B","C"))
Shared_111 <- subset_samples(AMR_data.ps, participant_id %in% c("111") & pool_timepoint %in% c("A", "B", "C"))
#Shared_111 <- subset_samples(AMR_data.ps, participant_id %in% c("111"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Shared_118_dist <- distance(Shared_118, method = "bray")
Shared_111_dist <- distance(Shared_111, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Shared_118_pcoa <- capscale(Shared_118_dist ~ 1)
Shared_111_pcoa <- capscale(Shared_111_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Shared_118_scores <- vegan::scores(Shared_118_pcoa, display = "sites")
Shared_111_scores <- vegan::scores(Shared_111_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Shared_118_111 <- protest(Shared_118_scores, Shared_111_scores)

# Print Procrustes result summary
print(procrustes_result_Shared_118_111)

#Procrustes Sum of Squares (m12 squared):       4.441e-16 
#Correlation in a symmetric Procrustes rotation: 1
#Significance:  1
#Number of permutations: 5

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Shared_118_111)







#### Participant 118_130 ####
# Subset the phyloseq objects to keep only samples within Shared_1 household

First_subset_118 <- subset_samples(AMR_data.ps, participant_id %in% c("118"))
Shared_118 <- subset_samples(First_subset_118, pool_timepoint %in% c("A","B","C"))
Shared_130 <- subset_samples(AMR_data.ps, participant_id %in% c("130") & pool_timepoint %in% c("A", "B", "C"))
#Shared_130 <- subset_samples(AMR_data.ps, participant_id %in% c("130"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Shared_118_dist <- distance(Shared_118, method = "bray")
Shared_130_dist <- distance(Shared_130, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Shared_118_pcoa <- capscale(Shared_118_dist ~ 1)
Shared_130_pcoa <- capscale(Shared_130_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Shared_118_scores <- vegan::scores(Shared_118_pcoa, display = "sites")
Shared_130_scores <- vegan::scores(Shared_130_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Shared_118_130 <- protest(Shared_118_scores, Shared_130_scores)

# Print Procrustes result summary
print(procrustes_result_Shared_118_130)

#Procrustes Sum of Squares (m12 squared):       0
#Correlation in a symmetric Procrustes rotation: 1
#Significance:  1
#Number of permutations: 5

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Shared_118_130)













