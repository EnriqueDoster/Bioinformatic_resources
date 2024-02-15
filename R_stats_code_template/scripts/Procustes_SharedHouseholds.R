## Procrustes analysis
##  Households

# Load required packages
#Load R packages
library(phyloseq)
library(vegan)



#### Start with  Households  117_118 ####
# Subset the phyloseq objects to keep only samples within Shared_1 household
Shared_117 <- subset_samples(AMR_data.ps, participant_id %in% c("117"))
Shared_118 <- subset_samples(AMR_data.ps, participant_id %in% c("118"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Shared_117_dist <- distance(Shared_117, method = "bray")
Shared_118_dist <- distance(Shared_118, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Shared_117_pcoa <- capscale(Shared_117_dist ~ 1)
Shared_118_pcoa <- capscale(Shared_118_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Shared_117_scores <- vegan::scores(Shared_117_pcoa, display = "sites")
Shared_118_scores <- vegan::scores(Shared_118_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Shared_117_118 <- protest(Shared_117_scores, Shared_118_scores)

# Print Procrustes result summary
print(procrustes_result_Shared_117_118)

#Procrustes Sum of Squares (m12 squared):       -4.441e-16 
#Correlation in a symmetric Procrustes rotation: 1
#Significance:  1
#Number of permutations: 5

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Shared_117_118)





#### Start with Shared Households  103_105 ####
# Subset the phyloseq objects to keep only samples within Shared_1 household
Shared_103 <- subset_samples(AMR_data.ps, participant_id %in% c("103"))
Shared_105 <- subset_samples(AMR_data.ps, participant_id %in% c("105"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Shared_103_dist <- distance(Shared_103, method = "bray")
Shared_105_dist <- distance(Shared_105, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Shared_103_pcoa <- capscale(Shared_103_dist ~ 1)
Shared_105_pcoa <- capscale(Shared_105_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Shared_103_scores <- vegan::scores(Shared_103_pcoa, display = "sites")
Shared_105_scores <- vegan::scores(Shared_105_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Shared_103_105 <- protest(Shared_103_scores, Shared_105_scores)

# Print Procrustes result summary
print(procrustes_result_Shared_103_105)

#Procrustes Sum of Squares (m12 squared):        0.7991  
#Correlation in a symmetric Procrustes rotation: 0.4483 
#Significance:  0.633 
#Number of permutations: 5039

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Shared_103_105)






#### Start with Control Households  108_130 ####
# Subset the phyloseq objects to keep only samples within Control Group
Control_108 <- subset_samples(AMR_data.ps, participant_id %in% c("108"))
Control_130 <- subset_samples(AMR_data.ps, participant_id %in% c("130"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Control_108_dist <- distance(Control_108, method = "bray")
Control_130_dist <- distance(Control_130, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Control_108_pcoa <- capscale(Control_108_dist ~ 1)
Control_130_pcoa <- capscale(Control_130_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Control_108_scores <- vegan::scores(Control_108_pcoa, display = "sites")
Control_130_scores <- vegan::scores(Control_130_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Control_108_130 <- protest(Control_108_scores, Control_130_scores)

# Print Procrustes result summary
print(procrustes_result_Control_108_130)

#Procrustes Sum of Squares (m12 squared):        0.8441  
#Correlation in a symmetric Procrustes rotation: 0.3948  
#Significance:  0.717 
#Number of permutations: 5039

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Control_108_130)





#### Start with Control Households  113_115 ####
# Subset the phyloseq objects to keep only samples within Control Group
Control_113 <- subset_samples(AMR_data.ps, participant_id %in% c("113"))
Control_115 <- subset_samples(AMR_data.ps, participant_id %in% c("115"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Control_113_dist <- distance(Control_113, method = "bray")
Control_115_dist <- distance(Control_115, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Control_113_pcoa <- capscale(Control_113_dist ~ 1)
Control_115_pcoa <- capscale(Control_115_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Control_113_scores <- vegan::scores(Control_113_pcoa, display = "sites")
Control_115_scores <- vegan::scores(Control_115_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Control_113_115 <- protest(Control_113_scores, Control_115_scores)

# Print Procrustes result summary
print(procrustes_result_Control_113_115)

#Procrustes Sum of Squares (m12 squared):         0.5222  
#Correlation in a symmetric Procrustes rotation: 0.6912  
#Significance:  0.66667 

#Number of permutations: 23

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Control_113_115)






#### Start with Control Households  108_104 ####
# Subset the phyloseq objects to keep only samples within Control Group
Control_108 <- subset_samples(AMR_data.ps, participant_id %in% c("108"))
Control_104 <- subset_samples(AMR_data.ps, participant_id %in% c("104"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Control_108_dist <- distance(Control_108, method = "bray")
Control_104_dist <- distance(Control_104, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Control_108_pcoa <- capscale(Control_108_dist ~ 1)
Control_104_pcoa <- capscale(Control_104_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Control_108_scores <- vegan::scores(Control_108_pcoa, display = "sites")
Control_104_scores <- vegan::scores(Control_104_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Control_108_104 <- protest(Control_108_scores, Control_104_scores)

# Print Procrustes result summary
print(procrustes_result_Control_108_104)

#Procrustes Sum of Squares (m12 squared):         0.7033  
#Correlation in a symmetric Procrustes rotation: 0.5447  
#Significance:  0.362 

#Number of permutations: 5039

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Control_108_104)





#### Start with Control Households  108_111 ####
# Subset the phyloseq objects to keep only samples within Control Group
Control_108 <- subset_samples(AMR_data.ps, participant_id %in% c("108"))
Control_111 <- subset_samples(AMR_data.ps, participant_id %in% c("111"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Control_108_dist <- distance(Control_108, method = "bray")
Control_111_dist <- distance(Control_111, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Control_108_pcoa <- capscale(Control_108_dist ~ 1)
Control_111_pcoa <- capscale(Control_111_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Control_108_scores <- vegan::scores(Control_108_pcoa, display = "sites")
Control_111_scores <- vegan::scores(Control_111_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Control_108_111 <- protest(Control_108_scores, Control_111_scores)

# Print Procrustes result summary
print(procrustes_result_Control_108_111)

#Procrustes Sum of Squares (m12 squared):         0.8987  
#Correlation in a symmetric Procrustes rotation: 0.3182  
#Significance:  0.907 

#Number of permutations: 5039

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Control_108_111)







#### Start with Control Households  104_130 ####
# Subset the phyloseq objects to keep only samples within Control Group
Control_104 <- subset_samples(AMR_data.ps, participant_id %in% c("104"))
Control_130 <- subset_samples(AMR_data.ps, participant_id %in% c("130"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Control_104_dist <- distance(Control_104, method = "bray")
Control_130_dist <- distance(Control_130, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Control_104_pcoa <- capscale(Control_104_dist ~ 1)
Control_130_pcoa <- capscale(Control_130_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Control_104_scores <- vegan::scores(Control_104_pcoa, display = "sites")
Control_130_scores <- vegan::scores(Control_130_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Control_104_130 <- protest(Control_104_scores, Control_130_scores)

# Print Procrustes result summary
print(procrustes_result_Control_104_130)

#Procrustes Sum of Squares (m12 squared):         0.9532   
#Correlation in a symmetric Procrustes rotation: 0.2163  
#Significance:  0.951 

#Number of permutations: 5039

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Control_104_130)







#### Start with Control Households  104_111 ####
# Subset the phyloseq objects to keep only samples within Control Group
Control_104 <- subset_samples(AMR_data.ps, participant_id %in% c("104"))
Control_111 <- subset_samples(AMR_data.ps, participant_id %in% c("111"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Control_104_dist <- distance(Control_104, method = "bray")
Control_111_dist <- distance(Control_111, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Control_104_pcoa <- capscale(Control_104_dist ~ 1)
Control_111_pcoa <- capscale(Control_111_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Control_104_scores <- vegan::scores(Control_104_pcoa, display = "sites")
Control_111_scores <- vegan::scores(Control_111_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Control_104_111 <- protest(Control_104_scores, Control_111_scores)

# Print Procrustes result summary
print(procrustes_result_Control_104_111)

#Procrustes Sum of Squares (m12 squared):        0.758   
#Correlation in a symmetric Procrustes rotation: 0.492  
#Significance:  0.489 

#Number of permutations: 5039

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Control_104_111)







#### Start with Control Households  130_111 ####
# Subset the phyloseq objects to keep only samples within Control Group
Control_130 <- subset_samples(AMR_data.ps, participant_id %in% c("130"))
Control_111 <- subset_samples(AMR_data.ps, participant_id %in% c("111"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Control_130_dist <- distance(Control_130, method = "bray")
Control_111_dist <- distance(Control_111, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Control_130_pcoa <- capscale(Control_130_dist ~ 1)
Control_111_pcoa <- capscale(Control_111_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Control_130_scores <- vegan::scores(Control_130_pcoa, display = "sites")
Control_111_scores <- vegan::scores(Control_111_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Control_130_111 <- protest(Control_130_scores, Control_111_scores)

# Print Procrustes result summary
print(procrustes_result_Control_130_111)

#Procrustes Sum of Squares (m12 squared):        0.7756    
#Correlation in a symmetric Procrustes rotation: 0.4737  
#Significance:  0.515  

#Number of permutations: 5039

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Control_130_111)







#### Start with Control Households  106_108 ####
# Subset the phyloseq objects to keep only samples within Control Group
Control_106 <- subset_samples(AMR_data.ps, participant_id %in% c("106"))
Control_108 <- subset_samples(AMR_data.ps, participant_id %in% c("108"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Control_106_dist <- distance(Control_106, method = "bray")
Control_108_dist <- distance(Control_108, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Control_106_pcoa <- capscale(Control_106_dist ~ 1)
Control_108_pcoa <- capscale(Control_108_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Control_106_scores <- vegan::scores(Control_106_pcoa, display = "sites")
Control_108_scores <- vegan::scores(Control_108_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Control_106_108 <- protest(Control_106_scores, Control_108_scores)

# Print Procrustes result summary
print(procrustes_result_Control_106_108)

#Procrustes Sum of Squares (m12 squared):        0.7131    
#Correlation in a symmetric Procrustes rotation: 0.5357  
#Significance:  0.385  

#Number of permutations: 5039

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Control_106_108)







#### Start with Control Households  106_130 ####
# Subset the phyloseq objects to keep only samples within Control Group
Control_106 <- subset_samples(AMR_data.ps, participant_id %in% c("106"))
Control_130 <- subset_samples(AMR_data.ps, participant_id %in% c("130"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Control_106_dist <- distance(Control_106, method = "bray")
Control_130_dist <- distance(Control_130, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Control_106_pcoa <- capscale(Control_106_dist ~ 1)
Control_130_pcoa <- capscale(Control_130_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Control_106_scores <- vegan::scores(Control_106_pcoa, display = "sites")
Control_130_scores <- vegan::scores(Control_130_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Control_106_130 <- protest(Control_106_scores, Control_130_scores)

# Print Procrustes result summary
print(procrustes_result_Control_106_130)

#Procrustes Sum of Squares (m12 squared):        0.8794   
#Correlation in a symmetric Procrustes rotation: 0.3473  
#Significance:  0.863  

#Number of permutations: 5039

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Control_106_130)








#### Start with Control Households  106_104 ####
# Subset the phyloseq objects to keep only samples within Control Group
Control_106 <- subset_samples(AMR_data.ps, participant_id %in% c("106"))
Control_104 <- subset_samples(AMR_data.ps, participant_id %in% c("104"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Control_106_dist <- distance(Control_106, method = "bray")
Control_104_dist <- distance(Control_104, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Control_106_pcoa <- capscale(Control_106_dist ~ 1)
Control_104_pcoa <- capscale(Control_104_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Control_106_scores <- vegan::scores(Control_106_pcoa, display = "sites")
Control_104_scores <- vegan::scores(Control_104_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Control_106_104 <- protest(Control_106_scores, Control_104_scores)

# Print Procrustes result summary
print(procrustes_result_Control_106_104)

#Procrustes Sum of Squares (m12 squared):        0.6995 
#Correlation in a symmetric Procrustes rotation: 0.5482 
#Significance:  0.362 


#Number of permutations: 5039

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Control_106_104)







#### Start with Control Households  106_111 ####
# Subset the phyloseq objects to keep only samples within Control Group
Control_106 <- subset_samples(AMR_data.ps, participant_id %in% c("106"))
Control_111 <- subset_samples(AMR_data.ps, participant_id %in% c("111"))


# Calculate distance matrices for both datasets using the same distance metric, such as Bray-Curtis
Control_106_dist <- distance(Control_106, method = "bray")
Control_111_dist <- distance(Control_111, method = "bray")


# Perform PCoA ordination using the vegan::capscale function
Control_106_pcoa <- capscale(Control_106_dist ~ 1)
Control_111_pcoa <- capscale(Control_111_dist ~ 1)

# Extract ordination scores (site scores) from both ordination objects
Control_106_scores <- vegan::scores(Control_106_pcoa, display = "sites")
Control_111_scores <- vegan::scores(Control_111_pcoa, display = "sites")

# Run Procrustes analysis
procrustes_result_Control_106_111 <- protest(Control_106_scores, Control_111_scores)

# Print Procrustes result summary
print(procrustes_result_Control_106_111)

#Procrustes Sum of Squares (m12 squared):        0.8377 
#Correlation in a symmetric Procrustes rotation: 0.4029 
#Significance:  0.757  


#Number of permutations: 5039

# Optional: Plot the Procrustes analysis
plot(procrustes_result_Control_106_111)



