# Core SNV analysis
library(microbiome)

# check only keep new gene groups that show up at the End of the trial (compared to the start)
# compare in the conventional beef and look for similarities


#core.taxa.standard <- core_members(SNV_data_noSNP, detection = 1, prevalence = 50/100)
#length(core.taxa.standard)

#
##
### Fecal samples ####
##
#

#### Conventional fecal samples ####

# Start CONV
Start_CONV_core <- subset_samples(SNV_data_noSNP, Combined_order=="Start_CONV")
Start_CONV_core.taxa.standard <- core_members(Start_CONV_core, detection = 1, prevalence = 50/100)
length(Start_CONV_core.taxa.standard) # 23688

# End CONV
End_CONV_core <- subset_samples(SNV_data_noSNP, Combined_order=="End_CONV")
End_CONV_core.taxa.standard <- core_members(End_CONV_core, detection = 1, prevalence = 50/100)
length(End_CONV_core.taxa.standard) # 24045
str(End_CONV_core.taxa.standard)

unique_to_End_CONV <- setdiff(End_CONV_core.taxa.standard, Start_CONV_core.taxa.standard)
print(unique_to_End_CONV)
length(unique_to_End_CONV) # 973

#### RWA fecal samples ####

# Start RWA
Start_RWA_core <- subset_samples(SNV_data_noSNP, Combined_order=="Start_RWA")
Start_RWA_core.taxa.standard <- core_members(Start_RWA_core, detection = 1, prevalence = 50/100)
length(Start_RWA_core.taxa.standard) # 23379

# End RWA
End_RWA_core <- subset_samples(SNV_data_noSNP, Combined_order=="End_RWA")
End_RWA_core.taxa.standard <- core_members(End_RWA_core, detection = 1, prevalence = 50/100)
length(End_RWA_core.taxa.standard) # 23209
str(End_RWA_core.taxa.standard)

unique_to_End_RWA <- setdiff(End_RWA_core.taxa.standard, Start_RWA_core.taxa.standard)
length(unique_to_End_RWA) # 793


### Meat samples ####

#### Conventional beef samples #####
beef_CONV_core <- subset_samples(SNV_data_noSNP, Combined_order=="Conventional Beef")
beef_CONV_core.taxa.standard <- core_members(beef_CONV_core, detection = 1, prevalence = 50/100)
length(beef_CONV_core.taxa.standard) # 403
str(beef_CONV_core.taxa.standard)

#### RWA beef samples #####
beef_RWA_core <- subset_samples(SNV_data_noSNP, Combined_order=="RWA Beef")
beef_RWA_core.taxa.standard <- core_members(beef_RWA_core, detection = 1, prevalence = 50/100)
length(beef_RWA_core.taxa.standard) # 1018
str(beef_RWA_core.taxa.standard)




### Set comparisons ####

#### CONV comparisons #####
# Compare unique to end with core CONV SNVs
beef_and_uniqueEnd_CONV <- intersect(unique_to_End_CONV, beef_CONV_core.taxa.standard)
length(beef_and_uniqueEnd_CONV) # 0

# Compare to start with core CONV SNVs
beef_and_Start_CONV <- intersect(Start_CONV_core.taxa.standard, beef_CONV_core.taxa.standard)
length(beef_and_Start_CONV) # 260

# Compare to end with core CONV SNVs
beef_and_End_CONV <- intersect(End_CONV_core.taxa.standard, beef_CONV_core.taxa.standard)
length(beef_and_End_CONV) # 249

#### RWA comparisons #####
# Compare unique to end with core RWA SNVs
beef_and_uniqueEnd_RWA <- intersect(unique_to_End_RWA, beef_RWA_core.taxa.standard)
length(beef_and_uniqueEnd_RWA) # 3

# Compare to start with core RWA SNVs
beef_and_Start_RWA <- intersect(Start_RWA_core.taxa.standard, beef_RWA_core.taxa.standard)
length(beef_and_Start_RWA) # 339

# Compare to end with core RWA SNVs
beef_and_End_RWA <- intersect(End_RWA_core.taxa.standard, beef_RWA_core.taxa.standard)
length(beef_and_End_RWA) # 291


