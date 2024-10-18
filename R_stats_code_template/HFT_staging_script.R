# Source script to load packages
source("scripts/R_load_packages.R")

# Microbiome analysis

source("scripts/Micro_1_load_data.R") # Main output is data_micro.ps phyloseq object
source("scripts/Micro_2_count_summary.R")
source("scripts/Micro_3_alpha_diversity.R")

# AMR analysis
# You'll also have to redo this analysis because we were previously using the "data_noSNP"
# object, but now that we have snp confirmed counts we'll use the object "AMR_data.ps".
source("scripts/AMR_1_load_data.R") # The main output from this is the "data_noSNP" phyloseq object
source("scripts/AMR_2_alpha_diversity.R")
source("scripts/AMR_3_rel_abun_and_beta_diversity.R")


# metaSNV analysis
# using megares database reduced by 90% sequence identity
source("scripts/SNV_1_load_data.R")

# This results in the object, "SNV_data_noSNP", which you should use to try and replicate the alpha diversity analysis
# and beta diversity/relative abundance tables like we did with the microbiome and resistome


#Panels
source("scripts/AMR_1_load_data.R")
source("scripts/AMR_2_alpha_diversity.R")
source("scripts/SNV_1_load_data.R")
source("scripts/SNV_2_AlphaDiversity.R")
source("scripts/Micro_1_load_data.R")
source("scripts/Micro_2_count_summary.R")
source("scripts/Micro_3_alpha_diversity.R")