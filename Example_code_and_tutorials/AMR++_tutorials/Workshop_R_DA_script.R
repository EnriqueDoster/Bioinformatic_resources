### SOURCE PREVIOUS SCRIPTS ####
source("/home/training/epi_on_the_island2024/shotgun/phyloseq/1_LOAD_ALL.R")

library(DAtest)
library(eulerr)
library(cowplot)
library(tidyverse)

##
# Confirm data we want to use ####
##
#
sample_data(microbiome_sample)$collection_day

# Use microbiome_sample, agglomerate to phylum, name it something else
microbiome_DA_phylum.ps <- tax_glom(microbiome_sample, taxrank = "Phylum")
# Remove taxa with 0
microbiome_DA_phylum.ps <- prune_taxa(taxa_sums(microbiome_DA_phylum.ps) > 0, microbiome_DA_phylum.ps)
any(taxa_sums(microbiome_DA_phylum.ps)==0) # Double checkingthere are no empty taxa

# Confirm order of factors in collection_day
sample_data(microbiome_DA_phylum.ps)$collection_day <- factor(sample_data(microbiome_DA_phylum.ps)$collection_day, levels = c("Pre_weaning","At_weaning","Post_weaning"))


#
##
# Run test of all DA models, with multiple permutations and spike-ins to calculate FDR, AUC,etc ####
##
# Increase R to something like ~100 for more stable results, for quick tests low values are OK
test_DA_methods <- testDA(microbiome_DA_phylum.ps, predictor = "collection_day",R = 5, cores = 20)

# Summarize the model fit across the multiple runs
summary(test_DA_methods)

# You can also run all tests and store the results, without calculating how well the test performs
all_test <- allDA(microbiome_DA_phylum.ps, predictor = "collection_day", cores = 20)

# An example of pulling out the results for the ZIG model
all_test$results$abc

# Create Venn diagram of shared results
vennDA(all_test, tests=c("zig","abc"))

# Test the power for a particular test
power.abc <- powerDA(microbiome_DA_phylum.ps, predictor = "collection_day", test = "abc", cores = 20)
plot(power.abc)
summary(power.abc)

power.zig <- powerDA(microbiome_DA_phylum.ps, predictor = "collection_day", test = "zig", cores = 20)
plot(power.zig)
summary(power.zig)

# Arrange the plots side by side with labels
combined_plot <- plot_grid(
  plot(power.abc),
  plot(power.zig),
  labels = c("abc", "zig"),
  ncol = 2
)

combined_plot
