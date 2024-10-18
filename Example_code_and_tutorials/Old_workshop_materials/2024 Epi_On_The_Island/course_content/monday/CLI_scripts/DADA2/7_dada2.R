library(dada2); library(ggplot2); library(gridExtra)

# load environment
load("dada2.RData")

#------------------------------------------------------------------------
# Perform taxonomic assignment using SILVA dataset
# lastest Reference data should be downloaded and stored
#https://benjjneb.github.io/dada2/training.html
# addSpecies function assigns species-level annotation
#------------------------------------------------------------------------
taxa_silva <- assignTaxonomy(asvtab.nochim, "/home/training/epi_on_the_island2024/DBs/SILVA_138.1/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE) 
taxa_silva <- addSpecies(taxa_silva, "/home/training/epi_on_the_island2024/DBs/SILVA_138.1/silva_species_assignment_v138.1.fa.gz")
#------------------------------------------------------------------------

# Save taxonomic assignments and ASV table
write.csv(taxa_silva, 'taxa.csv')
write.csv(asvtab.nochim, 'asvtab.nochim.csv')

#------------------------------------------------------------------------
# Final files-- for phyloseq object
asvtab.nochim  # ASV count table
taxa_silva     # Taxa table
summary_tab    # Track file
#------------------------------------------------------------------------

# save environment
save.image("dada2.RData")
