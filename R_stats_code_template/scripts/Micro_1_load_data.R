
# import data

# This is the results before the repreps were sequenced
qiimedata <- import_biom("16S/table-with-taxonomy.biom", "16S/tree.nwk", "16S/dna-sequences.fasta")
map_file <- import_qiime_sample_data("16S/HFT_Metadata_8_17_2022.txt") 

# combining sample data with the rest
data <- merge_phyloseq(qiimedata,map_file)

## DATA EXPLORATION
data # we have 229 samples, 4010 taxa

# changing the SILVA style naming (k__Bacteria, etc.)
tax.data <- data.frame(tax_table(data)) # extract the taxonomy table as a data frame

# Modify taxa
# check the names of our ranks
rank_names(data) # "Rank1" - "Rank10" not ideal, lets change em
drops <- c("Rank6","Rank8") # Remove rank 6 and 8 
tax.data <- tax.data[ , !(names(tax.data) %in% drops)]

colnames(tax.data) <- c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species")
colnames(tax.data)

tax.data.names <- changeSILVAtaxa_w_species(tax.data) # this gets rid of the SILVA format

# now to change the NAs to a better naming scheme
for (i in 1:8){ tax.data.names[,i] <- as.character(tax.data.names[,i])} # converting all columns to characters
tax.data.names[is.na(tax.data.names)] <- "" # replacing the NAs with an empty string

# now filling in the empty slots with the highest assigned taxonomy
for (i in 1:nrow(tax.data.names)){
  if (tax.data.names[i,2] == ""){
    domain <- paste("unclassified ", tax.data.names[i,1], sep = "")
    tax.data.names[i, 2:8] <- domain
  } else if (tax.data.names[i,3] == ""){
    kingdom <- paste("unclassified ", tax.data.names[i,2], sep = "")
    tax.data.names[i, 3:8] <- kingdom
  } else if (tax.data.names[i,4] == ""){
    phylum <- paste("unclassified ", tax.data.names[i,3], sep = "")
    tax.data.names[i, 4:8] <- phylum
  } else if (tax.data.names[i,5] == ""){
    class <- paste("unclassified ", tax.data.names[i,4], sep = "")
    tax.data.names[i, 5:8] <- class
  } else if (tax.data.names[i,6] == ""){
    order <- paste("unclassified ", tax.data.names[i,5], sep = "")
    tax.data.names[i, 6:8] <- order
  } else if (tax.data.names[i,7] == ""){
    tax.data.names$Genus[i] <- paste("unclassified ",tax.data.names$Family[i], sep = "")
  } else if (tax.data.names[i,8] == ""){
    tax.data.names$Species[i] <- paste("unclassified ",tax.data.names$Genus[i], sep = "")
  }
}

head(tax.data.names) # great, no more NAs and no more k__
tax_table(data) <- as.matrix(tax.data.names) # re-insert the taxonomy table into the phyloseq object
tail(tax_table(data), 20) # sweet, lookin good!



##### Clean up sample types ####

# lets split up the controls and samples
controls <- phyloseq::subset_samples(data, treatment=="Blank")
controls # 17 samples, 4010 taxa

# Keep only the samples (non-blanks)
samples <- phyloseq::subset_samples(data, pool_timepoint!="Blank")
any(sample_sums(samples)==0)
sum(taxa_sums(samples)==0) # 509 taxa gone
samples <- prune_taxa(taxa_sums(samples) > 0, samples)
samples # 212 samples, 3501  taxa left

## # lets look at the number of reads per sample and the distribution
sample_sum_df <- data.frame(sum = sample_sums(samples))
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) # pretty good

# some QC checks
min(sample_sums(samples)) # 2
max(sample_sums(samples)) # 126,935
sort(sample_sums(samples)) # all good


#samples <- prune_samples(sample_sums(samples) > 100, samples)
#samples # 175 samples, 3501  taxa left

data_micro.ps <- subset_taxa(samples, Domain!="Eukaryota")
data_micro.ps # 216  samples, 4365   taxa left

# Check sample sums
max(sample_sums(data_micro.ps))
sort(sample_sums(data_micro.ps))

# Modify Meet rinsate name
sample_data(data_micro.ps)$sample_type[sample_data(data_micro.ps)$sample_type == "Meat Rinsate"] <- "Beef sample"

sample_data(data_micro.ps)$sample_type

sample_data(data_micro.ps)$Trial
sample_data(data_micro.ps)$Trial[sample_data(data_micro.ps)$treatment_order == "AB" & sample_data(data_micro.ps)$pool_timepoint == "B"] <- "Group 1"
sample_data(data_micro.ps)$Trial[sample_data(data_micro.ps)$treatment_order == "AB" & sample_data(data_micro.ps)$pool_timepoint == "C"] <- "Group 1"
sample_data(data_micro.ps)$Trial[sample_data(data_micro.ps)$treatment_order == "BA" & sample_data(data_micro.ps)$pool_timepoint == "E"] <- "Group 1"
sample_data(data_micro.ps)$Trial[sample_data(data_micro.ps)$treatment_order == "BA" & sample_data(data_micro.ps)$pool_timepoint == "F"] <- "Group 1"

sample_data(data_micro.ps)$Trial[sample_data(data_micro.ps)$treatment_order == "AB" & sample_data(data_micro.ps)$pool_timepoint == "E"] <- "Group 2"
sample_data(data_micro.ps)$Trial[sample_data(data_micro.ps)$treatment_order == "AB" & sample_data(data_micro.ps)$pool_timepoint == "F"] <- "Group 2"
sample_data(data_micro.ps)$Trial[sample_data(data_micro.ps)$treatment_order == "BA" & sample_data(data_micro.ps)$pool_timepoint == "B"] <- "Group 2"
sample_data(data_micro.ps)$Trial[sample_data(data_micro.ps)$treatment_order == "BA" & sample_data(data_micro.ps)$pool_timepoint == "C"] <- "Group 2"

sample_data(data_micro.ps)$Trial[sample_data(data_micro.ps)$treatment == "Conventional"] <- "Conventional Beef"
sample_data(data_micro.ps)$Trial[sample_data(data_micro.ps)$treatment == "RWA"] <- "RWA Beef"

## Combined order
sample_data(data_micro.ps)$Combined_order
sample_data(data_micro.ps)$Combined_order[sample_data(data_micro.ps)$treatment_order == "AB" & sample_data(data_micro.ps)$pool_timepoint == "A"] <- "Start_RWA"
sample_data(data_micro.ps)$Combined_order[sample_data(data_micro.ps)$treatment_order == "AB" & sample_data(data_micro.ps)$pool_timepoint == "B"] <- "Midpoint_RWA"
sample_data(data_micro.ps)$Combined_order[sample_data(data_micro.ps)$treatment_order == "AB" & sample_data(data_micro.ps)$pool_timepoint == "C"] <- "End_RWA"
sample_data(data_micro.ps)$Combined_order[sample_data(data_micro.ps)$treatment_order == "AB" & sample_data(data_micro.ps)$pool_timepoint == "D"] <- "Start_CONV"
sample_data(data_micro.ps)$Combined_order[sample_data(data_micro.ps)$treatment_order == "AB" & sample_data(data_micro.ps)$pool_timepoint == "E"] <- "Midpoint_CONV"
sample_data(data_micro.ps)$Combined_order[sample_data(data_micro.ps)$treatment_order == "AB" & sample_data(data_micro.ps)$pool_timepoint == "F"] <- "End_CONV"
sample_data(data_micro.ps)$Combined_order[sample_data(data_micro.ps)$treatment_order == "AB" & sample_data(data_micro.ps)$pool_timepoint == "G"] <- "Final_washout"


sample_data(data_micro.ps)$Combined_order[sample_data(data_micro.ps)$treatment_order == "BA" & sample_data(data_micro.ps)$pool_timepoint == "A"] <- "Start_CONV"
sample_data(data_micro.ps)$Combined_order[sample_data(data_micro.ps)$treatment_order == "BA" & sample_data(data_micro.ps)$pool_timepoint == "B"] <- "Midpoint_CONV"
sample_data(data_micro.ps)$Combined_order[sample_data(data_micro.ps)$treatment_order == "BA" & sample_data(data_micro.ps)$pool_timepoint == "C"] <- "End_CONV"
sample_data(data_micro.ps)$Combined_order[sample_data(data_micro.ps)$treatment_order == "BA" & sample_data(data_micro.ps)$pool_timepoint == "D"] <- "Start_RWA"
sample_data(data_micro.ps)$Combined_order[sample_data(data_micro.ps)$treatment_order == "BA" & sample_data(data_micro.ps)$pool_timepoint == "E"] <- "Midpoint_RWA"
sample_data(data_micro.ps)$Combined_order[sample_data(data_micro.ps)$treatment_order == "BA" & sample_data(data_micro.ps)$pool_timepoint == "F"] <- "End_RWA"
sample_data(data_micro.ps)$Combined_order[sample_data(data_micro.ps)$treatment_order == "BA" & sample_data(data_micro.ps)$pool_timepoint == "G"] <- "Final_washout"


