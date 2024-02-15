# Load AMR++ resistome data into phyloseq object####

# Read in AMR count matrix ####
# V3 results with SNP confirmation
amr_data <- read.table('AMR/deduped_SNPconfirmed_AMR_analytic_matrix.csv', header=T, row.names=1, sep=',', quote = "")
# convert this into an 'otu_table' format required for phyloseq
amr_data <- otu_table(amr_data, taxa_are_rows = T)

# Read in gene annotations ####
annotations <- read.table('AMR/megares_annotations_v3.00.csv', header=T, row.names=1, sep=",", quote = "")
#convert this into a 'taxonomy table' format required for phyloseq
annotations <- phyloseq::tax_table(as.matrix(annotations))

# Read in Sample metadata ####
amr_metadata <- read.table('AMR/HFT_AMR_Metadata.csv', header=T, sep=',', row.names = 1)
# convert to 'sample_data' format for phyloseq
amr_metadata <- sample_data(amr_metadata)

# Merge annotations, count matrix, and metadata into a phyloseq object ####
AMR_data.ps <- merge_phyloseq(amr_data, annotations, amr_metadata)


# Example commands to modify phyloseq object and metadata ####

# Change the label in the metadata for a certain column to something else
sample_data(AMR_data.ps)$sample_type[sample_data(AMR_data.ps)$sample_type == "Beef sample"] <- "Meat sample"

# Example of making new column, then filling it out based on other conditions in the metadata
## Combined order
sample_data(AMR_data.ps)$Combined_order
sample_data(AMR_data.ps)$Combined_order[sample_data(AMR_data.ps)$treatment_order == "AB" & sample_data(AMR_data.ps)$pool_timepoint == "A"] <- "Start_RWA"
sample_data(AMR_data.ps)$Combined_order[sample_data(AMR_data.ps)$treatment_order == "AB" & sample_data(AMR_data.ps)$pool_timepoint == "B"] <- "Midpoint_RWA"
sample_data(AMR_data.ps)$Combined_order[sample_data(AMR_data.ps)$treatment_order == "AB" & sample_data(AMR_data.ps)$pool_timepoint == "C"] <- "End_RWA"
sample_data(AMR_data.ps)$Combined_order[sample_data(AMR_data.ps)$treatment_order == "AB" & sample_data(AMR_data.ps)$pool_timepoint == "D"] <- "Start_CONV"
sample_data(AMR_data.ps)$Combined_order[sample_data(AMR_data.ps)$treatment_order == "AB" & sample_data(AMR_data.ps)$pool_timepoint == "E"] <- "Midpoint_CONV"
sample_data(AMR_data.ps)$Combined_order[sample_data(AMR_data.ps)$treatment_order == "AB" & sample_data(AMR_data.ps)$pool_timepoint == "F"] <- "End_CONV"
sample_data(AMR_data.ps)$Combined_order[sample_data(AMR_data.ps)$treatment_order == "AB" & sample_data(AMR_data.ps)$pool_timepoint == "G"] <- "Final_washout"

sample_data(AMR_data.ps)$Combined_order[sample_data(AMR_data.ps)$treatment_order == "BA" & sample_data(AMR_data.ps)$pool_timepoint == "A"] <- "Start_CONV"
sample_data(AMR_data.ps)$Combined_order[sample_data(AMR_data.ps)$treatment_order == "BA" & sample_data(AMR_data.ps)$pool_timepoint == "B"] <- "Midpoint_CONV"
sample_data(AMR_data.ps)$Combined_order[sample_data(AMR_data.ps)$treatment_order == "BA" & sample_data(AMR_data.ps)$pool_timepoint == "C"] <- "End_CONV"
sample_data(AMR_data.ps)$Combined_order[sample_data(AMR_data.ps)$treatment_order == "BA" & sample_data(AMR_data.ps)$pool_timepoint == "D"] <- "Start_RWA"
sample_data(AMR_data.ps)$Combined_order[sample_data(AMR_data.ps)$treatment_order == "BA" & sample_data(AMR_data.ps)$pool_timepoint == "E"] <- "Midpoint_RWA"
sample_data(AMR_data.ps)$Combined_order[sample_data(AMR_data.ps)$treatment_order == "BA" & sample_data(AMR_data.ps)$pool_timepoint == "F"] <- "End_RWA"
sample_data(AMR_data.ps)$Combined_order[sample_data(AMR_data.ps)$treatment_order == "BA" & sample_data(AMR_data.ps)$pool_timepoint == "G"] <- "Final_washout"

sample_data(AMR_data.ps)$Combined_order[sample_data(AMR_data.ps)$treatment == "Conventional" & sample_data(AMR_data.ps)$sample_type == "Beef sample"] <- "Conventional Beef"
sample_data(AMR_data.ps)$Combined_order[sample_data(AMR_data.ps)$treatment == "RWA" & sample_data(AMR_data.ps)$sample_type == "Beef sample"] <- "RWA Beef"


## Longitudinal analysis
sample_data(AMR_data.ps)$Week
sample_data(AMR_data.ps)$Week[sample_data(AMR_data.ps)$Combined_order == "Start_RWA"] <- "1"
sample_data(AMR_data.ps)$Week[sample_data(AMR_data.ps)$Combined_order == "Midpoint_RWA"] <- "2"
sample_data(AMR_data.ps)$Week[sample_data(AMR_data.ps)$Combined_order == "End_RWA"] <- "3"
sample_data(AMR_data.ps)$Week[sample_data(AMR_data.ps)$Combined_order == "Start_CONV"] <- "1"
sample_data(AMR_data.ps)$Week[sample_data(AMR_data.ps)$Combined_order == "Midpoint_CONV"] <- "2"
sample_data(AMR_data.ps)$Week[sample_data(AMR_data.ps)$Combined_order == "End_CONV"] <- "3"
sample_data(AMR_data.ps)$Week[sample_data(AMR_data.ps)$Combined_order == "Final_washout"] <- "4"

sample_data(AMR_data.ps)$Week <- as.numeric(sample_data(AMR_data.ps)$Week)


# Can write out the metadata file
#write.csv(as.matrix(sample_data(AMR_data.ps)), file = "metadata.csv", row.names = TRUE)


# Remove Rinsate Pools
AMR_data.ps <- subset_samples(AMR_data.ps, sample_type != "Rinsate Pool")
AMR_data.ps # 216 samples


# #### splitting of the genes requiring SNP confirmation
SNPconfirm <- subset_taxa(AMR_data.ps, snp=="RequiresSNPConfirmation")
SNPconfirm # 319 of the 2912 genes require SNP confirmation

data_noSNP <- subset_taxa(AMR_data.ps, snp!="RequiresSNPConfirmation")
data_noSNP # 2593 genes remain


# Check sample sums
sort(sample_sums(AMR_data.ps))




# Make metadata file to use later
AMR_mapfile <- sample_data(AMR_data.ps)

