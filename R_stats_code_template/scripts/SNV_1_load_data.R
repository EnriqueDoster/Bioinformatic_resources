# load metaSNV results
# read in metaSNV count matrix for megaresV3 results with sNP confirmation
SNV_data <- read.table('AMR/metaSNV/metaSNV_analytic_matrix.csv', header=T, row.names=1, sep=',', quote = "")
# convert this into an 'otu_table' format required for phyloseq
SNV_data <- otu_table(SNV_data, taxa_are_rows = T)


#read in SNV_annotations
SNV_annotations <- read.table('AMR/metaSNV/metaSNV_annotations.csv', header=T, row.names=1, sep=",", quote = "")
#convert this into a 'taxonomy table' format required for phyloseq
SNV_annotations <- phyloseq::tax_table(as.matrix(SNV_annotations))

# read in TE metadata
SNV_metadata <- read.table('AMR/metaSNV/SNV_Metadata.txt', header=T, sep='\t', row.names = 1, quote = "")
# convert to 'sample_data' format for phyloseq
SNV_metadata <- sample_data(SNV_metadata)

# merge the annotations, the count matrix, and metadata into a phyloseq object
SNV_data.ps <- merge_phyloseq(SNV_data, SNV_annotations, SNV_metadata)

# #### splitting of the genes requiring SNP confirmation
SNV_SNPconfirm <- subset_taxa(SNV_data.ps, snp=="RequiresSNPConfirmation")
SNV_SNPconfirm # 71079 of the 160113 genes require SNP confirmation
SNV_data_noSNP <- subset_taxa(SNV_data.ps, snp!="RequiresSNPConfirmation")
SNV_data_noSNP # 89034 genes remain, these are what we will use from now on

# Check SNP removed samples
SNV_data_noSNP # 218 samples


# Remove Rinsate Pools
SNV_data_noSNP <- subset_samples(SNV_data_noSNP, sample_type != "Rinsate Pool")
SNV_data_noSNP # 216 samples

# Check sample sums
sort(sample_sums(SNV_data_noSNP))

# Modify Meat rinsate name
sample_data(SNV_data_noSNP)$sample_type[sample_data(SNV_data_noSNP)$sample_type == "Meat Rinsate"] <- "Beef sample"

sample_data(SNV_data_noSNP)$sample_type

sample_data(SNV_data_noSNP)$Trial
sample_data(SNV_data_noSNP)$Trial[sample_data(SNV_data_noSNP)$treatment_order == "AB" & sample_data(SNV_data_noSNP)$pool_timepoint == "B"] <- "RWA Fecal"
sample_data(SNV_data_noSNP)$Trial[sample_data(SNV_data_noSNP)$treatment_order == "AB" & sample_data(SNV_data_noSNP)$pool_timepoint == "C"] <- "RWA Fecal"
sample_data(SNV_data_noSNP)$Trial[sample_data(SNV_data_noSNP)$treatment_order == "BA" & sample_data(SNV_data_noSNP)$pool_timepoint == "E"] <- "RWA Fecal"
sample_data(SNV_data_noSNP)$Trial[sample_data(SNV_data_noSNP)$treatment_order == "BA" & sample_data(SNV_data_noSNP)$pool_timepoint == "F"] <- "RWA Fecal"

sample_data(SNV_data_noSNP)$Trial[sample_data(SNV_data_noSNP)$treatment_order == "AB" & sample_data(SNV_data_noSNP)$pool_timepoint == "E"] <- "Conventional Fecal"
sample_data(SNV_data_noSNP)$Trial[sample_data(SNV_data_noSNP)$treatment_order == "AB" & sample_data(SNV_data_noSNP)$pool_timepoint == "F"] <- "Conventional Fecal"
sample_data(SNV_data_noSNP)$Trial[sample_data(SNV_data_noSNP)$treatment_order == "BA" & sample_data(SNV_data_noSNP)$pool_timepoint == "B"] <- "Conventional Fecal"
sample_data(SNV_data_noSNP)$Trial[sample_data(SNV_data_noSNP)$treatment_order == "BA" & sample_data(SNV_data_noSNP)$pool_timepoint == "C"] <- "Conventional Fecal"

sample_data(SNV_data_noSNP)$Trial[sample_data(SNV_data_noSNP)$treatment == "Conventional"] <- "Conventional Beef"
sample_data(SNV_data_noSNP)$Trial[sample_data(SNV_data_noSNP)$treatment == "RWA"] <- "RWA Beef"

## Combined order
sample_data(SNV_data_noSNP)$Combined_order
sample_data(SNV_data_noSNP)$Combined_order[sample_data(SNV_data_noSNP)$treatment_order == "AB" & sample_data(SNV_data_noSNP)$pool_timepoint == "A"] <- "Start_RWA"
sample_data(SNV_data_noSNP)$Combined_order[sample_data(SNV_data_noSNP)$treatment_order == "AB" & sample_data(SNV_data_noSNP)$pool_timepoint == "B"] <- "Midpoint_RWA"
sample_data(SNV_data_noSNP)$Combined_order[sample_data(SNV_data_noSNP)$treatment_order == "AB" & sample_data(SNV_data_noSNP)$pool_timepoint == "C"] <- "End_RWA"
sample_data(SNV_data_noSNP)$Combined_order[sample_data(SNV_data_noSNP)$treatment_order == "AB" & sample_data(SNV_data_noSNP)$pool_timepoint == "D"] <- "Start_CONV"
sample_data(SNV_data_noSNP)$Combined_order[sample_data(SNV_data_noSNP)$treatment_order == "AB" & sample_data(SNV_data_noSNP)$pool_timepoint == "E"] <- "Midpoint_CONV"
sample_data(SNV_data_noSNP)$Combined_order[sample_data(SNV_data_noSNP)$treatment_order == "AB" & sample_data(SNV_data_noSNP)$pool_timepoint == "F"] <- "End_CONV"
sample_data(SNV_data_noSNP)$Combined_order[sample_data(SNV_data_noSNP)$treatment_order == "AB" & sample_data(SNV_data_noSNP)$pool_timepoint == "G"] <- "Final_washout"


sample_data(SNV_data_noSNP)$Combined_order[sample_data(SNV_data_noSNP)$treatment_order == "BA" & sample_data(SNV_data_noSNP)$pool_timepoint == "A"] <- "Start_CONV"
sample_data(SNV_data_noSNP)$Combined_order[sample_data(SNV_data_noSNP)$treatment_order == "BA" & sample_data(SNV_data_noSNP)$pool_timepoint == "B"] <- "Midpoint_CONV"
sample_data(SNV_data_noSNP)$Combined_order[sample_data(SNV_data_noSNP)$treatment_order == "BA" & sample_data(SNV_data_noSNP)$pool_timepoint == "C"] <- "End_CONV"
sample_data(SNV_data_noSNP)$Combined_order[sample_data(SNV_data_noSNP)$treatment_order == "BA" & sample_data(SNV_data_noSNP)$pool_timepoint == "D"] <- "Start_RWA"
sample_data(SNV_data_noSNP)$Combined_order[sample_data(SNV_data_noSNP)$treatment_order == "BA" & sample_data(SNV_data_noSNP)$pool_timepoint == "E"] <- "Midpoint_RWA"
sample_data(SNV_data_noSNP)$Combined_order[sample_data(SNV_data_noSNP)$treatment_order == "BA" & sample_data(SNV_data_noSNP)$pool_timepoint == "F"] <- "End_RWA"
sample_data(SNV_data_noSNP)$Combined_order[sample_data(SNV_data_noSNP)$treatment_order == "BA" & sample_data(SNV_data_noSNP)$pool_timepoint == "G"] <- "Final_washout"

sample_data(SNV_data_noSNP)$Combined_order[sample_data(SNV_data_noSNP)$treatment == "Conventional" & sample_data(SNV_data_noSNP)$sample_type == "Beef sample"] <- "Conventional Beef"
sample_data(SNV_data_noSNP)$Combined_order[sample_data(SNV_data_noSNP)$treatment == "RWA" & sample_data(SNV_data_noSNP)$sample_type == "Beef sample"] <- "RWA Beef"


# Make metadata file to use later
SNV_mapfile <- sample_data(SNV_data_noSNP)

