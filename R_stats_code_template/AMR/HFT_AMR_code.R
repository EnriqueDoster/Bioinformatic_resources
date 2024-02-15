## R analysis

# import data
AMRPPdata <- import_biom("SamDedup_AMR_analytic_matrix_phyloseq.biom")
#write.csv(sample_names(AMRPPdata),"sample-names.csv") # to build out my metadata

## import metadata
map_file <- import_qiime_sample_data("AMR_metadata.txt")
map_file

data <- merge_phyloseq(AMRPPdata, map_file)
data # 225 samples, 2905 taxa

# check the names of our ranks
rank_names(data) # "Rank1" - "Rank7" not ideal, lets change em
colnames(tax_table(data)) <- c("Type","Class","Mechanism","Group","Gene","SNP")
rank_names(data) # beauty, now they are named properly

### need to split up the genes needing SNP confirmation:
data_SNPconfirm <- subset_taxa(data, SNP=="yes", T) # 287 of the 2905 genes
data_noSNP <- subset_taxa(data, SNP=="no") 
data_noSNP # 2618 of the 1775 genes

## removing sequins
data_noSNP <- subset_taxa(data_noSNP, Class != "sequins")
sum(taxa_sums(data_noSNP)==0) # great
data_noSNP #2,552 taxa

## Extracting sequins
data_sequins <- subset_taxa(data, Class == "sequins")
sum(taxa_sums(data_sequins)==0) # great
data_sequins #66 taxa


data_melted <- psmelt(data_sequins)

data_melted <- data_melted %>%
  group_by(Sample) %>%
  dplyr::summarize(sum_sequins = sum(Abundance))

data_melted_nosequins <- psmelt(data_noSNP)
data_melted_nosequins <- data_melted_nosequins %>%
  group_by(Sample) %>%
  dplyr::summarize(sum_bacti = sum(Abundance))

merged_counts <- cbind(data_melted_nosequins,data_melted)
merged_counts <- merged_counts[,-3]
sequins_summary <- merged_counts %>%
  group_by(Sample) %>%
  dplyr::summarize(percent_sequins = sum_sequins / (sum_sequins + sum_bacti) * 100)
  
boxplot(sequins_summary$percent_sequins)

ggplot(sequins_summary, aes(x= Sample, y= percent_sequins)) + theme_bw() +
  geom_boxplot()

# some QC checks
min(sequins_summary$percent_sequins) # 25 samples with 0 sequins
max(sequins_summary$percent_sequins) # 38.11711
mean(sequins_summary$percent_sequins) # 3.190753
median(sequins_summary$percent_sequins) # 1.933205
sort(sequins_summary$percent_sequins)




# some QC checks
min(sample_sums(data_noSNP)) # 539
max(sample_sums(data_noSNP)) # 885,671
sort(sample_sums(data_noSNP)) # 539, 552, 3641, 5177, 5294, 6591, 16,532, 16,886...
# for now, going to cut off those below 3921 but I might be inclined to do 16,886

data_noSNP <- prune_samples(sample_sums(data_noSNP) > 3500, data_noSNP)
data_noSNP # lost two samples and now have 2618 taxa
any(taxa_sums(data_noSNP)==0) # no, good