# Summarize metadata counts
micro.summarized_counts <- sample_data(data_micro.ps) %>%
  group_by(sample_type) %>% 
  summarize(total_16S_reads = sum(raw_paired_reads), mean_16S_reads = mean(raw_paired_reads),
            min_16S_reads = min(raw_paired_reads),max_16S_reads = max(raw_paired_reads),
            median_16S_reads = median(raw_paired_reads),total_ASV_counts = sum(ASV_hits), mean_ASV_counts = mean(ASV_hits),
            min_ASV_counts = min(ASV_hits),max_ASV_counts = max(ASV_hits),
            median_ASV_counts = median(ASV_hits))

micro.summarized_counts
#write.csv(summarized_counts,"Final_R_analysis/microbiome_summarized_counts.csv")

sample_data(data_micro.ps) %>%
  dplyr::group_by(Combined_order) %>% 
  dplyr::summarize(total_16S_reads = sum(raw_paired_reads), mean_16S_reads = mean(raw_paired_reads),
            min_16S_reads = min(raw_paired_reads),max_16S_reads = max(raw_paired_reads),
            median_16S_reads = median(raw_paired_reads),total_ASV_counts = sum(ASV_hits), mean_ASV_counts = mean(ASV_hits),
            min_ASV_counts = min(ASV_hits),max_ASV_counts = max(ASV_hits),
            median_ASV_counts = median(ASV_hits))


# Script to summarize counts
### split NTC and EBs from samples
controls <- subset_samples(data_micro.ps, sample_type=="Extraction Blank" | sample_type=="Negative Template Control" | sample_type=="Rinsate Blank")
controls <- prune_taxa(taxa_sums(controls) > 0, controls)
controls # 168 taxa, 21 control samples
sort(sample_sums(controls)) # 
controls_family <- tax_glom(controls, taxrank = "Family", NArm = F)
controls_family_melt <- psmelt(controls_family)

ggplot(controls_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  labs(y= "Abundance (Number of ASVs per sample)") +
  geom_bar(stat = "summary", colour = "black") +
  #scale_fill_manual() +
  theme(axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 16, angle = 45, hjust = 0.95, vjust = 0.95))

### samples  ####

## # lets look at the number of reads per sample and the distribution
sample_sum_df <- data.frame(sum = sample_sums(samples))
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 5000) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) # looks good

# some QC checks
min(sample_sums(samples)) # 1
max(sample_sums(samples)) # 521,743
mean(sample_sums(samples)) # 66,091
median(sample_sums(samples)) # 39,762
sort(sample_sums(samples)) # 
# HFT-F-114D_S72_L001  631       HFT-F-115C_S78_L001 2163     HFT-F-115A_S76_L001 3579
                                                       
