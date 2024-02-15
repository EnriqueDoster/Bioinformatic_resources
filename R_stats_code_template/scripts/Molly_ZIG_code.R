# AMR ZIG model
library(ggrepel)
# Load in RData object
# load("C:/Users/enriq/Dropbox/Projects/HumanFeedingTrial/Final_R_analysis/AMR/HFT_ZIG_AMR-TE.RData")

# Calculate the top taxa at the Group level

data_noSNP.css <- tax_glom(data_noSNP.css, taxrank = "Group")
data_noSNP.rel_abun <- transform_sample_counts(data_noSNP.css, function(x) {x/sum(x)}*100)

melt_data_noSNP.rel_abun <- psmelt(data_noSNP.rel_abun)

Median_rel_abun <- melt_data_noSNP.rel_abun %>% 
  group_by(Group) %>%
  dplyr::summarise(median_rel_abun = median(Abundance)) %>%
  arrange(-median_rel_abun)

top_goups <- Median_rel_abun %>%
  filter(median_rel_abun > 0.05) %>%
  select(Group)

###
####
##### Prep input data #####
####
###


# Main data object, with SNP removed and CSS normalized
data_noSNP.css 

# Subset samples
Arm1_TrtA_data_noSNP.css <- subset_samples(data_noSNP.css, arm == "1" & blinded_treatment == "A")
Arm1_TrtB_data_noSNP.css <- subset_samples(data_noSNP.css, arm == "1" & blinded_treatment == "B")
Arm2_TrtA_data_noSNP.css <- subset_samples(data_noSNP.css, arm == "2" & blinded_treatment == "A")
Arm2_TrtB_data_noSNP.css <- subset_samples(data_noSNP.css, arm == "2" & blinded_treatment == "B")

fecal_data_noSNP.css <- subset_samples(data_noSNP.css, sample_type == "Fecal")


# Trim taxa
fecal_data_noSNP.css <- preDA(fecal_data_noSNP.css, min.samples=3, min.reads = 1) # out of 28 samples


###
####
##### Differential abundance #####
####
###

### 
#### Start here and replace the phyloseq object you want to analyze. 

## using the DA.zig to run metagenomeSeq's fitZIG to look for DA AMR Classes

## agglomerate our CSS normalized count table at the class level
Arm1_TrtA_group.css <- tax_glom(Arm1_TrtA_data_noSNP.css, taxrank = "Group")
# Prune groups with >0 counts
Arm1_TrtA_group.css <- prune_taxa(taxa_sums(Arm1_TrtA_group.css) > 0, Arm1_TrtA_group.css)
# For some reason, aggregation doesn't change the row names, so updating here with the Group level annotation
taxa_names(Arm1_TrtA_group.css) <- as.data.frame(tax_table(Arm1_TrtA_group.css))$Group

#ZIG model for DA classes between pool timepoints
DA_zig_Arm1_TrtA <- DA.zig(Arm1_TrtA_group.css,predictor = "pool_timepoint",paired = "participant_id", p.adj = "BH", allResults = T)
DA_DESEQ_Arm1_TrtA <- DA.ds2(Arm1_TrtA_group.css,predictor = "pool_timepoint",paired = "participant_id", p.adj = "BH", allResults = T)

# Extract the fit and design
fit_DA_zig_Arm1_TrtA = DA_zig_Arm1_TrtA@fit
design_DA_zig_Arm1_TrtA = DA_zig_Arm1_TrtA@fit$design

# Make contrasts for pairwise comparisons
# This will change the (Intercept) column to Intercept, we'll need to make the model match this naming convention
contrasts_DA_zig_Arm1_TrtA_AvsB = makeContrasts(Intercept-predictorB, levels=design_DA_zig_Arm1_TrtA)
contrasts_DA_zig_Arm1_TrtA_BvsC = makeContrasts(predictorB-predictorC, levels=design_DA_zig_Arm1_TrtA)
contrasts_DA_zig_Arm1_TrtA_AvsC= makeContrasts(Intercept-predictorB, levels=design_DA_zig_Arm1_TrtA)

# We need to change the names of the coefficients in the original model
colnames(fit_DA_zig_Arm1_TrtA$coefficients)[1] <- "Intercept"
colnames(fit_DA_zig_Arm1_TrtA$stdev.unscaled)[1] <- "Intercept"

# Run Ebayes and output table with results for AvsB
fitcontrasts_DA_zig_Arm1_TrtA_AvsB = contrasts.fit(fit_DA_zig_Arm1_TrtA, contrasts_DA_zig_Arm1_TrtA_AvsB)
resEB_fitcontrasts_DA_zig_Arm1_TrtA_AvsB = eBayes(fitcontrasts_DA_zig_Arm1_TrtA_AvsB )
table_DA_zig_Arm1_TrtA_AvsB <- topTable(resEB_fitcontrasts_DA_zig_Arm1_TrtA_AvsB, coef=1, adjust.method="BH",number = 1000)
write.csv(table_DA_zig_Arm1_TrtA_AvsB, "ZIGfit_results_Arm1_TrtA_AvsB.csv")

# Run Ebayes and output table with results for BvsC
fitcontrasts_DA_zig_Arm1_TrtA_BvsC = contrasts.fit(fit_DA_zig_Arm1_TrtA, contrasts_DA_zig_Arm1_TrtA_BvsC)
resEB_fitcontrasts_DA_zig_Arm1_TrtA_BvsC = eBayes(fitcontrasts_DA_zig_Arm1_TrtA_BvsC )
table_DA_zig_Arm1_TrtA_BvsC <- topTable(resEB_fitcontrasts_DA_zig_Arm1_TrtA_BvsC, coef=1, adjust.method="BH",number = 1000)
write.csv(table_DA_zig_Arm1_TrtA_BvsC, "ZIGfit_results_Arm1_TrtA_BvsC.csv")

# Run Ebayes and output table with results for AvsC
fitcontrasts_DA_zig_Arm1_TrtA_AvsC = contrasts.fit(fit_DA_zig_Arm1_TrtA, contrasts_DA_zig_Arm1_TrtA_AvsC)
resEB_fitcontrasts_DA_zig_Arm1_TrtA_AvsC = eBayes(fitcontrasts_DA_zig_Arm1_TrtA_AvsC )
table_DA_zig_Arm1_TrtA_AvsC <- topTable(resEB_fitcontrasts_DA_zig_Arm1_TrtA_AvsC, coef=1, adjust.method="BH",number = 1000)
write.csv(table_DA_zig_Arm1_TrtA_AvsC, "ZIGfit_results_Arm1_TrtA_AvsC.csv")


###
####
##### Volcano plot example #####
####
###


## Volcano plot
## Sort by ordered adj.P.Val
table_DA_zig_Arm1_TrtA_AvsB <- table_DA_zig_Arm1_TrtA_AvsB[order(table_DA_zig_Arm1_TrtA_AvsB$adj.P.Val), ] 

## Create a column to indicate which genes to label
table_DA_zig_Arm1_TrtA_AvsB <- table_DA_zig_Arm1_TrtA_AvsB %>%
  mutate(threshold = case_when(adj.P.Val <= 0.05 ~ 'TRUE',
                                adj.P.Val > 0.05 ~ 'FALSE'))

# Now, choose which genes to label based on relative abundance
table_DA_zig_Arm1_TrtA_AvsB$name <- row.names(table_DA_zig_Arm1_TrtA_AvsB)
table_DA_zig_Arm1_TrtA_AvsB <- table_DA_zig_Arm1_TrtA_AvsB %>%
  mutate(labels = case_when((adj.P.Val <= 0.05 & name %in% top_goups$Group) ~ 'TRUE'))

# Count how many groups had significant changes
table_DA_zig_Arm1_TrtA_AvsB %>% group_by(threshold) %>% dplyr::count(threshold)
# Count how many high abundance groups had significant changes
table_DA_zig_Arm1_TrtA_AvsB %>% group_by(labels) %>% dplyr::count(labels)


# Additionally, we can make the size of each point on the graph correspond with 
# the "average expression" of each feature (Phyla/AMR gene)
radius_phylum_node <- sqrt(table_DA_zig_Arm1_TrtA_AvsB$AveExpr/pi)

# Create plot
ggplot(table_DA_zig_Arm1_TrtA_AvsB) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), colour = threshold, size = radius_phylum_node)) +
  scale_color_manual(values=c("Black","Red")) +
  geom_text_repel(aes(x = logFC, y = -log10(adj.P.Val),size = 15 ,label = ifelse(labels == TRUE, rownames(table_DA_zig_Arm1_TrtA_AvsB),""))) +
  geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.5) +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=18),
    strip.text.y=element_text(size=18, angle=0),
    axis.text.x=element_text(size=18, angle=0, hjust=0),
    axis.text.y=element_text(size=18),
    axis.title=element_text(size=18),
    legend.position = "none",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    panel.background = element_rect(fill = "white")
  ) +
  #ggtitle("Rare resistome composition by collection time (only classes < 1%)") +
  xlab("logFC") + 
  ylab("-log10 adjusted p-value") 
#


#ggsave("manuscript_docs/Figures/supp-Pen_level_ZIG_by_Time.jpeg", width = 30, height = 20, units = "cm")


