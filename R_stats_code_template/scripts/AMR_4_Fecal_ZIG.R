library(DAtest)
library(metagenomeSeq)
library(metagMisc)
library(ggrepel)

# 2552 total taxa
data_noSNP.css <- phyloseq_transform_css(data_noSNP, log = F)
# Calculate the top taxa at the Group level

data_noSNP_group <- tax_glom(data_noSNP, taxrank = "group")
data_noSNP_group.rel_abun <- transform_sample_counts(data_noSNP_group, function(x) {x/sum(x)}*100)

melt_data_noSNP.rel_abun <- psmelt(data_noSNP_group.rel_abun)

Median_rel_abun <- melt_data_noSNP.rel_abun %>% 
  group_by(group) %>%
  dplyr::summarise(median_rel_abun = median(Abundance)) %>%
  arrange(-median_rel_abun)

top_goups <- Median_rel_abun %>%
  filter(median_rel_abun > 0.05) %>%
  select(group)

###
####
##### Differential abundance #####
####
###

### 
#### Start here and replace the phyloseq object you want to analyze. 


# 2552 taxa
fecal_data_noSNP.css <- subset_samples(data_noSNP.css, sample_type == "Fecal")


## agglomerate our CSS normalized count table at the class level
fecal_data_group.css <- tax_glom(fecal_data_noSNP.css, taxrank = "group")
# Prune groups with >0 counts
fecal_data_group.css <- prune_taxa(taxa_sums(fecal_data_group.css) > 0, fecal_data_group.css)
# For some reason, aggregation doesn't change the row names, so updating here with the Group level annotation
taxa_names(fecal_data_group.css) <- as.data.frame(tax_table(fecal_data_group.css))$Group


# Trim taxa
filtered_fecal_data_group.css <- preDA(fecal_data_group.css, min.samples=80, min.reads = 1) # 10% of samples (16/159)
# 203 taxa left


## using the DA.zig to run metagenomeSeq's fitZIG to look for DA AMR Classes
# Confirm factor order
sample_data(filtered_fecal_data_group.css)$Combined_order <- factor(sample_data(filtered_fecal_data_group.css)$Combined_order, levels = c("Start_CONV","Midpoint_CONV","End_CONV",
                                                                                                                        "Start_RWA","Midpoint_RWA","End_RWA","Final_washout"))

#ZIG model for DA classes between pool timepoints
DA_zig_fecal <- DA.zig(filtered_fecal_data_group.css,predictor = "Combined_order",paired = "participant_id", p.adj = "BH", allResults = T)

# Extract the fit and design
fit_DA_zig_fecal = DA_zig_fecal@fit
design_DA_zig_fecal = DA_zig_fecal@fit$design

# Make contrasts for pairwise comparisons
# This will change the (Intercept) column to Intercept, we'll need to make the model match this naming convention
contrasts_DA_zig_fecal_AvsB = makeContrasts(Intercept-predictorMidpoint_CONV, levels=design_DA_zig_fecal)
contrasts_DA_zig_fecal_BvsC = makeContrasts(predictorMidpoint_CONV-predictorEnd_CONV, levels=design_DA_zig_fecal)
contrasts_DA_zig_fecal_DvsE= makeContrasts(predictorStart_RWA-predictorMidpoint_RWA, levels=design_DA_zig_fecal)
contrasts_DA_zig_fecal_EvsF= makeContrasts(predictorMidpoint_RWA-predictorEnd_RWA, levels=design_DA_zig_fecal)
contrasts_DA_zig_fecal_CvsF= makeContrasts(predictorEnd_CONV-predictorEnd_RWA, levels=design_DA_zig_fecal)


# We need to change the names of the coefficients in the original model
colnames(fit_DA_zig_fecal$coefficients)[1] <- "Intercept"
colnames(fit_DA_zig_fecal$stdev.unscaled)[1] <- "Intercept"

# Run Ebayes and output table with results for AvsB
fitcontrasts_DA_zig_fecal_AvsB = contrasts.fit(fit_DA_zig_fecal, contrasts_DA_zig_fecal_AvsB)
resEB_fitcontrasts_DA_zig_fecal_AvsB = eBayes(fitcontrasts_DA_zig_fecal_AvsB )
table_DA_zig_fecal_AvsB <- topTable(resEB_fitcontrasts_DA_zig_fecal_AvsB, coef=1, adjust.method="BH",number = 1000)
#write.csv(table_DA_zig_fecal_AvsB, "ZIGfit_results_fecal_AvsB.csv")

# Run Ebayes and output table with results for BvsC
fitcontrasts_DA_zig_fecal_BvsC = contrasts.fit(fit_DA_zig_fecal, contrasts_DA_zig_fecal_BvsC)
resEB_fitcontrasts_DA_zig_fecal_BvsC = eBayes(fitcontrasts_DA_zig_fecal_BvsC )
table_DA_zig_fecal_BvsC <- topTable(resEB_fitcontrasts_DA_zig_fecal_BvsC, coef=1, adjust.method="BH",number = 1000)
#write.csv(table_DA_zig_fecal_BvsC, "ZIGfit_results_fecal_BvsC.csv")

# Run Ebayes and output table with results for DvsE
fitcontrasts_DA_zig_fecal_DvsE = contrasts.fit(fit_DA_zig_fecal, contrasts_DA_zig_fecal_DvsE)
resEB_fitcontrasts_DA_zig_fecal_DvsE = eBayes(fitcontrasts_DA_zig_fecal_DvsE )
table_DA_zig_fecal_DvsE <- topTable(resEB_fitcontrasts_DA_zig_fecal_DvsE, coef=1, adjust.method="BH",number = 1000)
#write.csv(table_DA_zig_fecal_DvsE, "ZIGfit_results_fecal_DvsE.csv")

# Run Ebayes and output table with results for EvsF
fitcontrasts_DA_zig_fecal_EvsF = contrasts.fit(fit_DA_zig_fecal, contrasts_DA_zig_fecal_EvsF)
resEB_fitcontrasts_DA_zig_fecal_EvsF = eBayes(fitcontrasts_DA_zig_fecal_EvsF )
table_DA_zig_fecal_EvsF <- topTable(resEB_fitcontrasts_DA_zig_fecal_EvsF, coef=1, adjust.method="BH",number = 1000)
#write.csv(table_DA_zig_fecal_EvsF, "ZIGfit_results_fecal_EvsF.csv")

# Run Ebayes and output table with results for CvsF
fitcontrasts_DA_zig_fecal_CvsF = contrasts.fit(fit_DA_zig_fecal, contrasts_DA_zig_fecal_CvsF)
resEB_fitcontrasts_DA_zig_fecal_CvsF = eBayes(fitcontrasts_DA_zig_fecal_CvsF )
table_DA_zig_fecal_CvsF <- topTable(resEB_fitcontrasts_DA_zig_fecal_CvsF, coef=1, adjust.method="BH",number = 1000)
#write.csv(table_DA_zig_fecal_CvsF, "ZIGfit_results_fecal_CvsF.csv")


###
####
##### Volcano plot start CONV to Midpoint CONV #####
####
###


## Volcano plot
## Sort by ordered adj.P.Val
table_DA_zig_fecal_AvsB <- table_DA_zig_fecal_AvsB[order(table_DA_zig_fecal_AvsB$adj.P.Val), ] 

## Create a column to indicate which genes to label
table_DA_zig_fecal_AvsB <- table_DA_zig_fecal_AvsB %>%
  mutate(threshold = case_when(adj.P.Val <= 0.05 ~ 'TRUE',
                               adj.P.Val > 0.05 ~ 'FALSE'))

# Now, choose which genes to label based on relative abundance
table_DA_zig_fecal_AvsB$name <- row.names(table_DA_zig_fecal_AvsB)

table_DA_zig_fecal_AvsB <- table_DA_zig_fecal_AvsB %>%
  mutate(labels = case_when((adj.P.Val <= 0.05 & name %in% top_goups$Group) ~ 'TRUE'))

# Count how many groups had significant changes
table_DA_zig_fecal_AvsB %>% group_by(threshold) %>% dplyr::count(threshold)
# Count how many high abundance groups had significant changes
table_DA_zig_fecal_AvsB %>% group_by(labels) %>% dplyr::count(labels)


# Additionally, we can make the size of each point on the graph correspond with 
# the "average expression" of each feature (Phyla/AMR gene)
radius_phylum_node <- sqrt(table_DA_zig_fecal_AvsB$AveExpr/pi)

# Create plot
ggplot(table_DA_zig_fecal_AvsB) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), colour = threshold, size = radius_phylum_node)) +
  scale_color_manual(values=c("Black","Red")) +
  geom_text_repel(aes(x = logFC, y = -log10(adj.P.Val),size = 15 ,label = ifelse(labels == TRUE, rownames(table_DA_zig_fecal_AvsB),""))) +
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

###
####
##### Volcano plot Midpoint CONV to End CONV #####
####
###


## Volcano plot
## Sort by ordered adj.P.Val
table_DA_zig_fecal_BvsC <- table_DA_zig_fecal_BvsC[order(table_DA_zig_fecal_BvsC$adj.P.Val), ] 

## Create a column to indicate which genes to label
table_DA_zig_fecal_BvsC <- table_DA_zig_fecal_BvsC %>%
  mutate(threshold = case_when(adj.P.Val <= 0.05 ~ 'TRUE',
                               adj.P.Val > 0.05 ~ 'FALSE'))

# Now, choose which genes to label based on relative abundance
table_DA_zig_fecal_BvsC$name <- row.names(table_DA_zig_fecal_BvsC)
table_DA_zig_fecal_BvsC <- table_DA_zig_fecal_BvsC %>%
  mutate(labels = case_when((adj.P.Val <= 0.05 & name %in% top_goups$Group) ~ 'TRUE'))

# Count how many groups had significant changes
table_DA_zig_fecal_BvsC %>% group_by(threshold) %>% dplyr::count(threshold)
# Count how many high abundance groups had significant changes
table_DA_zig_fecal_BvsC %>% group_by(labels) %>% dplyr::count(labels)


# Additionally, we can make the size of each point on the graph correspond with 
# the "average expression" of each feature (Phyla/AMR gene)
radius_phylum_node <- sqrt(table_DA_zig_fecal_BvsC$AveExpr/pi)

# Create plot
ggplot(table_DA_zig_fecal_BvsC) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), colour = threshold, size = radius_phylum_node)) +
  scale_color_manual(values=c("Black","Red")) +
  geom_text_repel(aes(x = logFC, y = -log10(adj.P.Val),size = 15 ,label = ifelse(labels == TRUE, rownames(table_DA_zig_fecal_BvsC),""))) +
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







###
####
##### Volcano plot Start RWA vs RWA midpoint #####
####
###


## Volcano plot
## Sort by ordered adj.P.Val
table_DA_zig_fecal_DvsE <- table_DA_zig_fecal_DvsE[order(table_DA_zig_fecal_DvsE$adj.P.Val), ] 

## Create a column to indicate which genes to label
table_DA_zig_fecal_DvsE <- table_DA_zig_fecal_DvsE %>%
  mutate(threshold = case_when(adj.P.Val <= 0.05 ~ 'TRUE',
                               adj.P.Val > 0.05 ~ 'FALSE'))

# Now, choose which genes to label based on relative abundance
table_DA_zig_fecal_DvsE$name <- row.names(table_DA_zig_fecal_DvsE)
table_DA_zig_fecal_DvsE <- table_DA_zig_fecal_DvsE %>%
  mutate(labels = case_when((adj.P.Val <= 0.05 & name %in% top_goups$Group) ~ 'TRUE'))

# Count how many groups had significant changes
table_DA_zig_fecal_DvsE %>% group_by(threshold) %>% dplyr::count(threshold)
# Count how many high abundance groups had significant changes
table_DA_zig_fecal_DvsE %>% group_by(labels) %>% dplyr::count(labels)


# Additionally, we can make the size of each point on the graph correspond with 
# the "average expression" of each feature (Phyla/AMR gene)
radius_phylum_node <- sqrt(table_DA_zig_fecal_DvsE$AveExpr/pi)

# Create plot
ggplot(table_DA_zig_fecal_DvsE) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), colour = threshold, size = radius_phylum_node)) +
  scale_color_manual(values=c("Black","Red")) +
  geom_text_repel(aes(x = logFC, y = -log10(adj.P.Val),size = 15 ,label = ifelse(labels == TRUE, rownames(table_DA_zig_fecal_DvsE),""))) +
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





###
####
##### Volcano plot Midpoint RWA vs End RWA  #####
####
###


## Volcano plot
## Sort by ordered adj.P.Val
table_DA_zig_fecal_EvsF <- table_DA_zig_fecal_EvsF[order(table_DA_zig_fecal_EvsF$adj.P.Val), ] 

## Create a column to indicate which genes to label
table_DA_zig_fecal_EvsF <- table_DA_zig_fecal_EvsF %>%
  mutate(threshold = case_when(adj.P.Val <= 0.05 ~ 'TRUE',
                               adj.P.Val > 0.05 ~ 'FALSE'))

# Now, choose which genes to label based on relative abundance
table_DA_zig_fecal_EvsF$name <- row.names(table_DA_zig_fecal_EvsF)
table_DA_zig_fecal_EvsF <- table_DA_zig_fecal_EvsF %>%
  mutate(labels = case_when((adj.P.Val <= 0.05 & name %in% top_goups$Group) ~ 'TRUE'))

# Count how many groups had significant changes
table_DA_zig_fecal_EvsF %>% group_by(threshold) %>% dplyr::count(threshold)
# Count how many high abundance groups had significant changes
table_DA_zig_fecal_EvsF %>% group_by(labels) %>% dplyr::count(labels)


# Additionally, we can make the size of each point on the graph correspond with 
# the "average expression" of each feature (Phyla/AMR gene)
radius_phylum_node <- sqrt(table_DA_zig_fecal_EvsF$AveExpr/pi)

# Create plot
ggplot(table_DA_zig_fecal_EvsF) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), colour = threshold, size = radius_phylum_node)) +
  scale_color_manual(values=c("Black","Red")) +
  geom_text_repel(aes(x = logFC, y = -log10(adj.P.Val),size = 15 ,label = ifelse(labels == TRUE, rownames(table_DA_zig_fecal_EvsF),""))) +
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



###
####
##### Volcano plot End CONV vs End RWA  #####
####
###


## Volcano plot
## Sort by ordered adj.P.Val
table_DA_zig_fecal_CvsF <- table_DA_zig_fecal_CvsF[order(table_DA_zig_fecal_CvsF$adj.P.Val), ] 

## Create a column to indicate which genes to label
table_DA_zig_fecal_CvsF <- table_DA_zig_fecal_CvsF %>%
  mutate(threshold = case_when(adj.P.Val <= 0.05 ~ 'TRUE',
                               adj.P.Val > 0.05 ~ 'FALSE'))

# Now, choose which genes to label based on relative abundance
table_DA_zig_fecal_CvsF$name <- row.names(table_DA_zig_fecal_CvsF)
table_DA_zig_fecal_CvsF <- table_DA_zig_fecal_CvsF %>%
  mutate(labels = case_when((adj.P.Val <= 0.05 & name %in% top_goups$Group) ~ 'TRUE'))

# Count how many groups had significant changes
table_DA_zig_fecal_CvsF %>% group_by(threshold) %>% dplyr::count(threshold)
# Count how many high abundance groups had significant changes
table_DA_zig_fecal_CvsF %>% group_by(labels) %>% dplyr::count(labels)


# Additionally, we can make the size of each point on the graph correspond with 
# the "average expression" of each feature (Phyla/AMR gene)
radius_phylum_node <- sqrt(table_DA_zig_fecal_CvsF$AveExpr/pi)

# Create plot
ggplot(table_DA_zig_fecal_CvsF) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), colour = threshold, size = radius_phylum_node)) +
  scale_color_manual(values=c("Black","Red")) +
  geom_text_repel(aes(x = logFC, y = -log10(adj.P.Val),size = 15 ,label = ifelse(labels == TRUE, rownames(table_DA_zig_fecal_CvsF),""))) +
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
