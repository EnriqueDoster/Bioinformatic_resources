library(phyloseq); library(ggplot2); library(btools); library(dplyr);
library(tidyr); library(stringr); library(randomcoloR); library(metagMisc);
library(metagenomeSeq); library(GUniFrac); library(randomForest); library(pairwiseAdonis)
library(knitr); library(readr); library(kableExtra); library(scales); library(vegan);
library(ggdendro)

#setwd
setwd("/Users/ljpinnell/Documents/VERO/HumanFeedingTrial/AMRPlusPlus/")

# source some scripts
source("/Users/ljpinnell/Documents/R_Functions_Scripts/MergeLowAbund.R")

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

# some QC checks
min(sample_sums(data_noSNP)) # 539
max(sample_sums(data_noSNP)) # 885,671
sort(sample_sums(data_noSNP)) # 539, 552, 3641, 5177, 5294, 6591, 16,532, 16,886...
# for now, going to cut off those below 3921 but I might be inclined to do 16,886

data_noSNP <- prune_samples(sample_sums(data_noSNP) > 3500, data_noSNP)
data_noSNP # lost two samples and now have 2618 taxa
any(taxa_sums(data_noSNP)==0) # no, good


#############################################################################################
##############################         ALPHA DIVERSITY         ##############################
#############################################################################################
#############################################################################################

alpha_div <- estimate_richness(data_noSNP, measures = c("Observed","Shannon","Simpson","InvSimpson"))
alpha_div
alpha_div.df <- as(sample_data(data_noSNP), "data.frame")
alpha_div_meta <- cbind(alpha_div, alpha_div.df)

alpha_div_feces <- alpha_div_meta[which(alpha_div_meta$sample_type=="Fecal"),]

ggplot(alpha_div_feces, aes(x= treatment_order, y= Observed, fill = treatment_order, colour = treatment_order)) +
  theme_bw() + facet_wrap(~pool_timepoint, ncol = 7) +
  labs(y= "Observed ARGs") +
  geom_boxplot(alpha=0.4) +
  geom_point() +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 28),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 28),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

ggplot(alpha_div_feces, aes(x= treatment_order, y= Shannon, fill = treatment_order, colour = treatment_order)) +
  theme_bw() + facet_wrap(~pool_timepoint, ncol = 7) +
  labs(y= "Observed ARGs") +
  geom_boxplot(alpha=0.4) +
  geom_point() +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 28),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 28),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
  

#stats
# timepoint1  
alpha_div_t1 <- alpha_div_meta[which(alpha_div_meta$pool_timepoint=="A"),]
alpha_div_t1
kruskal.test(alpha_div_t1$Observed, alpha_div_t1$treatment_order) # NS
kruskal.test(alpha_div_t1$Shannon, alpha_div_t1$treatment_order) #NS

# timepoint2  
alpha_div_t2 <- alpha_div_meta[which(alpha_div_meta$pool_timepoint=="B"),]
alpha_div_t2
kruskal.test(alpha_div_t2$Observed, alpha_div_t2$treatment_order) #NS
kruskal.test(alpha_div_t2$Shannon, alpha_div_t2$treatment_order) #NS

# timepoint3  
alpha_div_t3 <- alpha_div_meta[which(alpha_div_meta$pool_timepoint=="C"),]
alpha_div_t3
kruskal.test(alpha_div_t3$Observed, alpha_div_t3$treatment_order)#NS
kruskal.test(alpha_div_t3$Shannon, alpha_div_t3$treatment_order)#NS

# timepoint4  
alpha_div_t4 <- alpha_div_meta[which(alpha_div_meta$pool_timepoint=="D"),]
alpha_div_t4
kruskal.test(alpha_div_t4$Observed, alpha_div_t4$treatment_order)#NS
kruskal.test(alpha_div_t4$Shannon, alpha_div_t4$treatment_order)#NS

# timepoint5 
alpha_div_t5 <- alpha_div_meta[which(alpha_div_meta$pool_timepoint=="E"),]
alpha_div_t5
kruskal.test(alpha_div_t5$Observed, alpha_div_t5$treatment_order)#NS
kruskal.test(alpha_div_t5$Shannon, alpha_div_t5$treatment_order)#NS

# timepoint6
alpha_div_t6 <- alpha_div_meta[which(alpha_div_meta$pool_timepoint=="F"),]
alpha_div_t6
kruskal.test(alpha_div_t6$Observed, alpha_div_t6$treatment_order)#NS
kruskal.test(alpha_div_t6$Shannon, alpha_div_t6$treatment_order)#NS

# timepoint7
alpha_div_t7 <- alpha_div_meta[which(alpha_div_meta$pool_timepoint=="G"),]
alpha_div_t7
kruskal.test(alpha_div_t7$Observed, alpha_div_t7$treatment_order)#NS
kruskal.test(alpha_div_t7$Shannon, alpha_div_t7$treatment_order)#NS
  

#############################################################################################
##############################         BETA DIVERSITY         ###############################
#############################################################################################
#############################################################################################

data_noSNP.css <- phyloseq_transform_css(data_noSNP, log = F)


#############################################################################################
##############################         RELA ABUNDANCE         ###############################
#############################################################################################
#############################################################################################

rel_abund <- transform_sample_counts(data_noSNP.css, function(x) {x/sum(x)}*100)

# agglomerate
ra_class <- tax_glom(rel_abund, taxrank = "Class")
ra_class # 45 classes

ra_mechanism <- tax_glom(rel_abund, taxrank = "Mechanism")
ra_mechanism # 131 mechanisms

ra_group <- tax_glom(rel_abund, taxrank = "Group")
ra_group # 607 groups

######## SPLITTING SAMPLES TYPES #########

ra_class_feces <- subset_samples(ra_class, sample_type=="Fecal")
ra_class_feces <- prune_taxa(taxa_sums(ra_class_feces) > 0, ra_class_feces)
ra_class_feces # 41 class, 159 samples
ra_class_feces_melt <- psmelt(ra_class_feces)


#ra_class_palette <- distinctColorPalette(41)
#write.csv(ra_class_palette,"class_AMR_palette.csv")
#write.csv(tax_table(ra_class_feces), "class_tax.csv")
#write.csv(tax_table(ra_class_meat),"meat-taxa.csv")

ra_class_meat <- subset_samples(ra_class, sample_type=="Meat Rinsate")
ra_class_meat <- prune_taxa(taxa_sums(ra_class_meat) > 0, ra_class_meat)
ra_class_meat # 20 class, 62 samples
ra_class_meat_melt <- psmelt(ra_class_meat)
class_meat_palette <- c("#ECBFDF",	"#B5B4E1",	"#DCC79B",	"#65B962",	"#A2E13E",	"#D0E5C2",	"#76F12D",	"#90A45B",	"#59EC72",	"#C3EAE6",	"#81EEE7",	"#EBD4C9",	"#DEAF64",	"#6D7A9C",	"#DC394F",	"#DECD53",	"#69EA9D",	"#77D5E5",	"#D27967",	"#BBE97A")

#### HUMAN FECES
feces.css <- subset_samples(data_noSNP.css,sample_type=="Fecal")
sum(taxa_sums(feces.css)==0) # 504 not there
feces.css <- prune_taxa(taxa_sums(feces.css) > 0, feces.css)
feces.dist <- vegdist(t(otu_table(feces.css)), method = "bray")
feces.ord <- vegan::metaMDS(comm = t(otu_table(feces.css)), distance = "bray", try = 10, trymax = 100, autotransform = F)
plot_ordination(feces.css, feces.ord, color = "treatment_order") +
  facet_wrap(~pool_timepoint)

## TIMEPOINT 1
feces_t1.css <- subset_samples(feces.css, pool_timepoint=="A")
feces_t1.css <- prune_taxa(taxa_sums(feces_t1.css) > 0, feces_t1.css)
feces_t1.css.df <- as(sample_data(feces_t1.css),"data.frame")
feces_t1.dist <- vegdist(t(otu_table(feces_t1.css)), method = "bray")
feces_t1.ord <- vegan::metaMDS(comm = t(otu_table(feces_t1.css)), distance = "bray", try = 10, trymax = 100, autotransform = F)
plot_ordination(feces_t1.css, feces_t1.ord, color = "treatment_order") +
  theme_bw() +
  labs(title ="TIMEPOINT A") +
  stat_ellipse(aes(fill= treatment_order), geom="polygon", alpha = 0.25) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", size = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

feces_t1.adonis <- adonis2(feces_t1.dist ~ treatment_order, data = feces_t1.css.df, permutations = 9999)
feces_t1.adonis #NS

feces_t1.hclust <- hclust(feces_t1.dist)
plot(feces_t1.hclust)

feces_t1.dendro <- as.dendrogram(feces_t1.hclust)
feces_t1.dendro.data <- dendro_data(feces_t1.dendro, type = "rectangle")
feces_t1.dendro.data

## add in metadata columns
write.csv(feces_t1.dendro.data[["labels"]][["label"]], "feces_t1_dendro_sample_order.csv")

feces_t1_treatment_col <- data.frame(trt = c("AB",	"BA",	"AB",	"BA",	"BA",	"BA",	"AB",	"BA",	"BA",	"AB",	"BA",	"BA",	"AB",	"BA",	"BA",	"AB",	"BA",	"AB",	"AB",	"BA",	"BA",	"BA",	"AB",	"BA",	"BA"))

feces_t1.dendro.data$labels <- cbind(feces_t1.dendro.data$labels, feces_t1_treatment_col)
feces_t1_treatment_col

ggplot(feces_t1.dendro.data$segments) + theme_bw() +
  labs(y= "Ward's Distance", title = "TIMEPOINT A") + coord_flip() + scale_x_reverse() + scale_y_reverse() +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend), size =.75, lineend = "round", linejoin = "round") +
  geom_point(data = feces_t1.dendro.data$labels, aes(x,y, colour = trt, fill = trt, shape = trt),
             size = 7.2, stroke = 0.75, position = position_nudge(x = -0.02)) +
  #geom_text(data = all.dendro.data$labels, aes(x,y, label = day, colour= day), 
  #         size =3.8, position = position_nudge(y = -0.017)) +
  scale_shape_manual(values =c(15,15)) +
  theme(legend.position = "none",
    panel.border = element_blank(),
    plot.title = element_text(size = 30),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(size = 0.75),
    axis.title.x = element_text(size = 24),
    axis.text.x = element_text(size = 14, colour = "black"))

t1_sample_order <- c("HFTF-103A_S12",	"HFTF-123A_S103",	"HFTF-135A_S159",	"HFTF-120A_S82",	"HFTF-111A_S54",	"HFTF-115A_S72",	"HFTF-104A_S19",	"HFTF-128A_S124",	"HFTF-132A_S145",	"HFTF-102A_S5",	"HFTF-110A_S47",	"HFTF-106A_S33",	"HFTF-108B_S41",	"HFTF-101A_S1",	"HFTF-129A_S131",	"HFTF-126A_S117",	"HFTF-133A_S152",	"HFTF-105A_S26",	"HFTF-117A_S76",	"HFTF-118A_S79",	"HFTF-121A_S89",	"HFTF-130A_S138",	"HFTF-124A_S110",	"HFTF-113A_S61",	"HFTF-114A_S65")

ggplot(ra_class_feces_melt, aes(x=Sample, y= Abundance, fill = Class)) +
  theme_bw() + coord_flip() + scale_x_reverse() + scale_y_reverse() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = ra_class_palette) +
  scale_x_discrete(limits = rev(t1_sample_order)) +
  theme(#legend.position = "none",
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.75),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.75, lineend = "square", colour = "black"),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 14, colour = "black"))


## TIMEPOINT 2
feces_t2.css <- subset_samples(feces.css, pool_timepoint=="B")
feces_t2.css <- prune_taxa(taxa_sums(feces_t2.css) > 0, feces_t2.css)
feces_t2.css.df <- as(sample_data(feces_t2.css),"data.frame")
feces_t2.dist <- vegdist(t(otu_table(feces_t2.css)), method = "bray")
feces_t2.ord <- vegan::metaMDS(comm = t(otu_table(feces_t2.css)), distance = "bray", try = 10, trymax = 100, autotransform = F)
plot_ordination(feces_t2.css, feces_t2.ord, color = "treatment_order") +
  theme_bw() +
  labs(title ="TIMEPOINT B") +
  stat_ellipse(aes(fill= treatment_order), geom="polygon", alpha = 0.25) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", size = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

feces_t2.adonis <- adonis2(feces_t2.dist ~ treatment_order, data = feces_t2.css.df, permutations = 9999)
feces_t2.adonis #NS

feces_t2.hclust <- hclust(feces_t2.dist)
plot(feces_t2.hclust)

feces_t2.dendro <- as.dendrogram(feces_t2.hclust)
feces_t2.dendro.data <- dendro_data(feces_t2.dendro, type = "rectangle")
feces_t2.dendro.data

## add in metadata columns
write.csv(feces_t2.dendro.data[["labels"]][["label"]], "feces_t2_dendro_sample_order.csv")

feces_t2_treatment_col <- data.frame(trt = c("BA",	"BA",	"AB",	"BA",	"BA",	"BA",	"AB",	"BA",	"BA",	"BA",	"AB",	"AB",	"AB",	"BA",	"AB",	"AB",	"BA",	"BA",	"AB",	"BA",	"BA",	"BA",	"AB",	"BA",	"BA"))

feces_t2.dendro.data$labels <- cbind(feces_t2.dendro.data$labels, feces_t2_treatment_col)


ggplot(feces_t2.dendro.data$segments) + theme_bw() +
  labs(y= "Ward's Distance", title = "TIMEPOINT B") + coord_flip() + scale_x_reverse() + scale_y_reverse() +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend), size =.75, lineend = "round", linejoin = "round") +
  geom_point(data = feces_t2.dendro.data$labels, aes(x,y, colour = trt, fill = trt, shape = trt),
             size = 7.2, stroke = 0.75, position = position_nudge(x = -0.02)) +
  #geom_text(data = all.dendro.data$labels, aes(x,y, label = day, colour= day), 
  #         size =3.8, position = position_nudge(y = -0.017)) +
  scale_shape_manual(values =c(15,15)) +
  theme(legend.position = "none",
        panel.border = element_blank(),
        plot.title = element_text(size = 30),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(size = 0.75),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 14, colour = "black"))

t2_sample_order <- c("HFTF-111B_S55",	"HFTF-132B_S146",	"HFTF-103B_S13",	"HFTF-115B_S73",	"HFTF-123B_S104",	"HFTF-118B_S80",	"HFTF-135B_S160",	"HFTF-113B_S62",	"HFTF-114B_S66",	"HFTF-101B_S2",	"HFTF-105B_S27",	"HFTF-108C_S42",	"HFTF-126B_S118",	"HFTF-130B_S139",	"HFTF-117B_S77",	"HFTF-102B_S6",	"HFTF-106B_S34",	"HFTF-110B_S48",	"HFTF-124B_S111",	"HFTF-133B_S153",	"HFTF-129B_S132",	"HFTF-120B_S83",	"HFTF-104B_S20",	"HFTF-121B_S90",	"HFTF-128B_S125")

ggplot(ra_class_feces_melt, aes(x=Sample, y= Abundance, fill = Class)) +
  theme_bw() + coord_flip() + scale_x_reverse() + scale_y_reverse() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = ra_class_palette) +
  scale_x_discrete(limits = rev(t2_sample_order)) +
  theme(legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.75),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(size = 0.75, lineend = "square", colour = "black"),
    axis.title.x = element_text(size = 24),
    axis.text.x = element_text(size = 14, colour = "black"))

## TIMEPOINT 3
feces_t3.css <- subset_samples(feces.css, pool_timepoint=="C")
feces_t3.css <- prune_taxa(taxa_sums(feces_t3.css) > 0, feces_t3.css)
feces_t3.css.df <- as(sample_data(feces_t3.css),"data.frame")
feces_t3.dist <- vegdist(t(otu_table(feces_t3.css)), method = "bray")
feces_t3.ord <- vegan::metaMDS(comm = t(otu_table(feces_t3.css)), distance = "bray", try = 10, trymax = 100, autotransform = F)
plot_ordination(feces_t3.css, feces_t3.ord, color = "treatment_order") +
  theme_bw() +
  labs(title ="TIMEPOINT C") +
  stat_ellipse(aes(fill= treatment_order), geom="polygon", alpha = 0.25) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", size = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

feces_t3.adonis <- adonis2(feces_t3.dist ~ treatment_order, data = feces_t3.css.df, permutations = 9999)
feces_t3.adonis #NS

feces_t3.hclust <- hclust(feces_t3.dist)
plot(feces_t3.hclust)

feces_t3.dendro <- as.dendrogram(feces_t3.hclust)
feces_t3.dendro.data <- dendro_data(feces_t3.dendro, type = "rectangle")
feces_t3.dendro.data

## add in metadata columns
write.csv(feces_t3.dendro.data[["labels"]][["label"]], "feces_t3_dendro_sample_order.csv")

feces_t3_treatment_col <- data.frame(trt = c("BA",	"AB",	"AB",	"AB",	"AB",	"BA",	"AB",	"BA",	"BA",	"BA",	"BA",	"BA",	"AB",	"BA",	"BA",	"BA",	"AB",	"BA",	"AB",	"AB",	"BA",	"BA",	"BA",	"BA",	"BA",	"AB"))

feces_t3.dendro.data$labels <- cbind(feces_t3.dendro.data$labels, feces_t3_treatment_col)
feces_t3.dendro.data$labels

ggplot(feces_t3.dendro.data$segments) + theme_bw() +
  labs(y= "Ward's Distance", title = "TIMEPOINT C") + coord_flip() + scale_x_reverse() + scale_y_reverse() +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend), size =.75, lineend = "round", linejoin = "round") +
  geom_point(data = feces_t3.dendro.data$labels, aes(x,y, colour = trt, fill = trt, shape = trt),
             size = 7.2, stroke = 0.75, position = position_nudge(x = -0.02)) +
  #geom_text(data = all.dendro.data$labels, aes(x,y, label = day, colour= day), 
  #         size =3.8, position = position_nudge(y = -0.017)) +
  scale_shape_manual(values =c(15,15)) +
  theme(legend.position = "none",
        panel.border = element_blank(),
        plot.title = element_text(size = 30),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(size = 0.75),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 14, colour = "black"))

t3_sample_order <- c("HFTF-101C_S3",	"HFTF-102C_S7",	"HFTF-103C_S14",	"HFTF-104C_S21",	"HFTF-105C_S28",	"HFTF-106C_S35",	"HFTF-108D_S43",	"HFTF-110C_S49",	"HFTF-111C_S56",	"HFTF-113C_S63",	"HFTF-114C_S67",	"HFTF-115C_S74",	"HFTF-117C_S78",	"HFTF-118C_S81",	"HFTF-120C_S84",	"HFTF-121C_S91",	"HFTF-122C_S98",	"HFTF-123C_S105",	"HFTF-124C_S112",	"HFTF-126C_S119",	"HFTF-128C_S126",	"HFTF-129C_S133",	"HFTF-130C_S140",	"HFTF-132C_S147",	"HFTF-133C_S154",	"HFTF-135C_S161")

ggplot(ra_class_feces_melt, aes(x=Sample, y= Abundance, fill = Class)) +
  theme_bw() + coord_flip() + scale_x_reverse() + scale_y_reverse() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = ra_class_palette) +
  scale_x_discrete(limits = rev(t3_sample_order)) +
  theme(legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.75),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(size = 0.75, lineend = "square", colour = "black"),
    axis.title.x = element_text(size = 24),
    axis.text.x = element_text(size = 14, colour = "black"))

## TIMEPOINT 4
feces_t4.css <- subset_samples(feces.css, pool_timepoint=="D")
feces_t4.css <- prune_taxa(taxa_sums(feces_t4.css) > 0, feces_t4.css)
feces_t4.css.df <- as(sample_data(feces_t4.css),"data.frame")
feces_t4.dist <- vegdist(t(otu_table(feces_t4.css)), method = "bray")
feces_t4.ord <- vegan::metaMDS(comm = t(otu_table(feces_t4.css)), distance = "bray", try = 10, trymax = 100, autotransform = F)
plot_ordination(feces_t4.css, feces_t4.ord, color = "treatment_order") +
  theme_bw() +
  labs(title ="TIMEPOINT D") +
  stat_ellipse(aes(fill= treatment_order), geom="polygon", alpha = 0.25) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", size = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

feces_t4.adonis <- adonis2(feces_t4.dist ~ treatment_order, data = feces_t4.css.df, permutations = 9999)
feces_t4.adonis #NS

feces_t4.hclust <- hclust(feces_t4.dist)
plot(feces_t4.hclust)

feces_t4.dendro <- as.dendrogram(feces_t4.hclust)
feces_t4.dendro.data <- dendro_data(feces_t4.dendro, type = "rectangle")
feces_t4.dendro.data

## add in metadata columns
write.csv(feces_t4.dendro.data[["labels"]][["label"]], "feces_t4_dendro_sample_order.csv")

feces_t4_treatment_col <- data.frame(trt = c("BA",	"BA",	"BA",	"AB",	"BA",	"AB",	"AB",	"BA",	"BA",	"BA",	"BA",	"AB",	"BA",	"AB",	"AB",	"BA",	"BA",	"BA",	"BA",	"AB",	"BA",	"BA",	"AB"))

feces_t4.dendro.data$labels <- cbind(feces_t4.dendro.data$labels, feces_t4_treatment_col)
feces_t4.dendro.data$labels

ggplot(feces_t4.dendro.data$segments) + theme_bw() +
  labs(y= "Ward's Distance", title = "TIMEPOINT D") + coord_flip() + scale_x_reverse() + scale_y_reverse() +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend), size =.75, lineend = "round", linejoin = "round") +
  geom_point(data = feces_t4.dendro.data$labels, aes(x,y, colour = trt, fill = trt, shape = trt),
             size = 8, stroke = 0.75, position = position_nudge(x = -0.02)) +
  #geom_text(data = all.dendro.data$labels, aes(x,y, label = day, colour= day), 
  #         size =3.8, position = position_nudge(y = -0.017)) +
  scale_shape_manual(values =c(15,15)) +
  theme(legend.position = "none",
        panel.border = element_blank(),
        plot.title = element_text(size = 30),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(size = 0.75),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 14, colour = "black"))

t4_sample_order <- c("HFTF-114D_S68",	"HFTF-132D_S148",	"HFTF-133D_S155",	"HFTF-126D_S120",	"HFTF-123D_S106",	"HFTF-124D_S113",	"HFTF-102D_S8",	"HFTF-106D_S36",	"HFTF-101D_S4",	"HFTF-111D_S57",	"HFTF-115D_S75",	"HFTF-103D_S15",	"HFTF-113D_S64",	"HFTF-104D_S22",	"HFTF-105D_S29",	"HFTF-120D_S85",	"HFTF-129D_S134",	"HFTF-110D_S50",	"HFTF-128D_S127",	"HFTF-108E_S44",	"HFTF-130D_S141",	"HFTF-121D_S92",	"HFTF-122D_S99")

ggplot(ra_class_feces_melt, aes(x=Sample, y= Abundance, fill = Class)) +
  theme_bw() + coord_flip() + scale_x_reverse() + scale_y_reverse() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = ra_class_palette) +
  scale_x_discrete(limits = rev(t4_sample_order)) +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.75),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.75, lineend = "square", colour = "black"),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 14, colour = "black"))

## TIMEPOINT 5
feces_t5.css <- subset_samples(feces.css, pool_timepoint=="E")
feces_t5.css <- prune_taxa(taxa_sums(feces_t5.css) > 0, feces_t5.css)
feces_t5.css.df <- as(sample_data(feces_t5.css),"data.frame")
feces_t5.dist <- vegdist(t(otu_table(feces_t5.css)), method = "bray")
feces_t5.ord <- vegan::metaMDS(comm = t(otu_table(feces_t5.css)), distance = "bray", try = 10, trymax = 100, autotransform = F)
plot_ordination(feces_t5.css, feces_t5.ord, color = "treatment_order") +
  theme_bw() +
  labs(title ="TIMEPOINT E") +
  stat_ellipse(aes(fill= treatment_order), geom="polygon", alpha = 0.25) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", size = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

feces_t5.adonis <- adonis2(feces_t5.dist ~ treatment_order, data = feces_t5.css.df, permutations = 9999)
feces_t5.adonis #NS

feces_t5.hclust <- hclust(feces_t5.dist)
plot(feces_t5.hclust)

feces_t5.dendro <- as.dendrogram(feces_t5.hclust)
feces_t5.dendro.data <- dendro_data(feces_t5.dendro, type = "rectangle")
feces_t5.dendro.data

## add in metadata columns
write.csv(feces_t5.dendro.data[["labels"]][["label"]], "feces_t5_dendro_sample_order.csv")

feces_t5_treatment_col <- data.frame(trt = c("BA",	"BA",	"AB",	"AB",	"BA",	"BA",	"BA",	"AB",	"BA",	"BA",	"BA",	"BA",	"BA",	"AB",	"BA",	"BA",	"AB",	"AB",	"AB",	"AB"))

feces_t5.dendro.data$labels <- cbind(feces_t5.dendro.data$labels, feces_t5_treatment_col)
feces_t5.dendro.data$labels

ggplot(feces_t5.dendro.data$segments) + theme_bw() +
  labs(y= "Ward's Distance", title = "TIMEPOINT E") + coord_flip() + scale_x_reverse() + scale_y_reverse() +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend), size =.75, lineend = "round", linejoin = "round") +
  geom_point(data = feces_t5.dendro.data$labels, aes(x,y, colour = trt, fill = trt, shape = trt),
             size = 8, stroke = 0.75, position = position_nudge(x = -0.02)) +
  #geom_text(data = all.dendro.data$labels, aes(x,y, label = day, colour= day), 
  #         size =3.8, position = position_nudge(y = -0.017)) +
  scale_shape_manual(values =c(15,15)) +
  theme(legend.position = "none",
        panel.border = element_blank(),
        plot.title = element_text(size = 30),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(size = 0.75),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 14, colour = "black"))

t5_sample_order <- c("HFTF-114E_S69",	"HFTF-120E_S86",	"HFTF-122E_S100",	"HFTF-104E_S23",	"HFTF-133E_S156",	"HFTF-111E_S58",	"HFTF-123E_S107",	"HFTF-102E_S9",	"HFTF-130E_S142",	"HFTF-129E_S135",	"HFTF-110E_S51",	"HFTF-128E_S128",	"HFTF-132E_S149",	"HFTF-105E_S30",	"HFTF-106E_S37",	"HFTF-121E_S93",	"HFTF-126E_S121",	"HFTF-108F_S45",	"HFTF-103E_S16",	"HFTF-124E_S114")

ggplot(ra_class_feces_melt, aes(x=Sample, y= Abundance, fill = Class)) +
  theme_bw() + coord_flip() + scale_x_reverse() + scale_y_reverse() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = ra_class_palette) +
  scale_x_discrete(limits = rev(t5_sample_order)) +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.75),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.75, lineend = "square", colour = "black"),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 14, colour = "black"))



## TIMEPOINT 6
feces_t6.css <- subset_samples(feces.css, pool_timepoint=="F")
feces_t6.css <- prune_taxa(taxa_sums(feces_t6.css) > 0, feces_t6.css)
feces_t6.css.df <- as(sample_data(feces_t6.css),"data.frame")
feces_t6.dist <- vegdist(t(otu_table(feces_t6.css)), method = "bray")
feces_t6.ord <- vegan::metaMDS(comm = t(otu_table(feces_t6.css)), distance = "bray", try = 10, trymax = 100, autotransform = F)
plot_ordination(feces_t6.css, feces_t6.ord, color = "treatment_order") +
  theme_bw() +
  labs(title ="TIMEPOINT F") +
  stat_ellipse(aes(fill= treatment_order), geom="polygon", alpha = 0.25) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", size = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

feces_t6.adonis <- adonis2(feces_t6.dist ~ treatment_order, data = feces_t6.css.df, permutations = 9999)
feces_t6.adonis #NS

feces_t6.hclust <- hclust(feces_t6.dist)
plot(feces_t6.hclust)

feces_t6.dendro <- as.dendrogram(feces_t6.hclust)
feces_t6.dendro.data <- dendro_data(feces_t6.dendro, type = "rectangle")
feces_t6.dendro.data

## add in metadata columns
write.csv(feces_t6.dendro.data[["labels"]][["label"]], "feces_t6_dendro_sample_order.csv")

feces_t6_treatment_col <- data.frame(trt = c("AB",	"AB",	"AB",	"AB",	"BA",	"AB",	"BA",	"BA",	"BA",	"BA",	"BA",	"AB",	"BA",	"AB",	"AB",	"BA",	"BA",	"BA",	"BA",	"BA"))

feces_t6.dendro.data$labels <- cbind(feces_t6.dendro.data$labels, feces_t6_treatment_col)
feces_t6.dendro.data$labels

ggplot(feces_t6.dendro.data$segments) + theme_bw() +
  labs(y= "Ward's Distance", title = "TIMEPOINT F") + coord_flip() + scale_x_reverse() + scale_y_reverse() +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend), size =.75, lineend = "round", linejoin = "round") +
  geom_point(data = feces_t6.dendro.data$labels, aes(x,y, colour = trt, fill = trt, shape = trt),
             size = 8, stroke = 0.75, position = position_nudge(x = -0.02)) +
  #geom_text(data = all.dendro.data$labels, aes(x,y, label = day, colour= day), 
  #         size =3.8, position = position_nudge(y = -0.017)) +
  scale_shape_manual(values =c(15,15)) +
  theme(legend.position = "none",
        panel.border = element_blank(),
        plot.title = element_text(size = 30),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(size = 0.75),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 14, colour = "black"))

t6_sample_order <- c("HFTF-104F_S24",	"HFTF-124F_S115",	"HFTF-106F_S38",	"HFTF-122F_S101",	"HFTF-128F_S129",	"HFTF-103F_S17",	"HFTF-108A_S40",	"HFTF-120F_S87",	"HFTF-102F_S10",	"HFTF-132F_S150",	"HFTF-133F_S157",	"HFTF-123F_S108",	"HFTF-126F_S122",	"HFTF-114F_S70",	"HFTF-111F_S59",	"HFTF-105F_S31",	"HFTF-129F_S136",	"HFTF-110F_S52",	"HFTF-122A_S96",	"HFTF-130F_S143")

ggplot(ra_class_feces_melt, aes(x=Sample, y= Abundance, fill = Class)) +
  theme_bw() + coord_flip() + scale_x_reverse() + scale_y_reverse() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = ra_class_palette) +
  scale_x_discrete(limits = rev(t6_sample_order)) +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.75),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.75, lineend = "square", colour = "black"),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 14, colour = "black"))

## TIMEPOINT 7
feces_t7.css <- subset_samples(feces.css, pool_timepoint=="G")
feces_t7.css <- prune_taxa(taxa_sums(feces_t7.css) > 0, feces_t7.css)
feces_t7.css.df <- as(sample_data(feces_t7.css),"data.frame")
feces_t7.dist <- vegdist(t(otu_table(feces_t7.css)), method = "bray")
feces_t7.ord <- vegan::metaMDS(comm = t(otu_table(feces_t7.css)), distance = "bray", try = 10, trymax = 100, autotransform = F)
plot_ordination(feces_t7.css, feces_t7.ord, color = "treatment_order") +
  theme_bw() +
  labs(title ="TIMEPOINT G") +
  stat_ellipse(aes(fill= treatment_order), geom="polygon", alpha = 0.25) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", size = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

feces_t7.adonis <- adonis2(feces_t7.dist ~ treatment_order, data = feces_t7.css.df, permutations = 9999)
feces_t7.adonis #NS

feces_t7.hclust <- hclust(feces_t7.dist)
plot(feces_t7.hclust)

feces_t7.dendro <- as.dendrogram(feces_t7.hclust)
feces_t7.dendro.data <- dendro_data(feces_t7.dendro, type = "rectangle")
feces_t7.dendro.data

## add in metadata columns
write.csv(feces_t7.dendro.data[["labels"]][["label"]], "feces_t7_dendro_sample_order.csv")

feces_t7_treatment_col <- data.frame(trt = c("AB",	"AB",	"AB",	"AB",	"BA",	"AB",	"BA",	"BA",	"BA",	"BA",	"BA",	"AB",	"BA",	"AB",	"AB",	"BA",	"BA",	"BA",	"BA",	"BA"))

feces_t7.dendro.data$labels <- cbind(feces_t7.dendro.data$labels, feces_t7_treatment_col)
feces_t7.dendro.data$labels

ggplot(feces_t7.dendro.data$segments) + theme_bw() +
  labs(y= "Ward's Distance", title = "TIMEPOINT G") + coord_flip() + scale_x_reverse() + scale_y_reverse() +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend), size =.75, lineend = "round", linejoin = "round") +
  geom_point(data = feces_t7.dendro.data$labels, aes(x,y, colour = trt, fill = trt, shape = trt),
             size = 8, stroke = 0.75, position = position_nudge(x = -0.02)) +
  #geom_text(data = all.dendro.data$labels, aes(x,y, label = day, colour= day), 
  #         size =3.8, position = position_nudge(y = -0.017)) +
  scale_shape_manual(values =c(15,15)) +
  theme(legend.position = "none",
        panel.border = element_blank(),
        plot.title = element_text(size = 30),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(size = 0.75),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 14, colour = "black"))

t7_sample_order <- c("HFTF-122B_S97",	"HFTF-124G_S116",	"HFTF-108G_S46",	"HFTF-103G_S18",	"HFTF-105G_S32",	"HFTF-111G_S60",	"HFTF-104G_S25",	"HFTF-114G_S71",	"HFTF-102G_S11",	"HFTF-126G_S123",	"HFTF-120G_S88",	"HFTF-130G_S144",	"HFTF-110G_S53",	"HFTF-129G_S137",	"HFTF-122G_S102",	"HFTF-128G_S130",	"HFTF-132G_S151",	"HFTF-133G_S158",	"HFTF-106G_S39",	"HFTF-123G_S109")

ggplot(ra_class_feces_melt, aes(x=Sample, y= Abundance, fill = Class)) +
  theme_bw() + coord_flip() + scale_x_reverse() + scale_y_reverse() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = ra_class_palette) +
  scale_x_discrete(limits = rev(t7_sample_order)) +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.75),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.75, lineend = "square", colour = "black"),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 14, colour = "black"))


#### MEAT
meat.css <- subset_samples(data_noSNP.css, sample_type=="Meat Rinsate")
sum(taxa_sums(meat.css)==0)
meat.css <- prune_taxa(taxa_sums(meat.css) > 0, meat.css)
meat.css.df <- as(sample_data(meat.css),"data.frame")
meat.dist <- vegdist(t(otu_table(meat.css)), method = "bray")
meat.ord <- vegan::metaMDS(comm = t(otu_table(meat.css)), distance = "bray", try = 20, trymax = 50, autotransform = F)

plot_ordination(meat.css, meat.ord, color = "treatment") +
  theme_bw() + facet_wrap(~meat_type) +
  stat_ellipse(geom = "polygon", aes(fill= treatment), alpha = 0.25) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", size = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

meat.adonis <- adonis2(meat.dist ~ treatment*meat_type, meat.css.df, permutations = 9999)
meat.adonis # meat type different, treatment NS

meat.hclust <- hclust(meat.dist)
plot(meat.hclust)
meat.dendro <- as.dendrogram(meat.hclust)
meat.dendro.data <- dendro_data(meat.dendro, type = "rectangle")
meat.dendro.data

### need to add some metadata columns for plotting, easier to export and do in excel
#write.csv(meat.dendro.data[["labels"]][["label"]], "meat_dendro_sample_order.csv")

meat_type_col <- data.frame(meat_type =c ("Tloin",	"Tloin",	"Tloin",	"Tloin",	"Ground Beef",	"Tloin",	"Ground Beef",	"Tloin",	"Sirloin",	"Sirloin",	"Sirloin",	"Sirloin",	"Sirloin",	"Sirloin",	"Ground Beef",	"Ground Beef",	"Sirloin",	"Sirloin",	"Sirloin",	"Sirloin",	"Ground Beef",	"Ground Beef",	"Sirloin",	"Sirloin",	"Ground Beef",	"Sirloin",	"Tloin",	"Ground Beef",	"Ground Beef",	"Sirloin",	"Sirloin",	"Ground Beef",	"Sirloin",	"Ground Beef",	"Sirloin",	"Ground Beef",	"Ground Beef",	"Ground Beef",	"Sirloin",	"Ground Beef",	"Ground Beef",	"Ground Beef",	"Ground Beef",	"Ground Beef",	"Sirloin",	"Sirloin",	"Ground Beef",	"Ground Beef",	"Ground Beef",	"Sirloin",	"Sirloin",	"Ground Beef",	"Ground Beef",	"Ground Beef",	"Ground Beef",	"Ground Beef",	"Ground Beef",	"Ground Beef",	"Ground Beef",	"Tloin",	"Tloin",	"Tloin"))
meat_trt_col <- data.frame(meat_trt =c("Conventional",	"RWA",	"RWA",	"RWA",	"RWA",	"RWA",	"Conventional",	"RWA",	"Conventional",	"RWA",	"Conventional",	"Conventional",	"RWA",	"RWA",	"RWA",	"RWA",	"RWA",	"Conventional",	"Conventional",	"RWA",	"Conventional",	"Conventional",	"RWA",	"RWA",	"Conventional",	"RWA",	"Conventional",	"RWA",	"RWA",	"Conventional",	"Conventional",	"RWA",	"RWA",	"RWA",	"Conventional",	"RWA",	"Conventional",	"RWA",	"RWA",	"RWA",	"Conventional",	"RWA",	"Conventional",	"RWA",	"Conventional",	"Conventional",	"Conventional",	"Conventional",	"Conventional",	"Conventional",	"RWA",	"Conventional",	"RWA",	"Conventional",	"RWA",	"Conventional",	"RWA",	"RWA",	"Conventional",	"Conventional",	"Conventional",	"Conventional"))

meat.dendro.data$labels <- cbind(meat.dendro.data$labels, meat_type_col, meat_trt_col)
meat.dendro.data$labels

ggplot(meat.dendro.data$segments) + theme_bw() +
  labs(y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y, xend=xend,yend=yend),
               size = 0.8, lineend = "round", linejoin = "round") +
  geom_point(data = meat.dendro.data$labels, aes(x,y, colour = meat_type, shape = meat_trt),
             size = 8, position = position_nudge(y = -0.02)) +
  scale_colour_manual(values = c("dodgerblue3","darkorange3","forestgreen")) +
  scale_shape_manual(values =c (15,16)) +
  theme(#legend.position = "none",
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(size = 0.75),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black"))
  


dendro_sample_order <- c("HFTR-CT12_S206",	"HFTR-RT13_S213",	"HFTR-RT16_S225",	"HFTR-RT15_S221",	"HFTR-RG13_S211",	"HFTR-RT14_S217",	"HFTR-CG12_S208",	"HFTR-RT12_S209",	"HFTR-CS1_S162",	"HFTR-RS10_S200",	"HFTR-CS2_S166",	"HFTR-CS5_S178",	"HFTR-RS11_S204",	"HFTR-RS5_S180",	"HFTR-RG11_S205",	"HFTR-RG2_S169",	"HFTR-RS3_S172",	"HFTR-CS6_S182",	"HFTR-CS9_S194",	"HFTR-RS4_S176",	"HFTR-CG14_S216",	"HFTR-CG13_S212",	"HFTR-RS7_S188",	"HFTR-RS6_S184",	"HFTR-CG15_S220",	"HFTR-RS9_S196",	"HFTR-CT13_S210",	"HFTR-RG10_S201",	"HFTR-RG8_S193",	"HFTR-CS11_S202",	"HFTR-CS8_S190",	"HFTR-RG1_S165",	"HFTR-RS1_S164",	"HFTR-RG12_S207",	"HFTR-CS10_S198",	"HFTR-RG6_S185",	"HFTR-CG7_S187",	"HFTR-RG15_S219",	"HFTR-RS2_S168",	"HFTR-RG9_S197",	"HFTR-CG11_S203",	"HFTR-RG5_S181",	"HFTR-CG3_S171",	"HFTR-RG4_S177",	"HFTR-CS3_S170",	"HFTR-CS4_S174",	"HFTR-CG16_S224",	"HFTR-CG2_S167",	"HFTR-CG6_S183",	"HFTR-CS7_S186",	"HFTR-RS8_S192",	"HFTR-CG10_S199",	"HFTR-RG14_S215",	"HFTR-CG9_S195",	"HFTR-RG7_S189",	"HFTR-CG8_S191",	"HFTR-RG16_S223",	"HFTR-RG3_S173",	"HFTR-CG4_S175",	"HFTR-CT16_S222",	"HFTR-CT14_S214",	"HFTR-CT15_S218")

ggplot(ra_class_meat_melt, aes(x= SampleID, y= Abundance, fill = Class)) +
  geom_bar(stat = "summary", colour = "black") +
  theme_bw() +
  labs(y= "Relative Abundance (%)") +
  scale_fill_manual(values = class_meat_palette) +
  scale_x_discrete(limits = dendro_sample_order) +
  theme(#legend.position = "none",
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(size = 0.75),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black"))

