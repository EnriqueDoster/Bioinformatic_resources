# AMR beta diversity by shared household status
#############################################################################################
####################################   CSS TRANSFORM    ###################################
SNV_data_noSNP <- prune_taxa(taxa_sums(SNV_data_noSNP) > 0, SNV_data_noSNP)
any(taxa_sums(SNV_data_noSNP)==0) # QUADRUPLE CHECKING - nope good.

SNV_data_noSNP.css <- phyloseq_transform_css(SNV_data_noSNP, log = F)
SNV_data_noSNP.css.df <- as(sample_data(SNV_data_noSNP.css), "data.frame")




###
####
# Only subset "Habitating" ######
####
### 

Habitating_SNV_data_noSNP.css <- subset_samples(SNV_data_noSNP.css, Habitating!="")

Habitating_SNV_data_noSNP.css <- prune_taxa(taxa_sums(Habitating_SNV_data_noSNP.css) > 0, Habitating_SNV_data_noSNP.css)


Habitating_SNV_data_noSNP.css.df <- as(sample_data(Habitating_SNV_data_noSNP.css),"data.frame")
Habitating_SNV_data_noSNP.css.dist <- vegdist(t(otu_table(Habitating_SNV_data_noSNP.css)), method = "bray")
Habitating_SNV_data_noSNP.css.ord <- vegan::metaMDS(comm = t(otu_table(Habitating_SNV_data_noSNP.css)), distance = "bray", try = 10, trymax = 100, autotransform = F)
plt_Habitating <- plot_ordination(Habitating_SNV_data_noSNP.css, Habitating_SNV_data_noSNP.css.ord, color = "Habitating", shape = "Habitating") +
  theme_bw() +
  labs(title ="Habitating") +
  stat_ellipse(aes(fill= Habitating), geom="polygon", alpha = 0.25) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", size = 1.0),
        strip.text = element_text(size =22, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 24),
        axis.ticks = element_line(colour = "black", linewidth = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
plt_Habitating 

# permanova test of beta diversity, by Habitating
Habitating_SNV_data_noSNP.css.adonis <- adonis2(Habitating_SNV_data_noSNP.css.dist ~Habitating, data = Habitating_SNV_data_noSNP.css.df, permutations = 9999)
Habitating_SNV_data_noSNP.css.adonis #


## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
Habitating.permanova <- pairwise.adonis(Habitating_SNV_data_noSNP.css.dist, Habitating_SNV_data_noSNP.css.df$Habitating, perm = 9999, p.adjust.m = "BH")
Habitating.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
Habitating.disper <- betadisper(Habitating_SNV_data_noSNP.css.dist, Habitating_SNV_data_noSNP.css.df$Habitating)
Habitating.permdisp <- permutest(Habitating.disper, permutations = 9999, pairwise = T)
Habitating.permdisp # looks like a few are significant



###
####
# Only subset "Co.Habitating" #####
####
### 

Co.Habitating_SNV_data_noSNP.css <- subset_samples(SNV_data_noSNP.css, Co.Habitating!="")

Co.Habitating_SNV_data_noSNP.css <- prune_taxa(taxa_sums(Co.Habitating_SNV_data_noSNP.css) > 0, Co.Habitating_SNV_data_noSNP.css)


Co.Habitating_SNV_data_noSNP.css.df <- as(sample_data(Co.Habitating_SNV_data_noSNP.css),"data.frame")
Co.Habitating_SNV_data_noSNP.css.dist <- vegdist(t(otu_table(Co.Habitating_SNV_data_noSNP.css)), method = "bray")
Co.Habitating_SNV_data_noSNP.css.ord <- vegan::metaMDS(comm = t(otu_table(Co.Habitating_SNV_data_noSNP.css)), distance = "bray", try = 10, trymax = 100, autotransform = F)
plt_Co.Habitating <- plot_ordination(Co.Habitating_SNV_data_noSNP.css, Co.Habitating_SNV_data_noSNP.css.ord, color = "Co.Habitating", shape = "Co.Habitating") +
  theme_bw() +
  labs(title ="Co.Habitating") +
  stat_ellipse(aes(fill= Co.Habitating), geom="polygon", alpha = 0.25) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", size = 1.0),
        strip.text = element_text(size =22, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 24),
        axis.ticks = element_line(colour = "black", linewidth = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
plt_Co.Habitating 

# permanova test of beta diversity, by Co.Habitating
Co.Habitating_SNV_data_noSNP.css.adonis <- adonis2(Co.Habitating_SNV_data_noSNP.css.dist ~Co.Habitating, data = Co.Habitating_SNV_data_noSNP.css.df, permutations = 9999)
Co.Habitating_SNV_data_noSNP.css.adonis #


## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
Co.Habitating.permanova <- pairwise.adonis(Co.Habitating_SNV_data_noSNP.css.dist, Co.Habitating_SNV_data_noSNP.css.df$Co.Habitating, perm = 9999, p.adjust.m = "BH")
Co.Habitating.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
Co.Habitating.disper <- betadisper(Co.Habitating_SNV_data_noSNP.css.dist, Co.Habitating_SNV_data_noSNP.css.df$Co.Habitating)
Co.Habitating.permdisp <- permutest(Co.Habitating.disper, permutations = 9999, pairwise = T)
Co.Habitating.permdisp # looks like a few are significant


#
##
# Calculating average pairwise distances for each group ####
##
#

## Indices of samples from shared households
shared_indices <- which(Habitating_SNV_data_noSNP.css.df$Habitating == "Shared")
# Extract pairwise distances for shared households
shared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_indices, shared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_distances[lower.tri(shared_distances)])


## Indices of samples from Shared_1
shared_1_indices <- which(Habitating_SNV_data_noSNP.css.df$Co.Habitating == "Shared_1" )
# Extract pairwise distances for shared households
shared_1_indices <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_1_indices, shared_1_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_1_indices[lower.tri(shared_1_indices)])


## Indices of samples from Shared_2
shared_2_indices <- which(Habitating_SNV_data_noSNP.css.df$Co.Habitating == "Shared_2")
# Extract pairwise distances for shared households
shared_2_indices <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[shared_2_indices, shared_2_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(shared_2_indices[lower.tri(shared_2_indices)])

## Nonshared ###
# Indices of samples from nonshared households
nonshared_indices <- which(Habitating_SNV_data_noSNP.css.df$Habitating == "Control")
# Extract pairwise distances for nonshared households
nonshared_distances <- as.matrix(Habitating_SNV_data_noSNP.css.dist)[nonshared_indices, nonshared_indices]
# Calculate the average distance, excluding the diagonal (which are zeros)
mean(nonshared_distances[lower.tri(nonshared_distances)])




