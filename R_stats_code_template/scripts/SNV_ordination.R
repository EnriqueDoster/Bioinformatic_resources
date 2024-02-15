# Beta diversity
#############################################################################################
##############################         BETA DIVERSITY         #########################@#####
#############################################################################################
#############################################################################################

###
##### Combined order - All samples #####
###
##
#

#### CSS transformation of counts
SNV_data_noSNP.css <- phyloseq_transform_css(SNV_data_noSNP, log = F)

SNV_data_noSNP.css <-  prune_taxa(taxa_sums(SNV_data_noSNP.css) > 0, SNV_data_noSNP.css)

#### create d.f. for beta-diversity metadata
SNV_data_noSNP.css.df <- as(sample_data(SNV_data_noSNP.css),"data.frame")

# ordinate it based on Bray-Curtis
SNV_data_noSNP.ord <- vegan::metaMDS(comm = t(otu_table(SNV_data_noSNP.css)), try = 50, trymax = 500, distance = "bray", autotransform = F)

# plot the ordination with 95% confidence ellipses coloured by groups
plot_ordination(SNV_data_noSNP.css, SNV_data_noSNP.ord, type = "samples", color = "Combined_order") +
  theme_bw() +
  geom_point(size = 2.5) +
  stat_ellipse(geom = "polygon", aes(fill=Combined_order), level = 0.95, lty =2, alpha= 0.3) +
  #scale_fill_manual(values = group_palette) +
  #scale_colour_manual(values = group_palette) +
  theme(plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        panel.border = element_rect(colour = "black", size =0.75),
        legend.key.size = unit(2, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 24, vjust = 1.75),
        axis.title.x = element_text(size = 24, vjust = -1.5))

## STATS BETWEEN ALL GROUPS
# create distance matrix we use to ordinate (Bray-Curtis)
SNV_data_noSNP.css.dist <- vegdist(t(otu_table(SNV_data_noSNP.css)), method = "bray")

## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
all.groups.permanova <- pairwise.adonis(SNV_data_noSNP.css.dist, SNV_data_noSNP.css.df$Combined_order, perm = 9999, p.adjust.m = "BH")
all.groups.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
all.groups.disper <- betadisper(SNV_data_noSNP.css.dist, SNV_data_noSNP.css.df$Combined_order)
all.groups.permdisp <- permutest(all.groups.disper, permutations = 9999, pairwise = T)
all.groups.permdisp # looks like a few are significant

###
##### Combined order - Fecal samples #####
###
##
#

Fecal.SNV_data_noSNP.css <- subset_samples(SNV_data_noSNP.css, sample_type =="Fecal")

Fecal.SNV_data_noSNP.css <-  prune_taxa(taxa_sums(Fecal.SNV_data_noSNP.css) > 0, Fecal.SNV_data_noSNP.css)


#### create d.f. for beta-diversity metadata
Fecal.SNV_data_noSNP.css.df <- as(sample_data(Fecal.SNV_data_noSNP.css),"data.frame")

# ordinate it based on Bray-Curtis
Fecal.SNV_data_noSNP.ord <- vegan::metaMDS(comm = t(otu_table(Fecal.SNV_data_noSNP.css)), try = 50, trymax = 500, distance = "bray", autotransform = F)

# plot the ordination with 95% confidence ellipses coloured by groups
plot_ordination(Fecal.SNV_data_noSNP.css, Fecal.SNV_data_noSNP.ord, type = "samples", color = "Combined_order") +
  theme_bw() +
  geom_point(size = 2.5) +
  stat_ellipse(geom = "polygon", aes(fill=Combined_order), level = 0.95, lty =2, alpha= 0.3) +
  #scale_fill_manual(values = group_palette) +
  #scale_colour_manual(values = group_palette) +
  theme(plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        panel.border = element_rect(colour = "black", size =0.75),
        legend.key.size = unit(2, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 24, vjust = 1.75),
        axis.title.x = element_text(size = 24, vjust = -1.5))

## STATS BETWEEN ALL GROUPS
# create distance matrix we use to ordinate (Bray-Curtis)
Fecal.SNV_data_noSNP.css.dist <- vegdist(t(otu_table(Fecal.SNV_data_noSNP.css)), method = "bray")

## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
Fecal.groups.permanova <- pairwise.adonis(Fecal.SNV_data_noSNP.css.dist, SNV_data_noSNP.css.df$Combined_order, perm = 9999, p.adjust.m = "BH")
Fecal.groups.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
Fecal.groups.disper <- betadisper(Fecal.SNV_data_noSNP.css.dist, Fecal.SNV_data_noSNP.css.df$Combined_order)
Fecal.groups.permdisp <- permutest(Fecal.groups.disper, permutations = 9999, pairwise = T)
Fecal.groups.permdisp # looks like a few are significant