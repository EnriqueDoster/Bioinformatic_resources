# PCV beta diversity


#############################################################################################
####################################   CSS TRANSFORM    ###################################

mhpp_PSV.ps <- prune_taxa(taxa_sums(mhpp_PSV.ps) > 0, mhpp_PSV.ps)
any(taxa_sums(mhpp_PSV.ps)==0) # QUADRUPLE CHECKING - nope good.

mhpp_PSV.ps.css <- phyloseq_transform_css(mhpp_PSV.ps, log = F)
mhpp_PSV.ps.css.df <- as(sample_data(mhpp_PSV.ps.css), "data.frame")


###
####
##### All samples- Beta diversity #######
####
### 

# Bray curtis distance matrix did converge with all samples
# Had to use a different beta diversity index, euclidean, calculated on hellinger transformed counts
data.dist <- vegdist(decostand(t(otu_table(mhpp_PSV.ps.css)), "hell"), "euclidean") 
data.ord <- vegan::metaMDS(comm = t(data.dist), distance = "none", try = 10, trymax = 999, autotransform = F)

plt_ord_by_Sample_type <- plot_ordination(mhpp_PSV.ps.css, data.ord, color = "Sample_Type") +
  theme_bw() +
  labs(title ="Sample type") +
  stat_ellipse(aes(fill= Sample_Type), geom="polygon", alpha = 0.25) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", size = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", linewidth = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
plt_ord_by_Sample_type

# permanova test of beta diversity, by Dilution
data.adonis <- adonis2(data.dist ~ Sample_Type, data = mhpp_PSV.ps.css.df, permutations = 9999)
data.adonis # 7e-04, R2 =0.13205 

## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
Sample_Type.permanova <- pairwise.adonis(data.dist, mhpp_PSV.ps.css.df$Sample_Type, perm = 9999, p.adjust.m = "BH")
Sample_Type.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
Sample_Type.disper <- betadisper(data.dist, mhpp_PSV.ps.css.df$Sample_Type)
Sample_Type.permdisp <- permutest(Sample_Type.disper, permutations = 9999, pairwise = T)
Sample_Type.permdisp # looks like a few are significant


## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
Group.permanova <- pairwise.adonis(data.dist, mhpp_PSV.ps.css.df$Group, perm = 9999, p.adjust.m = "BH")
Group.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
Group.disper <- betadisper(data.dist, mhpp_PSV.ps.css.df$Group)
Group.permdisp <- permutest(Group.disper, permutations = 9999, pairwise = T)
Group.permdisp # looks like a few are significant



### Group plot ####
plt_ord_by_Group <- plot_ordination(mhpp_PSV.ps.css, data.ord, color = "Group") +
  theme_bw() +
  labs(title ="Group") +
  stat_ellipse(aes(fill= Group), geom="polygon", alpha = 0.25) +
  theme(legend.position = "right",
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
plt_ord_by_Group


