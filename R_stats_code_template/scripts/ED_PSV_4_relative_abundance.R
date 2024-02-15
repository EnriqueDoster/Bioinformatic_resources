# Relative abundance plots
library(microViz)
library(ggsci) # if not using this, make sure to switch out the palatte in the figures
library(metagMisc)
library(dplyr)
library(ggdendro)
library(cowplot)

#############################################################################################
####################################   CSS TRANSFORM    ###################################

mhpp_PSV.ps <- prune_taxa(taxa_sums(mhpp_PSV.ps) > 0, mhpp_PSV.ps)
any(taxa_sums(mhpp_PSV.ps)==0) # QUADRUPLE CHECKING - nope good.

mhpp_PSV.ps.css <- phyloseq_transform_css(mhpp_PSV.ps, log = F)
mhpp_PSV.ps.css.df <- as(sample_data(mhpp_PSV.ps.css), "data.frame")

#############################################################################################
##############################         RELA ABUNDANCE         ###############################
#############################################################################################
#############################################################################################

rel_abund <- transform_sample_counts(mhpp_PSV.ps.css, function(x) {x/sum(x)}*100)

rel_abund_melt <- psmelt(rel_abund)




######## Re-label low abundance taxa into single category

# convert Phylum to a character vector from a factor because R
rel_abund_melt$PSV <- as.character(rel_abund_melt$OTU)

# calculate percentage across all counts
all_sample_factor_by_abund <- rel_abund_melt %>%
  dplyr::group_by(PSV) %>%
  dplyr::summarize(median_PSV_perc = median(Abundance)) %>%
  arrange(-median_PSV_perc)

#write.csv(all_sample_factor_by_abund,"figures/all_samples_mhppclass_medians_bysample.csv")


# find Phyla whose rel. abund. is less than 1%
remainder <- all_sample_factor_by_abund[all_sample_factor_by_abund$median_PSV_perc <= 0.5,]$PSV

# change the name of taxa to whose rel. abund. is less than 0.5% to "Low abundance phyla (<1%)"
rel_abund_melt[rel_abund_melt$PSV %in% remainder,]$PSV <- 'Low abundance PSV (<0.5%)'

rel_abund_melt_agg <- rel_abund_melt %>%
  group_by(Sample, PSV) %>%  # Group by both Sample and PSV
  summarise(Abundance = sum(Abundance)) %>%  # Sum the Abundance for each group
  ungroup()


# Determine the number of taxa
length(unique(rel_abund_melt_agg$PSV))


##
### Clustering at ASV level ####
##

# Cluster using the function hclust() and default settings
ps_mhpp_PSV.dist <- vegdist(t(otu_table(mhpp_PSV.ps.css)), method = "bray")
ps_mhpp_PSV.hclust <- hclust(ps_mhpp_PSV.dist)
plot(ps_mhpp_PSV.hclust) # example plot

# Extract data as dendrogram
ps_mhpp_PSV.dendro <- as.dendrogram(ps_mhpp_PSV.hclust)
ps_mhpp_PSV.dendro.data <- dendro_data(ps_mhpp_PSV.dendro, type = "rectangle")
ps_mhpp_PSV.dendro.data #  this object contains $segments and $labels

# Sample names in order based on clustering
sample_names_ps_mhpp_PSV <- ps_mhpp_PSV.dendro.data$labels$label

  
# Add metadata
ps_mhpp_PSV.dendro.data$labels <- left_join(ps_mhpp_PSV.dendro.data$labels, sample_data(mhpp_PSV.ps.css), by = c("label" = "sample_id"))

# Setup the data, so that the layout is inverted (this is more 
# "clear" than simply using coord_flip())
segment_data_ps_mhpp_PSV <- with(
  segment(ps_mhpp_PSV.dendro.data),
  data.frame(x = y, y = x, xend = yend, yend = xend))

# Use the dendrogram label data to position the gene labels
gene_pos_table_ps_mhpp_PSV <- with(
  ps_mhpp_PSV.dendro.data$labels, 
  data.frame(y_center = x, gene = as.character(label), x = y , Group = as.character(Group), Sample_Type = as.character(Sample_Type) , height = 1))


# Table to position the samples
sample_pos_table_ps_mhpp_PSV <- data.frame(sample = sample_names_ps_mhpp_PSV) %>%
  dplyr::mutate(x_center = (1:dplyr::n()),  width = 1)

##
######## Relative abundance bar plot #########
##

# Use class melted data and add gene locations
joined_ra_PSV_ps_PSV_melt <- rel_abund_melt_agg %>%
  left_join(gene_pos_table_ps_mhpp_PSV, by = c("Sample" = "gene")) %>%
  left_join(sample_pos_table_ps_mhpp_PSV, by = c("Sample" = "sample")) 

# Calculate the mean relative abundance of each class taxa, sort by most abundant to least


factor_by_abund <- rel_abund_melt %>%
  dplyr::group_by(PSV) %>%
  dplyr::summarize(median_PSV = median(Abundance)) %>%
  arrange(-median_PSV)

# Use the sorted class levels to clean up the order of taxa in the relative abundance plots
joined_ra_PSV_ps_PSV_melt$PSV <- factor(joined_ra_PSV_ps_PSV_melt$PSV, levels = as.character(factor_by_abund$PSV))

# Limits for the vertical axes
gene_axis_limits_ps_class <- with(
  gene_pos_table_ps_mhpp_PSV, 
  c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))) +   0.1 * c(-1, 1) # extra spacing: 0.1


# Useful 20 colors + black and white (https://sashamaps.net/docs/resources/20-colors/)
#'#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#ffffff', '#000000'

# Option to change name
#joined_ra_PSV_ps_PSV_melt <- joined_ra_PSV_ps_PSV_melt %>%
#  mutate(PSV = str_replace_all(PSV, "Mannheimia_haemolytica_strain_", "Mh "))

factor_by_abund <- joined_ra_PSV_ps_PSV_melt %>%
  dplyr::group_by(PSV) %>%
  dplyr::summarize(median_PSV = median(Abundance)) %>%
  arrange(-median_PSV)

# Use the sorted class levels to clean up the order of taxa in the relative abundance plots
joined_ra_PSV_ps_PSV_melt$PSV <- factor(joined_ra_PSV_ps_PSV_melt$PSV, levels = as.character(factor_by_abund$PSV))



#
## Relative abundance plot
#
plt_rel_PSV <- ggplot(joined_ra_PSV_ps_PSV_melt, 
                                 aes(x = x_center, y = Abundance, fill = PSV, 
                                     height = height, width = width)) + 
  coord_flip() +
  geom_bar(stat = "identity", colour = "black") +
  #scale_fill_manual(values = col_vector) + #use this if not using color palette from ggthemes of ggsci
  scale_fill_tableau(palette = "Tableau 20")+
  scale_x_continuous(breaks = sample_pos_table_ps_mhpp_PSV$x_center, 
                     labels = sample_pos_table_ps_mhpp_PSV$sample, 
                     expand = c(0, 0)) + 
  # For the y axis, alternatively set the labels as: gene_position_table$gene
  #scale_y_continuous(breaks = gene_pos_table[, "y_center"], labels = rep("", nrow(gene_pos_table)),limits = gene_axis_limits, expand = c(0, 0)) + 
  labs(x = "", y = "Relative abundance") +
  theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.75),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.75, lineend = "square", colour = "black"),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, colour = "black"),
        plot.margin = unit(c(1, 0.2, 0.2, -0.7), "cm"), 
        panel.grid.minor = element_blank())
plt_rel_PSV

# Dendrogram plot
plt_dendr_PSV <- ggplot(segment_data_ps_mhpp_PSV) + 
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend), size=.75, lineend="round", linejoin="round") +
  geom_point(data=gene_pos_table_ps_mhpp_PSV, aes(x, y_center, colour=Sample_Type, shape=Sample_Type),
             size=8, stroke=0.75, position=position_nudge(x=-0.03)) +
  scale_y_continuous(breaks=gene_pos_table_ps_mhpp_PSV$y_center, 
                     labels=gene_pos_table_ps_mhpp_PSV$gene, 
                     limits=gene_axis_limits_ps_class, 
                     expand=c(0, 0)) + 
  labs(x="Ward's Distance", y="", title="") +
  scale_x_reverse() + 
  scale_colour_manual(name="Sample type", values=c("#B2182B","#F4A582","#4393C3")) +  # Set legend title and colors
  scale_shape_manual(values=c(15,15,15,15,15,15), guide="none") +  
  theme_bw() + 
  theme(legend.position="right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        panel.border=element_blank(),
        plot.title=element_text(size=30),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x=element_line(linewidth=0.75),
        axis.title.x=element_text(size=14),
        axis.text.x=element_text(size=12, colour="black"))
plt_dendr_PSV

#plots_dendro_rel_PSV <- plot_grid(plt_dendr_PSV, plt_rel_PSV, align = 'h', rel_widths = c(0.5, 1.5))


# Extract the legends from each plot
legend_dendr_PSV <- get_legend(plt_dendr_PSV)
legend_rel_PSV <- get_legend(plt_rel_PSV)

# Remove the legends from the original plots
plt_dendr_PSV <- plt_dendr_PSV + theme(legend.position = "none")
plt_rel_PSV <- plt_rel_PSV + theme(legend.position = "none")

combined_legend <- ggdraw() +
  draw_plot(legend_dendr_PSV, x = -0.182, y = 0.7, width = 0.9, height = 0.3) +  # Adjust the x value to shift the top legend to the left
  draw_plot(legend_rel_PSV, x = 0, y = 0, width = 1, height = 0.7)

plots_dendro_rel_PSV <- plot_grid(
  plt_dendr_PSV, 
  plt_rel_PSV, 
  combined_legend,
  align = 'h', 
  ncol = 3,
  rel_widths = c(0.3, 1.1, 0.3)
)

plots_dendro_rel_PSV

#ggsave("mh_PSV_dendrogram_rel_abun_by_Sample_type.jpg", plot=last_plot(), device="jpeg", path="figures/", scale = 1, units = "in", dpi = 600)



# # This requires a directory called figures to be in your working directory
# # All 3 lines have to be run for the file to be created.
png("figures/mh_PSV_dendrogram_rel_abun_by_Sample_type.png", width = 1200, height = 800)
plots_dendro_rel_PSV
dev.off()

