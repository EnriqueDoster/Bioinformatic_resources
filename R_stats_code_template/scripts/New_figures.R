
#
## Ordination feces_ endtrial
#

# Subset samples from pool_timepoint=="A" , using the CSS normalized feces counts

meat.css <- prune_taxa(taxa_sums(meat.css) > 0, meat.css)
meat.css.df <- as(sample_data(meat.css),"data.frame")
meat.dist <- vegdist(t(otu_table(meat.css)), method = "bray")
meat.ord <- vegan::metaMDS(comm = t(otu_table(meat.css)), distance = "bray", try = 10, trymax = 100, autotransform = F)

plt_ord_meat <- plot_ordination(meat.css, meat.ord, color = "specific_sample_type") +
  theme_bw() +
  labs(title ="Meat samples") +
  stat_ellipse(lty = 2) +
  #stat_ellipse(aes(fill= specific_sample_type), geom="polygon", alpha = 0.25) +
  geom_point(size = 4) +
  #labs(color="Meat type") +
  theme(legend.position = "right",
        legend.title = element_blank(),
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
plt_ord_meat 

# permanova test of beta diversity, by End_trial
meat.adonis <- adonis2(meat.dist ~ Meat_type, data = meat.css.df, permutations = 9999)
meat.adonis #NS

pairwise.adonis2(meat.dist ~ specific_sample_type, meat.css.df)
pairwise.adonis2(meat.dist ~ treatment, meat.css.df)

# This requires a directory called figures to be in your working directory
# All 3 lines have to be run for the file to be created.
png("Figures/meat_NMDS_fig.png", width = 1200, height = 800)
plt_ord_meat
dev.off()

# svg format is nice because it stores all graphic data and can be requested by some manuscripts
svg("Figures/meat_NMDS_fig.svg", width = 500, height = 500)
plt_ord_meat
dev.off()



###### All sample dendrogram and phylum barplots ######

# Bray curtis distance matrix did converge with all samples
# Had to use a different beta diversity index, euclidean, calculated on hellinger transformed counts
data.dist <- vegdist(decostand(t(otu_table(data.css)), "hell"), "bray") 
data.ord <- vegan::metaMDS(comm = t(data.dist), distance = "none", try = 10, trymax = 999, autotransform = F)
#plot_ordination(data.css, data.ord, color = "sample_type") 

plt_ord_by_sampletype <- plot_ordination(data.css, data.ord, color = "sample_type") +
  theme_bw() +
  labs(title ="Sample type") +
  stat_ellipse(aes(fill= sample_type), geom="polygon", alpha = 0.25) +
  geom_point(size = 4) +
  theme(legend.position = "bottom",
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
plt_ord_by_sampletype 

# permanova test of beta diversity, by End_trial
data.adonis <- adonis2(data.dist ~ sample_type, data = data.css.df, permutations = 9999)
data.adonis #NS

# This requires a directory called figures to be in your working directory
# All 3 lines have to be run for the file to be created.
png("Figures/Allsamples_NMDS_fig.png", width = 1200, height = 800)
plt_ord_by_sampletype
dev.off()

# svg format is nice because it stores all graphic data and can be requested by some manuscripts
svg("Figures/Allsamples_NMDS_fig.svg", width = 500, height = 500)
plt_ord_by_sampletype
dev.off()


##
### Clustering at ASV level
##

# Cluster using the function hclust() and default settings
data.hclust <- hclust(data.dist)
plot(data.hclust) # example plot

# Extract data as dendrogram
data.dendro <- as.dendrogram(data.hclust)
data.dendro.data <- dendro_data(data.dendro, type = "rectangle")
data.dendro.data #  this object contains $segments and $labels

# Sample names in order based on clustering
sample_names_data <- data.dendro.data$labels$label

# Add metadata
data.dendro.data$labels <- left_join(data.dendro.data$labels, map_file, by = c("label" = "sample.id"))

# Setup the data, so that the layout is inverted (this is more 
# "clear" than simply using coord_flip())
segment_data_data <- with(
  segment(data.dendro.data),
  data.frame(x = y, y = x, xend = yend, yend = xend))

# Use the dendrogram label data to position the gene labels
gene_pos_table_data <- with(
  data.dendro.data$labels, 
  data.frame(y_center = x, gene = as.character(label), x = y , trt = as.character(specific_sample_type), height = 1))

# Table to position the samples
sample_pos_table_data <- data.frame(sample = sample_names_data) %>% 
  dplyr::mutate(x_center = (1:dplyr::n()),  width = 1)

##
### Relative abundance bar plot
##

# Use phylum melted data and add gene locations
joined_ra_phylum_data_melt <- ra_phylum_melt %>%
  left_join(gene_pos_table_data, by = c("Sample" = "gene")) %>%
  left_join(sample_pos_table_data, by = c("Sample" = "sample")) 

# Calculate the mean relative abundance of each phylum taxa, sort by most abundant to least
factor_by_abund <- ra_phylum_melt %>%
  group_by(Phylum) %>%
  summarize_at( vars(Abundance), list(mean_phylum = mean)) %>%
  arrange(-mean_phylum)

# Use the sorted phylum levels to clean up the order of taxa in the relative abundance plots
joined_ra_phylum_data_melt$Phylum <- factor(joined_ra_phylum_data_melt$Phylum, levels = as.character(factor_by_abund$Phylum))

# Limits for the vertical axes
gene_axis_limits_data <- with(
  gene_pos_table_data, 
  c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))) +   0.1 * c(-1, 1) # extra spacing: 0.1


# Useful 20 colors + black and white (https://sashamaps.net/docs/resources/20-colors/)
#'#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#ffffff', '#000000'

#
## Relative abundance plot
#
plt_rel_phylum_data <- ggplot(joined_ra_phylum_data_melt, 
                                  aes(x = x_center, y = Abundance, fill = Phylum, 
                                      height = height, width = width)) + 
  coord_flip() +
  scale_fill_jama() +
  geom_bar(stat = "identity", colour = "black") +
  #scale_fill_manual(values = col_vector) + #use this if not using library(ggsci)
  scale_x_continuous(breaks = sample_pos_table_data$x_center, 
                     labels = sample_pos_table_data$sample, 
                     expand = c(0, 0)) + 
  labs(x = "", y = "Relative abundance") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.75),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.75, lineend = "square", colour = "black"),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, colour = "black"),
        plot.margin = unit(c(1, 0.2, 0.2, -0.7), "cm"), 
        panel.grid.minor = element_blank())
plt_rel_phylum_data

# Dendrogram plot
plt_dendr_phylum_data <- ggplot(segment_data_data) + 
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend), size =.75, lineend = "round", linejoin = "round") +
  geom_point(data = gene_pos_table_data, aes(x,y_center, colour = trt, fill = trt, shape = trt),
             size = 8, stroke = 0.75, position = position_nudge(x = -0.03)) +
  scale_y_continuous(breaks = gene_pos_table_data$y_center, 
                     labels = gene_pos_table_data$gene, 
                     limits = gene_axis_limits_data, 
                     expand = c(0, 0)) + 
  labs(x = "Ward's Distance", y = "", colour = "", size = "", title = "All samples") +
  scale_x_reverse() + 
  scale_shape_manual(values =c(15,15,15,15,15)) +
  theme_bw() + 
  guides(shape = "none", fill = "none") + 
  theme(legend.position = "bottom",
        panel.border = element_blank(),
        plot.title = element_text(size = 30),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(size = 0.75),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, colour = "black"))
plt_dendr_phylum_data

phylum_data_plots <- plot_grid(plt_dendr_phylum_data, plt_rel_phylum_data, align = 'h', rel_widths = c(0.75, 1.25))

# This requires a directory called figures to be in your working directory
# All 3 lines have to be run for the file to be created.
png("Figures/Large_All_Phylum_dendro_relabun_fig.png", width = 1600, height = 1000)
phylum_data_plots
dev.off()

# svg format is nice because it stores all graphic data and can be requested by some manuscripts
svg("Figures/All_Phylum_dendro_relabun_fig.svg", width = 500, height = 500)
phylum_data_plots
dev.off()

