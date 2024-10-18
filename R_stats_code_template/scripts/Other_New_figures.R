# PLOT THE DENDROGRAMS
# now we can use ggplot and can add colours, shapes, etc based on our metadata
microbiome_dendro_plot <- ggplot(data.dendro.data$segments) +
  theme_minimal() +
  labs(title = "Microbiome Hierarchal Clustering", y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = data.dendro.data$labels, aes(x=x,y=y, colour =sample_type), size = 1, shape = 15) +
  theme(plot.title = element_text(size = 20),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
microbiome_dendro_plot

# use the dendrogram label data to position the samples in the right order
# make object of microbiome sample names
microbiome_sample_names <- data.dendro.data$labels$label
# plot the RA barplot by individual samples
microbiome_ra_dendro_plot <- ggplot(ra_phylum_melt, aes(x= sample.id, y= Abundance, fill = Phylum)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = microbiome_sample_names) +
  theme(plot.title = element_text(size = 20),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
microbiome_ra_dendro_plot
# combine this plot with the dendrogram we made earlier
plot_grid(microbiome_dendro_plot, microbiome_ra_dendro_plot, align = "v", ncol = 1)
