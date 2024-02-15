# PSV 

# Calculate alpha diversity and update metadata file in phyloseq object with new results
mhpp_alpha_div <- estimate_richness(mhpp_PSV.ps, measures = c("Observed","Shannon","Simpson","InvSimpson"))

# Making sure row names match with phyloseq objects
rownames(mhpp_alpha_div) <- sample_names(mhpp_PSV.ps)

# Extract the existing sample data from the phyloseq object
existing_sample_data <- sample_data(mhpp_PSV.ps)

# Merge alpha diversity data with the existing sample data
new_sample_data <- merge(existing_sample_data, mhpp_alpha_div, by = "row.names", all = TRUE)
row.names(new_sample_data) <- new_sample_data$Row.names
new_sample_data$Row.names <- NULL

# Update the sample data in the phyloseq object
sample_data(mhpp_PSV.ps) <- sample_data(new_sample_data)



# Plot richness and diversity
# Richness by sample type
ggplot(sample_data(mhpp_PSV.ps), mapping=aes(x= Sample_Type, y = Observed, fill = Sample_Type, color = Sample_Type)) +
  theme_bw() + 
  labs(y = "Observed PSVs") +
  #labs(x = "Sample Type") +
  #labs(title = "Unique Mh PSVs Identified") +
  #scale_fill_discrete(labels = c("Nasal", "Rope", "Waterbowl")) +
  #scale_color_discrete(labels = c("Nasal", "Rope", "Waterbowl")) +
  geom_boxplot(alpha=0.35, show.legend = F) +
  geom_jitter(width = .15) +
  geom_point(alpha = .05, shape = ".") +
  geom_smooth() +
  theme(legend.position = "none",
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        legend.title = element_text("Sample Type", size = 10),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 18),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
ggsave("figures/mh_PSV_richness_barplot.jpg", plot=last_plot(), device="jpeg", path="figures/", scale = 1, units = "in", dpi = 600)

# Richness by group
ggplot(sample_data(mhpp_PSV.ps), mapping=aes(x= Group, y = Observed, fill = Group, color = Group)) +
  theme_bw() + 
  labs(y = "Observed PSVs") +
  #labs(x = "Sample Type") +
  #labs(title = "Unique Mh PSVs Identified") +
  #scale_fill_discrete(labels = c("Nasal", "Rope", "Waterbowl")) +
  #scale_color_discrete(labels = c("Nasal", "Rope", "Waterbowl")) +
  geom_boxplot(alpha=0.35, show.legend = F) +
  geom_jitter(width = .15) +
  geom_point(alpha = .05, shape = ".") +
  geom_smooth() +
  theme(legend.position = "none",
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        legend.title = element_text("Sample Type", size = 10),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 18),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
ggsave("figures/mh_PSV_richness_groups_barplot.jpg", plot=last_plot(), device="jpeg", path="figures/", scale = 1, units = "in", dpi = 600)

# PW wilcoxon tests of richness by type
mhpp_pw_sample_type <- pairwise.wilcox.test(sample_data(mhpp_PSV.ps)$Observed, sample_data(mhpp_PSV.ps)$Sample_Type, p.adjust.method = "BH")
mhpp_pw_sample_type$p.value # p-values in the matrix

# PW wilcoxon tests of richness by group
mhpp_pw_group <- pairwise.wilcox.test(sample_data(mhpp_PSV.ps)$Observed, sample_data(mhpp_PSV.ps)$Group, p.adjust.method = "BH")
mhpp_pw_group$p.value # p-values in the matrix


## Shannon's ####

# Plot richness and diversity
# Richness by sample type
ggplot(sample_data(mhpp_PSV.ps), mapping=aes(x= Sample_Type, y = Shannon, fill = Sample_Type, color = Sample_Type)) +
  theme_bw() + 
  labs(y = "Shannon") +
  #labs(x = "Sample Type") +
  #labs(title = "Unique Mh PSVs Identified") +
  #scale_fill_discrete(labels = c("Nasal", "Rope", "Waterbowl")) +
  #scale_color_discrete(labels = c("Nasal", "Rope", "Waterbowl")) +
  geom_boxplot(alpha=0.35, show.legend = F) +
  geom_jitter(width = .15) +
  geom_point(alpha = .05, shape = ".") +
  geom_smooth() +
  theme(legend.position = "right",
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        legend.title = element_text("Sample Type", size = 10),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 18),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
ggsave("figures/mh_PSV_shannon_barplot.jpg", plot=last_plot(), device="jpeg", path="figures/", scale = 1, units = "in", dpi = 600)

# Richness by group
ggplot(sample_data(mhpp_PSV.ps), mapping=aes(x= Group, y = Shannon, fill = Group, color = Group)) +
  theme_bw() + 
  labs(y = "Shannon") +
  #labs(x = "Sample Type") +
  #labs(title = "Unique Mh PSVs Identified") +
  #scale_fill_discrete(labels = c("Nasal", "Rope", "Waterbowl")) +
  #scale_color_discrete(labels = c("Nasal", "Rope", "Waterbowl")) +
  geom_boxplot(alpha=0.35, show.legend = F) +
  geom_jitter(width = .15) +
  geom_point(alpha = .05, shape = ".") +
  geom_smooth() +
  theme(legend.position = "right",
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        legend.title = element_text("Sample Type", size = 10),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 18),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
ggsave("figures/mh_PSV_shannon_groups_barplot.jpg", plot=last_plot(), device="jpeg", path="figures/", scale = 1, units = "in", dpi = 600)

# PW wilcoxon tests of richness by type
mhpp_pw_sample_type <- pairwise.wilcox.test(sample_data(mhpp_PSV.ps)$Shannon, sample_data(mhpp_PSV.ps)$Sample_Type, p.adjust.method = "BH")
mhpp_pw_sample_type$p.value # p-values in the matrix

# PW wilcoxon tests of richness by group
mhpp_pw_group <- pairwise.wilcox.test(sample_data(mhpp_PSV.ps)$Shannon, sample_data(mhpp_PSV.ps)$Group, p.adjust.method = "BH")
mhpp_pw_group$p.value # p-values in the matrix
