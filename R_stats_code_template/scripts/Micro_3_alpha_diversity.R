#############################################################################################
##############################         ALPHA DIVERSITY         ##############################
#############################################################################################
#############################################################################################

micro.alpha_div <- estimate_richness(data_micro.ps, measures = c("Observed","Shannon","Simpson","InvSimpson"))
micro.alpha_div

micro.alpha_div.df <- as(sample_data(data_micro.ps), "data.frame")
micro.alpha_div_meta <- cbind(micro.alpha_div, micro.alpha_div.df)


ggplot(micro.alpha_div_meta, aes(x= sample_type, y= Observed, fill = sample_type, colour = sample_type)) +
  theme_bw() + 
  labs(y= "Observed ASVs") +
  geom_boxplot(alpha=0.4) +
  geom_point() +
  theme(legend.position = "right",
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


#### pairwise Wilxocon rank-sum with Benjamini-Hochberg correction for mult. comps.
micro.sample_type.richness.pw <- pairwise.wilcox.test(micro.alpha_div_meta$Observed, micro.alpha_div_meta$sample_type, p.adjust.method = "BH")
micro.sample_type.richness.pw # p-values in the matrix


ggplot(micro.alpha_div_meta, aes(x= sample_type, y= Shannon, fill = sample_type, colour = sample_type)) +
  theme_bw() + 
  labs(y= "Shannon's diversity") +
  geom_boxplot(alpha=0.4) +
  geom_point() +
  theme(legend.position = "right",
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

#### pairwise Wilxocon rank-sum with Benjamini-Hochberg correction for mult. comps.
micro.sample_type.shannon.pw <- pairwise.wilcox.test(micro.alpha_div_meta$Shannon, micro.alpha_div_meta$sample_type, p.adjust.method = "BH")
micro.sample_type.shannon.pw # p-values in the matrix


##
#### Now, only fecal samples #####
##

# Extract fecal samples
micro.alpha_div_fecal <- micro.alpha_div_meta[which(micro.alpha_div_meta$sample_type=="Fecal"),]

# Change ordre of "Combined_order"
micro.alpha_div_fecal$Combined_order <- factor(micro.alpha_div_fecal$Combined_order, levels = c("Start_CONV","Midpoint_CONV","End_CONV", "Start_RWA",
                                                                                    "Midpoint_RWA","End_RWA","Final_washout"))


ggplot(micro.alpha_div_fecal, aes(x= Combined_order, y= Observed, fill = Combined_order, colour = Combined_order)) +
  theme_bw() + 
  labs(y= "Observed ARGs") +
  geom_boxplot(alpha=0.4) +
  geom_point() +
  theme(legend.position = "right",
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

#### pairwise Wilxocon rank-sum with Benjamini-Hochberg correction for mult. comps.
micro.Combined_order.richness.pw <- pairwise.wilcox.test(micro.alpha_div_fecal$Observed, micro.alpha_div_fecal$Combined_order, p.adjust.method = "BH")
micro.Combined_order.richness.pw # p-values in the matrix


## Shannon's #
ggplot(micro.alpha_div_fecal, aes(x= Combined_order, y= Shannon, fill = Combined_order, colour = Combined_order)) +
  theme_bw() + 
  labs(y= "Shannon's diversity") +
  geom_boxplot(alpha=0.4) +
  geom_point() +
  theme(legend.position = "right",
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

#### pairwise Wilxocon rank-sum with Benjamini-Hochberg correction for mult. comps.
micro.Combined_order.shannon.pw <- pairwise.wilcox.test(micro.alpha_div_fecal$Shannon, micro.alpha_div_fecal$Combined_order, p.adjust.method = "BH")
micro.Combined_order.shannon.pw # p-values in the matrix

