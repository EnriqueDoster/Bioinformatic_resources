#############################################################################################
##############################         ALPHA DIVERSITY         ##############################
#############################################################################################
#############################################################################################

#Remove Rinsate Pools
#Is this even necessary? 
#SNV_data_noSNP<- subset_samples(SNV_data_noSNP, sample_type != "Rinsate Pool")
#SNV_data_noSNP # 216 samples  YAY! 216 samples does not contain rinsate pools

#I began this script using snv_dat.ps object, but it wasn't working because it did not contain the Combined order we made
#in the previous script. Changed it to SNV_data_noSNP and it is working. 


SNV_alpha_div <- estimate_richness(SNV_data_noSNP, measures = c("Observed","Shannon","Simpson","InvSimpson"))
SNV_alpha_div
SNV_alpha_div.df <- as(sample_data(SNV_data_noSNP), "data.frame")
SNV_alpha_div_meta <- cbind(alpha_div, alpha_div.df)

#I spent too much time trying to figure out how to remove the RP from the ggplot below. Moved these two lines to above and it all works. Duh. 
#Remove Rinsate Pools
#SNV_data_noSNP<- subset_samples(SNV_data_noSNP, sample_type != "Rinsate Pool")
#SNV_data_noSNP # 216 samples  YAY! 216 samples does not contain rinsate pools


ggplot(SNV_alpha_div_meta, aes(x= sample_type, y= Observed, fill = sample_type, colour = sample_type)) +
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
sample_type.richness.pw <- pairwise.wilcox.test(SNV_alpha_div_meta$Observed, SNV_alpha_div_meta$sample_type, p.adjust.method = "BH")
sample_type.richness.pw # p-values in the matrix


ggplot(SNV_alpha_div_meta, aes(x= sample_type, y= Shannon, fill = sample_type, colour = sample_type)) +
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
sample_type.shannon.pw <- pairwise.wilcox.test(SNV_alpha_div_meta$Shannon, SNV_alpha_div_meta$sample_type, p.adjust.method = "BH")
sample_type.shannon.pw # p-values in the matrix


##
#### Now, only fecal samples #####
##

# Extract fecal samples
SNV_alpha_div_fecal <- SNV_alpha_div_meta[which(SNV_alpha_div_meta$sample_type=="Fecal"),]

# Change order of "Combined_order"
SNV_alpha_div_fecal$Combined_order <- factor(SNV_alpha_div_fecal$Combined_order, levels = c("Start_CONV","Midpoint_CONV","End_CONV", "Start_RWA",
                                                                                    "Midpoint_RWA","End_RWA","Final_washout"))


ggplot(SNV_alpha_div_fecal, aes(x= Combined_order, y= Observed, fill = Combined_order, colour = Combined_order)) +
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
Combined_order.richness.pw <- pairwise.wilcox.test(SNV_alpha_div_fecal$Observed, SNV_alpha_div_fecal$Combined_order, p.adjust.method = "BH")
Combined_order.richness.pw # p-values in the matrix


## Shannon's #
ggplot(SNV_alpha_div_fecal, aes(x= Combined_order, y= Shannon, fill = Combined_order, colour = Combined_order)) +
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
Combined_order.shannon.pw <- pairwise.wilcox.test(SNV_alpha_div_fecal$Shannon, SNV_alpha_div_fecal$Combined_order, p.adjust.method = "BH")
Combined_order.shannon.pw # p-values in the matrix

