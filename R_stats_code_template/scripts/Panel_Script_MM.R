#Script for Panels in Cowplot#
#Richness and Diversity for AMR, SNV, and Microbiome

#Change theme to none
#Change axis.title.y=element_text(size=12)
#Change labs name


#AMR Alpha Sample Type
p1 <- ggplot(AMR_alpha_div_meta, aes(x= sample_type, y= Observed, fill = sample_type, colour = sample_type)) +
  theme_bw() + 
  labs(y= "") + 
  geom_violin(alpha=0.4) +
  geom_point() +
    theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 28),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
p1
#p1+geom_text(data = labs, aes(x = sample_type, y = -1, label = label), vjust = -1)

#SNV Alpha Sample Type
p2 <- ggplot(SNV_alpha_div_meta, aes(x= sample_type, y= Observed, fill = sample_type, colour = sample_type)) +
  theme_bw() + 
  labs(y= "") +
  geom_violin(alpha=0.4) +
  geom_point() +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 28),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
p2




#Micro Alpha Sample Type
p3 <- ggplot(micro.alpha_div_meta, aes(x= sample_type, y= Observed, fill = sample_type, colour = sample_type)) +
  theme_bw() + 
  labs(y= "") +
  geom_violin(alpha=0.4) +
  geom_point() +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 28),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
p3

#AMR Shannons
p1_1 <- ggplot(AMR_alpha_div_meta, aes(x= sample_type, y= Shannon, fill = sample_type, colour = sample_type)) +
  theme_bw() + 
  labs(y= "") +
  geom_violin(alpha=0.4) +
  geom_point() +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 28),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
p1_1


p2_2 <- ggplot(SNV_alpha_div_meta, aes(x= sample_type, y= Shannon, fill = sample_type, colour = sample_type)) +
  theme_bw() + 
  labs(y= "") +
  geom_violin(alpha=0.4) +
  geom_point() +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 28),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
p2_2


p3_3 <- ggplot(micro.alpha_div_meta, aes(x= sample_type, y= Shannon, fill = sample_type, colour = sample_type)) +
  theme_bw() + 
  labs(y= "") +
  geom_violin(alpha=0.4) +
  geom_point() +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 28),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
p3_3


#Plot final grid
final_plot <- plot_grid(
  plot_grid(p1, p1_1, labels = c("a", "b"), ncol = 2),
  plot_grid(p2, p2_2, labels = c("c", "d"), ncol = 2),
  plot_grid(p3, p3_3, labels = c("e","f"), ncol = 2),
  nrow = 3,
  align = "hv")
final_plot




