############    SNV FIGS ################
#Scripts for figs only, 9/5/2023 MM


#Alpha Diversity
#Observed ARGs by sample type 

ggplot(alpha_div_meta, aes(x= sample_type, y= Observed, fill = sample_type, colour = sample_type)) +
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


#Shannons Diversity
#By sample type

ggplot(alpha_div_meta, aes(x= sample_type, y= Shannon, fill = sample_type, colour = sample_type)) +
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

#Observed ARGs
#Fecal, by timepoint

ggplot(alpha_div_fecal, aes(x= Combined_order, y= Observed, fill = Combined_order, colour = Combined_order)) +
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


#Shannons
#Fecal by timepoint


ggplot(alpha_div_fecal, aes(x= Combined_order, y= Shannon, fill = Combined_order, colour = Combined_order)) +
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


#NMDS
#By sample type

plt_ord_by_sampletype <- plot_ordination(ps_class.css, data.ord, color = "sample_type") +
  theme_bw() +
  labs(title ="Sample type") +
  stat_ellipse(aes(fill= sample_type), geom="polygon", alpha = 0.25) +
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
plt_ord_by_sampletype 


#NMDS
#By sample group, trial (conv beef, conv fecal, rwa beef, rwa fecal)

plt_ord_ps_class <- plot_ordination(ps_class.css, ps_class.ord, color = "Trial", shape = "Trial") +
  theme_bw() +
  labs(title ="Sample group") +
  stat_ellipse(aes(fill= Trial), geom="polygon", alpha = 0.25) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", size = 1.0),
        strip.text = element_text(size =22, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 24),
        axis.ticks = element_line(colour = "black", size = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
plt_ord_ps_class 


#Relative Abundance


plt_rel_class_ps_class <- ggplot(joined_ra_class_ps_class_melt, 
                                 aes(x = x_center, y = Abundance, fill = class, 
                                     height = height, width = width)) + 
  coord_flip() +
  geom_bar(stat = "identity", colour = "black") +
  #scale_fill_manual(values = col_vector) + #use this if not using color palette from ggthemes of ggsci
  scale_fill_tableau(palette = "Tableau 20")+
  scale_x_continuous(breaks = sample_pos_table_ps_AMR$x_center, 
                     labels = sample_pos_table_ps_AMR$sample, 
                     expand = c(0, 0)) + 
  # For the y axis, alternatively set the labels as: gene_position_table$gene
  #scale_y_continuous(breaks = gene_pos_table[, "y_center"], labels = rep("", nrow(gene_pos_table)),limits = gene_axis_limits, expand = c(0, 0)) + 
  labs(x = "", y = "Relative abundance") +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
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
plt_rel_class_ps_class

#Dendrogram

plt_dendr_class_ps_class <- ggplot(segment_data_ps_AMR) + 
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend), linewidth =.75, lineend = "round", linejoin = "round") +
  geom_point(data = gene_pos_table_ps_AMR, aes(x,y_center, colour = trt, fill = trt, shape = trt),
             size = 8, stroke = 0.75, position = position_nudge(x = -0.03)) +
  scale_y_continuous(breaks = gene_pos_table_ps_AMR$y_center, 
                     labels = gene_pos_table_ps_AMR$gene, 
                     limits = gene_axis_limits_ps_class, 
                     expand = c(0, 0)) + 
  labs(x = "Ward's Distance", y = "", colour = "", size = "", title = "Sample type") +
  scale_x_reverse() + 
  scale_color_manual(values=c( "yellow", 'green',"#ff0000", "#cc0000", "#990000", "#0000ff", "#0000cc", "#000099", "#808080")) +
  scale_shape_manual(values =c(15,15,15,15,15,15,15,15,15)) +
  theme_bw() + 
  theme(legend.position = "right",
        panel.border = element_blank(),
        plot.title = element_text(size = 30),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(size = 0.75),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, colour = "black"))
plt_dendr_class_ps_class

#Combine

class_ps_class_plots <- plot_grid(plt_dendr_class_ps_class, plt_rel_class_ps_class, align = 'h', rel_widths = c(0.5, 1.5))
class_ps_class_plots

#NMDS
#Sample type, Combined order (Conv beef, end conv, end rwa, mid conv, mid rwa, rwa beef)

plt_ord_by_sampletype <- plot_ordination(ps_group.css, ps_group.ord, color = "Combined_order") +
  theme_bw() +
  labs(title ="Sample type") +
  stat_ellipse(aes(fill= Combined_order), geom="polygon", alpha = 0.25) +
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
plt_ord_by_sampletype 



#NMDS
#ALL Samples
# plot the ordination with 95% confidence ellipses coloured by groups
# (conv beef, end conv, end rwa, final washout, mid conv, mid rwa, rwa beef, start conv, start rwa)


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



#NMDS
#FECAL ONLY
# plot the ordination with 95% confidence ellipses coloured by groups
# (conv beef, end conv, end rwa, final washout, mid conv, mid rwa, rwa beef, start conv, start rwa)



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







