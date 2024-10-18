feces.css.df <- as(feces.css@sam_data,"data.frame") # make dataframe from metadata
feces.gunifrac.dist <- gunifrac(feces.css) # calculate generalized UniFrac distances
feces.gunifrac.ord <- metaMDS(comm = feces.gunifrac.dist, try = 20, trymax = 500, autotransform = F) # ordinate distances

## create dataframe with the points for individuals and for the centroids
feces.gunifrac.plot <- ordiplot(feces.gunifrac.ord$points)
feces.gunifrac.scrs <- scores(feces.gunifrac.plot, display = "sites")
feces.gunifrac.scrs <- cbind(as.data.frame(feces.gunifrac.scrs), treatment= feces.css.df$treatment)
feces.gunifrac.cent <- aggregate(cbind(MDS1,MDS2) ~ treatment, data = feces.gunifrac.scrs, FUN = mean) # calculate the centroids for variable of interest
feces.gunifrac.segs <- merge(feces.gunifrac.scrs, setNames(feces.gunifrac.cent, c("treatment","cMDS1","cMDS2")), by = 'treatment', sort = F)
feces.gunifrac.segs

## plot it
ggplot(feces.gunifrac.segs, aes(fill = treatment, colour = treatment)) + theme_bw() +
  labs(x="NMDS1", y="NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(aes(x=MDS1,y=MDS2, colour= treatment), alpha = 0.5, shape = 18, size = 5) +
  stat_ellipse(geom = "polygon", aes(x=MDS1,y=MDS2), alpha = c(0.05), lty=2, level = 0.95, linewidth = 0.2) +
  geom_segment(aes(x=MDS1,y=MDS2, xend=cMDS1, yend=cMDS2), alpha = 0.1, linewidth = 1) +
  geom_point(aes(x=cMDS1, y=cMDS2), size = 18, shape = 18) +
  geom_text(aes(x=cMDS1, y=cMDS2, label = treatment), colour = "white", size = 5) +
  scale_fill_manual(values = treatment_palette) +
  scale_colour_manual(values = treatment_palette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))