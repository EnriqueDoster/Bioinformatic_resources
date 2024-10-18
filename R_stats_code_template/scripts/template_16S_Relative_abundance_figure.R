# Most of this code is from Valeria
data_micro.ps <- prune_taxa(taxa_sums(data_micro.ps) > 0, data_micro.ps)
any(taxa_sums(data_micro.ps)==0) # QUADRUPLE CHECKING - nope good.

data_micro.ps.css <- phyloseq_transform_css(data_micro.ps, log = F)
data_micro.ps.css.df <- as(sample_data(data_micro.ps.css), "data.frame")


rel_abund.ps <- transform_sample_counts(data_micro.ps.css, function(x) {x/sum(x)}*100)

samples.family <- tax_glom(rel_abund.ps, taxrank = "Family", NArm = F) # 184 families

##I changed the merge_low_abundance function for the output to be "Other" instead of "zzzOther"
merge_low_abundance_ra <- function(data, threshold=1){
  otu.table <- as.data.frame(otu_table(data))
  otu.list <- row.names(otu.table[rowMeans(otu.table) < threshold,])
  merged <- merge_taxa(data, otu.list, 1)
  for (i in 1:dim(phyloseq::tax_table(merged))[1]){
    if (is.na(phyloseq::tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- paste0("Other")
      phyloseq::tax_table(merged)[i,1:7] <- paste0("Other")}
  }
  return(merged)
}

##Merging low abundance families
samples.family.filt <- merge_low_abundance_ra(samples.family, threshold = 0.6)
samples.family.filt #20 families over 0.5% RA
samples.family.filt.melt <- psmelt(samples.family.filt)   # melting to long format
samples.family.filt.melt <- samples.family.filt.melt %>%
  mutate(Family = factor(Family, levels = c(setdiff(Family, "Other"), "Other")))##Factoring the Family column so that "Other" is the last category

##Color palette
family.filt.palette <- distinctColorPalette(20) ##create a random color palette
names(family.filt.palette) <- unique(samples.family.filt.melt$Family) ##Assign colors from the palette to the taxa(family) names
family.filt.palette$'Other' <- "grey95" ##Making sure "Other is "grey95"

##Another function (top_taxa_legend)
##Selecting the top most abundant taxa to include in the RA legend
top_taxa_legend <- function(data, taxlevel = "Family",  n = 10) { #Taxlevel takes a string ("Family", "Genus", etc), and data[[level]] gets the corresponding column from data (data should be a melted df)
  # If you've already "factored" the taxlevel column, this takes that column and turns it back into a character vector (I did this so later on I'd get actual taxa names, not the factor numbers)
  taxlevel_column <- as.character(data[[taxlevel]])
  #Then, aggregate abundances by tax level (getting the mean RA across samples)
  taxlevel_abundance <- aggregate(Abundance ~ taxlevel_column, data = data, mean) #Here, the result is a df with two columns, first one are taxa names, second one is the average abundance
  # I'll rename the first column of taxlevel_abundance. Instead of it being "data[[taxlevel]]", this changes it to the "taxlevel" provided.
  colnames(taxlevel_abundance)[1] <- taxlevel
  taxlevel_abundance <- taxlevel_abundance[order(-taxlevel_abundance$Abundance), ]# This orders taxa by abundance
  top_taxa <- head(taxlevel_abundance[[taxlevel]], n)# This selects the top (n) taxa and returns only the names
  # Now, to handle Other (just making sure it is at the bottom of the list, but may not be necessary)
  if ("Other" %in% top_taxa) {
    # If "Other" is in top taxa, this removes it
    top_taxa <- top_taxa[top_taxa != "Other"]
    # Then, it places "Other " at the bottom of the character list
    top_taxa <- c(top_taxa, "Other")
  } else {
    # If "Other " is not in top taxa, it just adds "Other " at the bottom of the list
    top_taxa <- c(top_taxa, "Other")
  }
  return(top_taxa)
}

##Apply the function to obtain top families
top_families <- top_taxa_legend(samples.family.filt.melt, taxlevel = "Family", n = 10)
top_families

##If you didn't "factor" before using the top_taxa_legend function, you could do it here before you plot
#samples.family.filt.melt <- samples.family.filt.melt %>%
#mutate(Family = factor(Family, levels = c(setdiff(Family, "Other"), "Other")))

#Plotting- make sure to add "breaks"
dendroRA.family.plot <- ggplot(samples.family.filt.melt, aes(x=Sample, y= Abundance, fill = Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  #scale_x_discrete(limits = dendro_bray_sample_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values = family.filt.palette, breaks = top_families) +
  theme(#legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank())
dendroRA.family.plot
