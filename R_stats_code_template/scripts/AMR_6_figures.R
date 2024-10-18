# Figures
# make sure this object is loaded
data_noSNP.css
rel_abund_data_noSNP.css <- transform_sample_counts(data_noSNP.css, function(x) {x/sum(x)}*100)


# One way is to aggregate to group level first
group_data_noSNP.css <- tax_glom(data_noSNP.css, taxrank = "Group")
# Remove taxa with no counts
group_data_noSNP.css <- prune_taxa(taxa_sums(group_data_noSNP.css) > 0, group_data_noSNP.css)
# Optional check
any(taxa_sums(group_data_noSNP.css)==0) # QUADRUPLE CHECKING - nope good.


#
##
### Plot tetracycline counts, by AMR group ####
##
#
tet_class.css <- subset_taxa(data_noSNP.css, Class == "Tetracyclines")
plot_bar(tet_class.css, x = "sample_type", fill = "Group")



#
##
### Plot a subset of features,  ####
##
#
tet_class_melted.css <- subset_taxa(data_noSNP.css, Class == "Tetracyclines") %>%
  psmelt() # you should be familiar with the melt function and how that data looks. Its the most common format for plotting figures
ggplot(tet_class_melted.css, aes(x= sample_type, y= Abundance, fill = Group)) +
  theme_bw() +
  #facet_grid(~sample_type) +
  geom_bar(stat = "summary", colour = "black") 

