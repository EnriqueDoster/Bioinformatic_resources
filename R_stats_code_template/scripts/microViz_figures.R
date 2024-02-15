# MicroViz figures

# Here are all the packages you'll need.
# Follow install options https://david-barnett.github.io/microViz/index.html
options(width = 100)
library(microViz)
library(phyloseq)
library(ggplot2)
library(DT)
library(ggraph)
library(patchwork) # for arranging groups of plots
library(metagMisc)
library(microbiome)
#knitr::opts_chunk$set(fig.height = 6, fig.width = 9)

  
# Transform counts to CSS normalized (we also do this for ordination, relative abundance plots, and testing with ZIG)
data_noSNP.css <- phyloseq_transform_css(data_noSNP, log = F)


# SNV data
SNV_data_noSNP.css <- phyloseq_transform_css(SNV_data_noSNP, log = F)


ord_explore(SNV_data_noSNP)


#
### Example of barplot ####
#

plot_list <- data_noSNP.css %>%
  ps_filter(sample_type == "Beef sample") %>%             # Filter out just Beef samples from the "sample_type" column
  tax_fix() %>%
  comp_barplot( group_by = "specific_sample_type",        # What variable you want to group them by
                tax_level = "Class",                      # What level to view results in
                n_taxa = 15,                              # How many taxa to show, the rest will be grouped into an "Other" category
                sample_order = "bray",                    # Use bray-cutis clustering to sort samples within each group
                label = "treatment")                      # Change the label to the treatment variable

# Plot them side by side with the patchwork package.
patch <- patchwork::wrap_plots(plot_list, nrow = 1, guides = "collect")
patch & coord_flip() # make all plots horizontal (note: use & instead of +)


#
### Grouping - at the group level for Tetracycline resistance ####
#

# Just subset out the Tetracycline group
tet_class.css <- subset_taxa(data_noSNP.css, Class == "Tetracyclines")

tet_plot_list <- tet_class.css %>%
  ps_filter(sample_type == "Beef sample") %>%
  tax_fix() %>%
  comp_barplot( group_by = "specific_sample_type",
                tax_level = "Group", # At which level to show results
                n_taxa = 15,
                sample_order = "bray", # Cluster samples using bray-curtis distance
                label = "treatment") # Change label seen for each sample

# Plot them side by side with the patchwork package.
tet_patch <- patchwork::wrap_plots(tet_plot_list, nrow = 1, guides = "collect")
tet_patch & coord_flip() # make all plots horizontal (note: use & instead of +)


#
### Grouping - at the group level for fecal by participant ID ####
#

# participant_id should be coded as a factor
sample_data(data_noSNP.css)$participant_id <- as.factor(sample_data(data_noSNP.css)$participant_id)

# Select the order of factors for "Combined_order"
sample_data(data_noSNP.css)$Combined_order <- factor(sample_data(data_noSNP.css)$Combined_order, levels = c("Start_CONV","Midpoint_CONV","End_CONV",
                                                                                                            "Start_RWA","Midpoint_RWA","End_RWA","Final_washout"))

data_noSNP.css %>%
  ps_filter(sample_type == "Fecal" & Combined_order != "Final_washout") %>%
  ps_arrange(participant_id) %>%
  #ps_mutate(
 #   participant_id = paste("Participant", participant_id), # better labels
 #   participant_id = factor(participant_id, rev(unique(participant_id))) # fix plot order
#  ) %>%
  comp_barplot(
    tax_level = "Class", bar_width = 0.7, sample_order = "asis",n_taxa = 10, # don't sort
    x = "participant_id" # x argument is available since microViz 0.9.7
  ) +
  facet_wrap(
    facets = vars(Combined_order), #labeller = as_labeller(~ paste("Age", ., "days")
    scales = "fixed", ncol = 6
  ) +
  coord_flip() +
  labs(x = "Participant ID", y = "Relative abundance") +
  scale_y_continuous(expand = expansion(add = c(0, 0.05))) + # axis starts exactly at 0
  theme_bw() + # slightly clearer axes for facets
  theme(panel.spacing.x = unit(6, "mm"))


###
##### Now trying interactive ordination plot #####
###
ord_explore(data_noSNP.css)
