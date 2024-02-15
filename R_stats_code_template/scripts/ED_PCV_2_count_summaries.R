# PSV 

PSV_rel_abund <- transform_sample_counts(mhpp_PSV.ps, function(x) {x/sum(x)}*100)

PSV_rel_abund_melt <- psmelt(PSV_rel_abund)

# calculate percentage across all counts
all_sample_PSV_by_abund <- PSV_rel_abund_melt %>%
  dplyr::group_by(OTU) %>%
  dplyr::summarize(median_PSV = median(Abundance)) %>%
  arrange(-median_PSV)

all_sample_PSV_by_abund

# 
# # Example for doing at other taxa level, this doesn't work for PSVs
# ra_class <- tax_glom(PSV_rel_abund, taxrank = "Class")
# ra_class 
# 
# ra_class_melt <- psmelt(ra_class)
# 
# # calculate percentage across all counts
# all_sample_PSV_class_by_abund <- ra_class_melt %>%
#   dplyr::group_by(class) %>%
#   dplyr::summarize(median_PSV = median(Abundance)) %>%
#   arrange(-median_PSV)
# 
# all_sample_PSV_class_by_abund

