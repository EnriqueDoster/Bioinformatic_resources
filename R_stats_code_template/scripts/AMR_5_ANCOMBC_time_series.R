# ANCOMBC time series analysis
library(ANCOMBC)
library(dplyr)

###
##### Start ANCOMBC - Fecal samples, comparing by Combined_order #####
###

CONV_fecal_data.ps <- subset_samples(AMR_data.ps, sample_type == "Fecal" & treatment_order == "AB" & Combined_order != "Final_washout")

colnames(phyloseq::tax_table(CONV_fecal_data.ps))[4] <- "Species"

sample_data(CONV_fecal_data.ps)$Week <- as.integer(sample_data(CONV_fecal_data.ps)$Week)

CONV_ANCOM_feces_output = ancombc2(data = CONV_fecal_data.ps, 
                                   assay_name = "counts", 
                                   tax_level = "Species",
                                   fix_formula = "participant_id + Week",
                                   p_adj_method = "holm",  
                                   struc_zero = FALSE, 
                                   neg_lb = FALSE,
                                   verbose = TRUE,
                                   global = FALSE, 
                                   pairwise = FALSE, 
                                   n_cl = 3)

# Extract the results from the ANCOMBC2 output
res_prim = CONV_ANCOM_feces_output$res

# Select the relevant columns for the "Week" variable
df_week = res_prim %>%
  dplyr::select(taxon, ends_with("Week"))

# Filter, arrange, and mutate the data for the figure
df_fig_week = df_week %>%
  dplyr::filter(diff_Week == 1, passed_ss_Week == 1) %>% 
  dplyr::arrange(desc(lfc_Week)) %>%
  dplyr::mutate(direct = ifelse(lfc_Week > 0, "Positive LFC", "Negative LFC"))

# Convert the "taxon" and "direct" columns to factors
df_fig_week$taxon = factor(df_fig_week$taxon, levels = df_fig_week$taxon)
df_fig_week$direct = factor(df_fig_week$direct, levels = c("Positive LFC", "Negative LFC"))

# Create the figure
fig_week = df_fig_week %>%
  ggplot(aes(x = taxon, y = lfc_Week, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_Week - se_Week, ymax = lfc_Week + se_Week), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes as one unit increase of Week") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))

# Print the figure
print(fig_week)

