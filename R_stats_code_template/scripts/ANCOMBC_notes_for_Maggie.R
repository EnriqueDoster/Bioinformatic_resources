# HFT ANCOMBC
set.seed(123)
library(ANCOMBC)
# library(DAtest) # optional, for advanced pruning
library(microbiome)
library(DT)

## The types of tests that ANCOMBC performs
# 1) Global test = "any difference between groups"
# 2) Pairwise comparisons = "difference between specific groups"
# 3) Trend test = "tests for difference in taxa across numeric variable (eg. age, days, dose)"
# 4) Sensitivity analysis = basically tests for FDR for each taxa
# 5) Structural zero determination = look it up


###
##### Start ANCOMBC - Fecal samples, comparing by Combined_order #####
###

# Subset the fecal samples from raw data in a phyloseq object
# Here you'll subset the individuals 
fecal_AMR_data.ps <- subset_samples(AMR_data.ps, sample_type == "Fecal")

# Easy prune of taxa with 0 counts
fecal_AMR_data.ps = prune_taxa(taxa_sums(fecal_AMR_data.ps) > 0, fecal_AMR_data.ps) # 2242 taxa, down from 2912 in all samples

# Trim taxa using DAtest
# filtered_fecal_data_group <- preDA(fecal_AMR_data.ps, min.samples=2, min.reads = 1) # 10% of samples (77/154)
#


# Specify the order of variables in a factor
# your first variable here will be used as the "reference" in the model
sample_data(fecal_AMR_data.ps)$Combined_order <- factor(sample_data(fecal_AMR_data.ps)$Combined_order, levels = c("Start_CONV","Midpoint_CONV",
                                                                                                                  "End_CONV","Start_RWA","Midpoint_RWA",
                                                                                                                  "End_RWA","Final_washout"))
# This is optional, but need to remove 1 factor level to make a square matrix possible for the trend test
#fecal_AMR_data.ps <- subset_samples(AMR_data.ps, Combined_order != "Final_washout")

# Change variable in metadata to factor
# You might do this with age, or treat it numerically and try the "trend test"
sample_data(fecal_AMR_data.ps)$participant_id <- as.factor(sample_data(fecal_AMR_data.ps)$participant_id)

# ANCOMBC-specific issue (might not be required for microbiome data):
# Change label for most specific annotation level, in this case AMR gene group, to "Species"
colnames(phyloseq::tax_table(fecal_AMR_data.ps))[4] <- "Species"

# ANCOMBC model code
# Adjust the "fix_formula" to match your study design".
# Change the "group" variable to whatever you want to run pairwise comparisons on.
# We can play around with the "rand_formula" later for repeated measures.
ancom_output_1 = ancombc2(data = fecal_AMR_data.ps, assay_name = "counts", tax_level = "Species",
                          fix_formula = "Combined_order + participant_id + blinded_treatment", rand_formula = NULL,
                          p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                          group = "Combined_order", struc_zero = TRUE, neg_lb = TRUE,
                          alpha = 0.05, n_cl = 2, verbose = TRUE,
                          global = TRUE, pairwise = TRUE,
                          dunnet = TRUE, trend = FALSE,
                          iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE),
                          em_control = list(tol = 1e-5, max_iter = 100),
                          lme_control = lme4::lmerControl(),
                          mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                          trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                      nrow = 2,
                                                                      byrow = TRUE),matrix(c(-1, 0, 1, -1),
                                                                                           nrow = 2,
                                                                                           byrow = TRUE)),
                                               node = list(2, 2),
                                               solver = "ECOS",
                                               B = 100))



# Extract results from pairwise comparisons
res_pair = ancom_output_1$res_pair

#### Calculate sensitivity scores for pool_timepoint ####
## Rename pairwise comparisons
df_Combined_order = res_pair %>%
  dplyr::select(taxon, contains("Combined_order")) 

# Extract the pairwise comparison results
df_fig_Combined_order = df_Combined_order %>%
  filter(diff_Combined_orderMidpoint_CONV == 1 | diff_Combined_orderEnd_CONV_Combined_orderMidpoint_CONV == 1 
         | diff_Combined_orderStart_RWA_Combined_orderEnd_CONV  == 1 | diff_Combined_orderMidpoint_RWA_Combined_orderStart_RWA  == 1 
         | diff_Combined_orderEnd_RWA_Combined_orderMidpoint_RWA == 1 | diff_Combined_orderFinal_washout_Combined_orderEnd_RWA == 1 ) %>%
  mutate(lfc_Midpoint_CONV = ifelse(diff_Combined_orderMidpoint_CONV == 1,  # This compares the 2nd time point to Start_CONV
                                    lfc_Combined_orderMidpoint_CONV, 0),
         lfc_End_CONV = ifelse(diff_Combined_orderEnd_CONV_Combined_orderMidpoint_CONV== 1, 
                               lfc_Combined_orderEnd_CONV_Combined_orderMidpoint_CONV, 0),
         lfc_Start_RWA_CONV = ifelse(diff_Combined_orderStart_RWA == 1, 
                                     lfc_Combined_orderStart_RWA, 0),
         lfc_Midpoint_RWA = ifelse(diff_Combined_orderMidpoint_RWA_Combined_orderStart_RWA== 1, 
                                   lfc_Combined_orderMidpoint_RWA_Combined_orderStart_RWA, 0),
         lfc_End_RWA = ifelse(diff_Combined_orderEnd_RWA_Combined_orderMidpoint_RWA== 1, 
                              lfc_Combined_orderEnd_RWA_Combined_orderMidpoint_RWA, 0)
  ) %>%
  transmute(taxon, 
            `Midpoint_CONV vs. Start_CONV` = round(lfc_Midpoint_CONV, 2),
            `End_CONV vs. Midpoint_CONV` = round(lfc_End_CONV, 2),
            `Start_RWA vs. Start_CONV` = round(lfc_Start_RWA_CONV, 2),
            `Midpoint_RWA vs. Start_RWA` = round(lfc_Midpoint_RWA, 2),
            `End_RWA vs. Midpoint_RWA` = round(lfc_End_RWA, 2)) %>%
  pivot_longer(cols = `Midpoint_CONV vs. Start_CONV`:`End_CONV vs. Midpoint_CONV`:`Start_RWA vs. Start_CONV`:`Midpoint_RWA vs. Start_RWA`:`End_RWA vs. Midpoint_RWA`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon)

# Extract sensitivity analysis for the pairwise comparison results
df_fig_Combined_order_color = df_Combined_order %>%
  filter(diff_Combined_orderMidpoint_CONV == 1 | diff_Combined_orderEnd_CONV_Combined_orderMidpoint_CONV == 1 
         | diff_Combined_orderStart_RWA_Combined_orderEnd_CONV  == 1 | diff_Combined_orderMidpoint_RWA_Combined_orderStart_RWA  == 1 
         | diff_Combined_orderEnd_RWA_Combined_orderMidpoint_RWA == 1 | diff_Combined_orderFinal_washout_Combined_orderEnd_RWA == 1 ) %>%
  mutate(lfc_Midpoint_CONV = ifelse(diff_Combined_orderMidpoint_CONV == 1 & passed_ss_Combined_orderMidpoint_CONV == 1,  # This compares the 2nd time point to Start_CONV
                                    "aquamarine3", "black"),
         lfc_End_CONV = ifelse(diff_Combined_orderEnd_CONV_Combined_orderMidpoint_CONV== 1 & passed_ss_Combined_orderEnd_CONV_Combined_orderMidpoint_CONV == 1, 
                               "aquamarine3", "black"),
         lfc_Start_RWA_CONV = ifelse(diff_Combined_orderStart_RWA == 1 & passed_ss_Combined_orderStart_RWA == 1, 
                                     "aquamarine3", "black"),
         lfc_Midpoint_RWA = ifelse(diff_Combined_orderMidpoint_RWA_Combined_orderStart_RWA== 1 & passed_ss_Combined_orderMidpoint_RWA_Combined_orderStart_RWA == 1, 
                                   "aquamarine3", "black"),
         lfc_End_RWA = ifelse(diff_Combined_orderEnd_RWA_Combined_orderMidpoint_RWA== 1 & passed_ss_Combined_orderEnd_RWA_Combined_orderMidpoint_RWA == 1, 
                              "aquamarine3", "black")
  ) %>%
  transmute(taxon, 
            `Midpoint_CONV vs. Start_CONV` = lfc_Midpoint_CONV,
            `End_CONV vs. Midpoint_CONV` = lfc_End_CONV,
            `Start_RWA vs. Start_CONV` = lfc_Start_RWA_CONV,
            `Midpoint_RWA vs. Start_RWA` = lfc_Midpoint_RWA,
            `End_RWA vs. Midpoint_RWA` = lfc_End_RWA) %>%
  pivot_longer(cols = `Midpoint_CONV vs. Start_CONV`:`End_CONV vs. Midpoint_CONV`:`Start_RWA vs. Start_CONV`:`Midpoint_RWA vs. Start_RWA`:`End_RWA vs. Midpoint_RWA`, 
               names_to = "group", values_to = "color") %>%
  arrange(taxon)


df_fig_pair = df_fig_Combined_order %>%
  dplyr::left_join(df_fig_Combined_order_color, by = c("taxon", "group"))



# Change order of X-value
df_fig_pair
df_fig_pair$group <- factor(df_fig_pair$group, levels = c("Start_RWA vs. Start_CONV","Midpoint_CONV vs. Start_CONV",
                                                                              "End_CONV vs. Midpoint_CONV","Midpoint_RWA vs. Start_RWA","End_RWA vs. Midpoint_RWA"))


# Plot heatmap
lo = floor(min(df_fig_pair$value))
up = ceiling(max(df_fig_pair$value))
mid = (lo + up)/2


fig_pair = df_fig_pair %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value, color = color), size = 4) +
  scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = "Log fold changes by combined order") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_pair



## Testing new code #####



##
### Sensitivity scores #####
##

tab_sens = ancom_output_1$ss_tab

tab_sens %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(tab_sens)[-1], digits = 2)

