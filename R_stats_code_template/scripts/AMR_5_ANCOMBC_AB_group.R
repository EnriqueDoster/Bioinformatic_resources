# HFT ANCOMBC
set.seed(123)
library(ANCOMBC)
library(DAtest) # optional, for advanced pruning
library(microbiome)

###
##### Start ANCOMBC - Fecal samples, comparing by Combined_order #####
###

# Subset the fecal samples from raw data
fecal_AMR_data.ps <- subset_samples(AMR_data.ps, sample_type == "Fecal")

AB_fecal_AMR_data.ps <- subset_samples(fecal_AMR_data.ps, treatment_order == "AB")

# Easy prune of taxa with 0 counts
AB_fecal_AMR_data.ps = prune_taxa(taxa_sums(AB_fecal_AMR_data.ps) > 0, AB_fecal_AMR_data.ps) # 2242 taxa, down from 2912 in all samples

# Trim taxa using DAtest
# filtered_fecal_data_group <- preDA(AB_fecal_AMR_data.ps, min.samples=2, min.reads = 1) # 10% of samples (77/154)
#


# Specify the order of variables in a factor
sample_data(AB_fecal_AMR_data.ps)$Combined_order <- factor(sample_data(AB_fecal_AMR_data.ps)$Combined_order, levels = c("Start_CONV","Midpoint_CONV",
                                                                                                                  "End_CONV","Start_RWA","Midpoint_RWA",
                                                                                                                  "End_RWA","Final_washout"))
# This is optional, but need to remove 1 factor level to make a square matrix possible for the trend test
#AB_fecal_AMR_data.ps <- subset_samples(AMR_data.ps, Combined_order != "Final_washout")

# Change variable in metadata to factor
sample_data(AB_fecal_AMR_data.ps)$participant_id <- as.factor(sample_data(AB_fecal_AMR_data.ps)$participant_id)

# Change label for most specific annotation level, in this case AMR gene group, to "Species"
colnames(phyloseq::tax_table(AB_fecal_AMR_data.ps))[4] <- "Species"


ANCOM_feces_output = ancombc2(data = AB_fecal_AMR_data.ps, assay_name = "counts", tax_level = "Species",
                              fix_formula = "Combined_order + participant_id + blinded_treatment",
                              p_adj_method = "holm", 
                              group = "Combined_order", struc_zero = TRUE, neg_lb = TRUE,verbose = TRUE,
                              global = TRUE, pairwise = TRUE, n_cl = 3)





AB_ancom_output_1 = ancombc2(data = AB_fecal_AMR_data.ps, assay_name = "counts", tax_level = "Species",
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




#res_prim = output$res
#tab_sens = output$pseudo_sens_tab
res_pair = AB_ancom_output_1$res_pair

#res_global = output$res_global
#res_dunn = output$res_dunn
#res_trend = output$res_trend

#### Calculate sensitivity scores for pool_timepoint ####
## Rename pairwise comparisons
df_Combined_order = res_pair %>%
  dplyr::select(taxon, contains("Combined_order")) 
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

# Change order of X-value
df_fig_Combined_order
df_fig_Combined_order$group <- factor(df_fig_Combined_order$group, levels = c("Start_RWA vs. Start_CONV","Midpoint_CONV vs. Start_CONV",
                                                                              "End_CONV vs. Midpoint_CONV","Midpoint_RWA vs. Start_RWA","End_RWA vs. Midpoint_RWA"))




# Plot heatmap
lo = floor(min(df_fig_Combined_order$value))
up = ceiling(max(df_fig_Combined_order$value))
mid = (lo + up)/2


AB_fig_organism = df_fig_Combined_order %>%
  ggplot(aes(x = group, y = reorder(taxon, value), fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "LogFC by combined order (AB)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
AB_fig_organism



