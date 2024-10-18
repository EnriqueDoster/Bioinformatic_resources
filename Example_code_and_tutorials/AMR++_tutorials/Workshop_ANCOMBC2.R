# Make sure we load this with the latest R packages
#ml purge
#ml GCC/13.2.0 R_tamu/4.3.3-gfbf-2023b GSL/2.7 GMP/6.3.0 MPFR/4.2.1

source("/home/training/epi_on_the_island2024/shotgun/phyloseq/1_LOAD_ALL.R")
library(ANCOMBC)
library(tidyr)
#
##
# Run ANCOMBC2 separately - Only keep this part for your dataset #####
##
#

# Here, we'll use "microbiome_sample" to 

# Specify the order of variables in a factor
# your first variable here will be used as the "reference" in the model

hist(sort(sample_sums(microbiome_sample)))
plot(sort(taxa_sums(microbiome_sample)))


microbiome_DA_family.ps <- tax_glom(microbiome_sample, taxrank = "Family")
# Remove taxa with 0
microbiome_DA_family.ps <- prune_taxa(taxa_sums(microbiome_DA_family.ps) > 0, microbiome_DA_family.ps)
any(taxa_sums(microbiome_DA_family.ps)==0) # Double checkingthere are no empty taxa

# Confirm order of factors in collection_day
sample_data(microbiome_DA_family.ps)$collection_day <- factor(sample_data(microbiome_DA_family.ps)$collection_day, levels = c("Pre_weaning","At_weaning","Post_weaning"))



# ANCOMBC model code
# Adjust the "fix_formula" to match your study design".
# Change the "group" variable to whatever you want to run pairwise comparisons on.
# We can play around with the "rand_formula" later for repeated measures.
ancom_model_output = ancombc2(data = microbiome_DA_family.ps, assay_name = "counts", 
                              fix_formula = "collection_day", rand_formula = NULL,
                              p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                              group = "collection_day", struc_zero = TRUE, neg_lb = TRUE,
                              alpha = 0.05, n_cl = 20, verbose = TRUE,
                              global = TRUE, pairwise = TRUE)


# Extract results from pairwise comparisons
res_pair = ancom_model_output$res_pair

# extract only results matching the variable used for pairwise comparisons
df_collection_day = res_pair %>%
  dplyr::select(taxon, contains("collection_day")) 

# Format data into long form for plotting
df_fig_collection_day <- df_collection_day %>%
  filter(`diff_collection_dayAt_weaning` == 1 | `diff_collection_dayPost_weaning_collection_dayAt_weaning` == 1 ) %>%
  mutate(across(starts_with("lfc_collection_day"), ~ifelse(get(gsub("lfc_", "diff_", cur_column())) == 1, .x, 0)),
         across(starts_with("lfc_collection_day"), ~round(.x, 2), .names = "{.col}_rounded"),
         across(starts_with("diff_collection_day"), ~ifelse(.x == 1 & get(gsub("diff_", "passed_ss_", cur_column())) == 1, "darkred", "darkgrey"), .names = "{.col}_color")) %>%
  pivot_longer(cols = contains("_rounded"), names_to = "group", values_to = "value", names_prefix = "lfc_") %>%
  pivot_longer(cols = contains("_color"), names_to = "color_group", values_to = "color", names_prefix = "diff_") %>%
  filter(gsub("_rounded", "", group) == gsub("_color", "", color_group)) %>%
  select(-color_group) %>%
  arrange(taxon) 

# Ensure 'group' column matches between LFC values and colors for direct comparison
df_fig_collection_day$group <- sub("_rounded", "", df_fig_collection_day$group)

# Adjust specific group names
df_fig_collection_day <- df_fig_collection_day %>%
  mutate(group = case_when(
    group == "diff_collection_dayAt_weaning" ~ "At weaning vs PreWeaning",
    group == "diff_collection_dayPost_weaning_collection_dayAt_weaning" ~ "PostWeaning vs At weaning",
    TRUE ~ group # Keeps the original name for groups not specified
  ))

#
## Heatmap plot ####
#
# Calculate LFC value ranges
lo = floor(min(df_fig_collection_day$value))
up = ceiling(max(df_fig_collection_day$value))
mid = (lo + up)/2

# Plot figure
fig_pair = df_fig_collection_day %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(lo, up),
                       na.value = "white", name = NULL) +
  geom_text(aes(group, taxon, label = ifelse(value == 0, "", value),color = color), size = 6) + # don't show the 0s , #, color = color
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = "Log fold changes between groups") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_pair

save.image(file = "ANCOMBC2_results.RData")

