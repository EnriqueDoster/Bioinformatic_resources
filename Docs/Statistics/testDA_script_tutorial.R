# Load packages ####
library(phyloseq)
library(ggplot2)
library(DAtest)
library(ANCOMBC)
library(dplyr)
library(tidyr)
library(eulerr)
#install.packages("eulerr")

# Run command to load RData file
# Change this to the path to your full dataset after the demonstration
load("/scratch/group/vlcs_689/lee_test/subsampled/phyloseq/SPX1.RData")

# Another option, is to run the script from here (your script has to be optimized for this):
# Unifrac can take hours, for example. 

# source("/scratch/group/vlcs_689/lee_test/subsampled/phyloseq/SPX1.R")


#
##
# Confirm data we want to use ####
##
#
sample_data(data4)$animals

# Use data4, agglomerate to family, name it something else
data4_DA_family <- tax_glom(data4, taxrank = "Family")

#
##
# Run test of all DA models, with multiple permutations and spike-ins to calculate FDR, AUC,etc ####
##
# Increase R to something like ~100 for more stable results, for quick tests low values are OK
test_DA_methods <- testDA(data4_DA_family, predictor = "animals",R = 2, cores = 5)

# Summarize the model fit across the multiple runs
summary(test_DA_methods)

# You can also run all tests and store the results, without calculating how well the test performs
all_test <- allDA(data4_DA_family, predictor = "animals", cores = 5)

# An example of pulling out the results for the ZIG model
all_test$results$zig

# Create Venn diagram of shared results
vennDA(all_test, tests=c("zig","kru"))

# Test the power for a particular test
power.zig <- powerDA(data4_DA_family, predictor = "animals", test = "zig", cores = 5)
plot(power.zig)
summary(power.zig)

#
##
# Run ANCOMBC2 separately - Only keep this part for your dataset #####
##
#

# Here, we'll use "data4" to 

# Specify the order of variables in a factor
# your first variable here will be used as the "reference" in the model
sample_data(data4)$animals <- factor(sample_data(data4)$animals ,levels = c("none", "nauplii+mysids","nauplii+mysids+seadragons"))


# ANCOMBC model code
# Adjust the "fix_formula" to match your study design".
# Change the "group" variable to whatever you want to run pairwise comparisons on.
# We can play around with the "rand_formula" later for repeated measures.
ancom_model_output = ancombc2(data = data4, assay_name = "counts", tax_level = "Family", 
                              fix_formula = "animals", rand_formula = NULL,
                              p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                              group = "animals", struc_zero = TRUE, neg_lb = TRUE,
                              alpha = 0.05, n_cl = 2, verbose = TRUE,
                              global = TRUE, pairwise = TRUE)


# Extract results from pairwise comparisons
res_pair = ancom_model_output$res_pair

# extract only results matching the variable used for pairwise comparisons
df_animals = res_pair %>%
  dplyr::select(taxon, contains("animals")) 

# Format data into long form for plotting
df_fig_animals <- df_animals %>%
  filter(`diff_animalsnauplii+mysids` == 1 | `diff_animalsnauplii+mysids+seadragons` == 1 | `diff_animalsnauplii+mysids+seadragons_animalsnauplii+mysids` == 1 ) %>%
  mutate(across(starts_with("lfc_animals"), ~ifelse(get(gsub("lfc_", "diff_", cur_column())) == 1, .x, 0)),
         across(starts_with("lfc_animals"), ~round(.x, 2), .names = "{.col}_rounded"),
         across(starts_with("diff_animals"), ~ifelse(.x == 1 & get(gsub("diff_", "passed_ss_", cur_column())) == 1, "darkred", "darkgrey"), .names = "{.col}_color")) %>%
  pivot_longer(cols = contains("_rounded"), names_to = "group", values_to = "value", names_prefix = "lfc_") %>%
  pivot_longer(cols = contains("_color"), names_to = "color_group", values_to = "color", names_prefix = "diff_") %>%
  filter(gsub("_rounded", "", group) == gsub("_color", "", color_group)) %>%
  select(-color_group) %>%
  arrange(taxon) 

# Ensure 'group' column matches between LFC values and colors for direct comparison
df_fig_animals$group <- sub("_rounded", "", df_fig_animals$group)

# Adjust specific group names
df_fig_animals <- df_fig_animals %>%
  mutate(group = case_when(
    group == "animalsnauplii+mysids" ~ "N+M vs none",
    group == "animalsnauplii+mysids+seadragons" ~ "N+M+S vs none",
    group == "animalsnauplii+mysids+seadragons_animalsnauplii+mysids" ~ "N+M+S vs NM",
    TRUE ~ group # Keeps the original name for groups not specified
  ))

  #
## Heatmap plot ####
#
# Calculate LFC value ranges
lo = floor(min(df_fig_animals$value))
up = ceiling(max(df_fig_animals$value))
mid = (lo + up)/2

# Plot figure
fig_pair = df_fig_animals %>%
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

#
##
# Compare ANCOMBC2 to ZIG model ####
##
# This part is optional!!!

# Extract results from testDA
# Assuming p-values are stored and we consider p < 0.05 as significant
significant_zig <- all_test$results$zig %>% filter(pval.adj < 0.05) %>% .$Feature

# Extract ANCOMBC2 results
ancombc_sig_tax <- unique(df_fig_animals$taxon)

# Convert the tax_table to a data frame
tax_table_df <- as.data.frame(tax_table(data4_DA_family))

# Make sure the row names (Feature IDs) are a column
tax_table_df$FeatureID <- rownames(tax_table_df)

# Assuming the family names are stored in a column named "Family"
# Extract only the FeatureID and Family columns
feature_family_df <- tax_table_df %>%
  dplyr::select(FeatureID, Family)

# Filter for significant features
significant_features_df <- feature_family_df %>%
  dplyr::filter(FeatureID %in% significant_zig)

# Convert them to unique sets if they aren't already
significant_families <- unique(significant_features_df$Family)
ancombc_families <- unique(ancombc_sig_tax)

# Prepare a list for eulerr
sets <- list(
  ZIG_Families = significant_families,
  ANCOMBC_Families = ancombc_families
)

# Generate Euler diagram
fit <- euler(sets)
plot(fit, quantities=TRUE)

#
## Compare again, but with filtering of low abundance and sparse features
#

# use preDA to filter out 
data4_DA_family.filtered <- preDA(data4_DA_family, min.samples = 9, min.reads = 10, min.abundance = 0) 

# Run all tests again
all_test.filtered <- allDA(data4_DA_family.filtered, predictor = "animals", cores = 5)

# Assuming p-values are stored and we consider p < 0.05 as significant
significant_zig.filtered  <- all_test.filtered$results$zig %>% filter(pval.adj < 0.05) %>% .$Feature

# Filter for significant features
significant_features_df.filtered  <- feature_family_df %>%
  dplyr::filter(FeatureID %in% significant_zig.filtered)

# Convert them to unique sets if they aren't already
significant_families.filtered <- unique(significant_features_df.filtered$Family)

# Prepare a list for eulerr
sets.filtered <- list(
  ZIG_Families_filtered = significant_families.filtered,
  ANCOMBC_Families = ancombc_families
)

# Generate Euler diagram
fit.filtered <- euler(sets.filtered)
plot(fit.filtered, quantities=TRUE)

