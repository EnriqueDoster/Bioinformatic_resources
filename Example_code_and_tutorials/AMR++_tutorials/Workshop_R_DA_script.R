#### SET WD ####
setwd("/home/training10/phyloseq_test/")

### SOURCE PREVIOUS SCRIPTS ####
source("1_load_libraries_data.R")
library(DAtest)
library(eulerr)


##
# Confirm data we want to use ####
##
#
sample_data(microbiome)$collection_day

#remove Na's
microbiome_DA.ps <- subset_samples(microbiome, collection_day != "NA")

microbiome_DA.ps <- prune_taxa(taxa_sums(microbiome_DA) > 0, microbiome_DA)
any(taxa_sums(microbiome_DA.ps)==0) # Double checkingthere are no empty taxa


# Use microbiome, agglomerate to family, name it something else
microbiome_DA_family.ps <- tax_glom(microbiome_DA.ps, taxrank = "Family")

# Confirm order of factors in collection_day
sample_data(microbiome_DA_family.ps)$collection_day <- factor(sample_data(microbiome_DA_family.ps)$collection_day, levels = c("Pre_weaning","At_weaning","Post_weaning"))


#
##
# Run test of all DA models, with multiple permutations and spike-ins to calculate FDR, AUC,etc ####
##
# Increase R to something like ~100 for more stable results, for quick tests low values are OK
test_DA_methods <- testDA(microbiome_DA_family.ps, predictor = "collection_day",R = 2, cores = 20)

# Summarize the model fit across the multiple runs
summary(test_DA_methods)

# You can also run all tests and store the results, without calculating how well the test performs
all_test <- allDA(microbiome_DA_family.ps, predictor = "collection_day", cores = 5)

# An example of pulling out the results for the ZIG model
all_test$results$zig

# Create Venn diagram of shared results
vennDA(all_test, tests=c("abc","kru"))

# Test the power for a particular test
power.abc <- powerDA(microbiome_DA_family.ps, predictor = "collection_day", test = "abc", cores = 20)
plot(power.abc)
summary(power.abc)

power.kru <- powerDA(microbiome_DA_family.ps, predictor = "collection_day", test = "kru", cores = 20)
plot(power.kru)
summary(power.kru)

# Arrange the plots side by side with labels
combined_plot <- plot_grid(
  plot(power.abc),
  plot(power.kru),
  labels = c("abc", "kru"),
  ncol = 2
)

combined_plot

#
##
# Run ANCOMBC2 separately - Only keep this part for your dataset #####
##
#

# Here, we'll use "microbiome" to 

# Specify the order of variables in a factor
# your first variable here will be used as the "reference" in the model
sample_data(microbiome)$collection_day <- factor(sample_data(microbiome)$collection_day ,levels = c("none", "nauplii+mysids","nauplii+mysids+seadragons"))


# ANCOMBC model code
# Adjust the "fix_formula" to match your study design".
# Change the "group" variable to whatever you want to run pairwise comparisons on.
# We can play around with the "rand_formula" later for repeated measures.
ancom_model_output = ancombc2(data = microbiome, assay_name = "counts", tax_level = "Family", 
                              fix_formula = "collection_day", rand_formula = NULL,
                              p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                              group = "collection_day", struc_zero = TRUE, neg_lb = TRUE,
                              alpha = 0.05, n_cl = 2, verbose = TRUE,
                              global = TRUE, pairwise = TRUE)


# Extract results from pairwise comparisons
res_pair = ancom_model_output$res_pair

# extract only results matching the variable used for pairwise comparisons
df_collection_day = res_pair %>%
  dplyr::select(taxon, contains("collection_day")) 

# Format data into long form for plotting
df_fig_collection_day <- df_collection_day %>%
  filter(`diff_collection_daynauplii+mysids` == 1 | `diff_collection_daynauplii+mysids+seadragons` == 1 | `diff_collection_daynauplii+mysids+seadragons_collection_daynauplii+mysids` == 1 ) %>%
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
    group == "collection_daynauplii+mysids" ~ "N+M vs none",
    group == "collection_daynauplii+mysids+seadragons" ~ "N+M+S vs none",
    group == "collection_daynauplii+mysids+seadragons_collection_daynauplii+mysids" ~ "N+M+S vs NM",
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

############################ DIFFERENTIAL ABUNDANCE (ANCOM-BC) ####
######## MICROBIOME
#### COLLECTION DAY AT GENUS LEVEL (At_weaning v. Pre-weaning, Post_weaning v. Pre_weaning)
ancombc_data <- tax_glom(microbiome_DA.ps, taxrank = "Genus") # need raw counts at the genus level
## MAKE Pre_weaning our 'reference' variable
ancombc_data@sam_data$collection_day <- factor(ancombc_data@sam_data$collection_day, levels = c("Pre_weaning","At_weaning","Post_weaning"))
ancombc_data@sam_data$collection_day # in the proper order now (Pre_weaning first)

## run ANCOM-BC2 for collection day
ancombc_collection_day <- ancombc2(ancombc_data, assay_name = "counts", tax_level = "Genus",
                                   fix_formula = "collection_day", rand_forumla = NULL, p_adj_method = "BH",
                                   group = "collection_day", prv_cut = 0.1, neg_lb = T, struc_zero = T,
                                   alpha = 0.05, n_cl = 1, verbose = T, global = T, pairwise = T)

## extract results from the pairwise comparison
res_pair <- ancombc_collection_day$res_pair

## melt res_pair into long-form table for plotting (### ENRIQUE THIS FILTER LINE WILL NEED TWEAKING BECAUSE I'M NOT SURE I KNOW WHAT THE NOMECLATURE WILL BE)
ancom_coll_day_fig <- res_pair %>%
  filter(`diff_collection_daynauplii+mysids` == 1 | `diff_collection_daynauplii+mysids+seadragons` == 1 | `diff_collection_daynauplii+mysids+seadragons_collection_daynauplii+mysids` == 1 ) %>%
  mutate(across(starts_with("lfc_collection_day"), ~ifelse(get(gsub("lfc_", "diff_", cur_column())) == 1, .x, 0)),
         across(starts_with("lfc_collection_day"), ~round(.x, 2), .names = "{.col}_rounded"),
         across(starts_with("diff_collectioan_day"), ~ifelse(.x == 1 & get(gsub("diff_", "passed_ss_", cur_column())) == 1, "black", "grey65"), .names = "{.col}_color")) %>%
  pivot_longer(cols = contains("_rounded"), names_to = "group", values_to = "value", names_prefix = "lfc_") %>%
  pivot_longer(cols = contains("_color"), names_to = "color_group", values_to = "color", names_prefix = "diff_") %>%
  filter(gsub("_rounded", "", group) == gsub("_color", "", color_group)) %>%
  select(-color_group) %>%
  arrange(taxon)
ancom_coll_day_fig

### Change the 'group' column to match LFC and colours
ancom_coll_day_fig$group # got 'rounded' at the end currently
ancom_coll_day_fig$group <- sub("_rounded","",ancom_coll_day_fig$group)
ancom_coll_day_fig$group # gone

## Rework our group names so they are shorter and more manageable
ancom_coll_day_fig <- ancom_coll_day_fig %>%
  mutate(group = case_when(
    group == 'collection_daynauplii+mysids' ~ "PrW vs AW", #Pre weaning v. at weaning
    group == "collection_daynauplii+mysids+seadragons" ~ "PoW vs AW", #Post weanin v. at weaning
    group == "collection_daynauplii+mysids+seadragons_collection_daynauplii+mysids" ~ "PoW vs PrW", # postweaning v. pre weaning
    TRUE ~ group # Keeps the original name for groups not specified
  ))
ancom_coll_day_fig$group # now they are much better

#### LET'S MAKE A HEATMAP!

## Calculate the range of log fold change for our colour gradient
low = min(ancom_coll_day_fig$value)
low # lowest log fold change is -5.41
high = max(ancom_coll_day_fig$value)
high # +5.39 is highest
mid = (low + high)/2
mid # -0.01 is the midpoint

### make the figure
diff_abund_hm <- ancom_coll_day_fig %>%
  ggplot(aes(x= group, y= taxon, fill = value)) +
  theme_minimal() +
  geom_tile(colour = "black") +
  scale_fill_gradient2('LFC', low = "dodgerblue3", high = "orangered3", mid = "white", midpoint = 0, limit = c(low,high)) +
  geom_text(aes(group, taxon, label = ifelse(value== 0, "", value), colour = color), size = 5, fontface = "bold") + # this turns off the cells with 0
  scale_color_identity(guide = "none") +
  theme(legend.title = element_text(face = "bold", size = 24),
        legend.text = element_text(colour = "black", size = 16),
        legend.key.size = unit(3,"lines"),
        legend.ticks = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        axis.title = element_blank(),
        axis.text.y= element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 24, colour = "black"))
diff_abund_hm # looks pretty good!