# Make sure we load this with the latest R packages

# Remember to purge modules and load the following 
#ml purge
#ml GCC/13.2.0 R_tamu/4.3.3-gfbf-2023b GSL/2.7 GMP/6.3.0 MPFR/4.2.1

# Run the next line if you want to load the results directly
#load(file = "/home/training/epi_on_the_island2024/shotgun/phyloseq/ANCOMBC2_results.RData")

# Otherwise, we can run this script to make sure we all start in the same place
#source("/home/training/epi_on_the_island2024/shotgun/phyloseq/1_LOAD_ALL.R")
library(ANCOMBC)
library(tidyr)

#
##
# Run ANCOMBC2 separately - Only keep this part for your dataset #####
##
#

# Here, we'll use "microbiome_sample" for analysis

# Specify the order of variables in a factor
# your first variable here will be used as the "reference" in the model
microbiome_DA_phylum.ps <- tax_glom(microbiome_sample, taxrank = "Phylum")
# Remove taxa with 0
microbiome_DA_phylum.ps <- prune_taxa(taxa_sums(microbiome_DA_phylum.ps) > 0, microbiome_DA_phylum.ps)
any(taxa_sums(microbiome_DA_phylum.ps)==0) # Double checkingthere are no empty taxa

# Confirm order of factors in collection_day
sample_data(microbiome_DA_phylum.ps)$collection_day <- factor(sample_data(microbiome_DA_phylum.ps)$collection_day, levels = c("Pre_weaning","At_weaning","Post_weaning"))

sample_data(microbiome_DA_phylum.ps)$collection_day

# Identify low abundance features
sort(taxa_sums(microbiome_DA_phylum.ps))

# Assuming `microbiome_DA_phylum.ps` is your phyloseq object
otu_table <- otu_table(microbiome_DA_phylum.ps)

# Create a logical matrix where TRUE indicates a count greater than zero
non_zero_counts <- otu_table > 0

# Sum the TRUE values across samples (assuming taxa are rows)
sample_presence_count <- apply(non_zero_counts, 1, sum)

# If you want to sort the results to see which taxa are present in more samples
sorted_sample_presence_count <- sort(sample_presence_count, decreasing = TRUE)
 
# ANCOMBC model code ####
# Adjust the "fix_formula" to match your study design".
# Change the "group" variable to whatever you want to run pairwise comparisons on.
# We can play around with the "rand_formula" later for repeated measures.
ancom_model_output = ancombc2(data = microbiome_DA_phylum.ps, assay_name = "counts", tax_level = "Phylum",
                              fix_formula = "collection_day", rand_formula = NULL,
                              p_adj_method = "holm", prv_cut = 0.10, lib_cut = 100, s0_perc = 0.05,
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
  #filter(`diff_collection_dayAt_weaning` == 1 | `diff_collection_dayPost_weaning_collection_dayAt_weaning` == 1 ) %>%
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
    group == "collection_dayAt_weaning" ~ "At weaning vs PreWeaning",
    group == "collection_dayPost_weaning_collection_dayAt_weaning" ~ "PostWeaning vs At weaning",
    FALSE ~ group # Keeps the original name for groups not specified
  )) %>%
  drop_na(group)  # Removes rows with NA in the 'group' column

df_fig_collection_day$group
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
  geom_text(aes(group, taxon, label = ifelse(value == 0, "", value),color = color), size = 3) + # don't show the 0s , #, color = color
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = "Log fold changes between groups") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_pair

#View(df_fig_collection_day)

#save.image(file = "/home/training/epi_on_the_island2024/shotgun/phyloseq/ANCOMBC2_results.RData")
 
# Try tornado plot ####

# Filter the results to only include significant features
significant_features = df_fig_collection_day$taxon

# Extract the names of the significant features
significant_feature_names = as.character(significant_features)

# Extract the bias_correct_log_table from the ancombc output
log_table = ancom_model_output$bias_correct_log_table

# Subset the table to include only the significant features
log_table_significant = log_table[rownames(log_table) %in% significant_feature_names, ]

# Convert the row names into a column
log_table_significant$Feature = rownames(log_table_significant)

# Reshape the data into a format suitable for ggplot2
df = reshape2::melt(log_table_significant, id.vars = "Feature")


# Calculating the average value for each Feature
median_values <- df %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(median_value = median(value, na.rm = TRUE)) # na.rm = TRUE to handle any NA values
median_values

df$value


# Merge the df data frame with the sample data
merged_df = merge(median_values, df_fig_collection_day, by.x = "Feature", by.y = "taxon")

colnames(merged_df)

merged_df$Feature <- factor(merged_df$Feature, levels = unique(merged_df$Feature))



# Assuming `sig_color` is already defined and you want to preserve existing colors for other groups
merged_df <- merged_df %>%
  mutate(sig_color = ifelse(
    group == "At weaning vs PreWeaning" & q_collection_dayAt_weaning < 0.05, "red",
    ifelse(group == "PostWeaning vs At weaning" & q_collection_dayPost_weaning_collection_dayAt_weaning < 0.05, "red", "black")
  ))

# Update the tornado_plot with new specifications
tornado_plot <- ggplot(merged_df, aes(x = value, y = Feature)) +
  geom_point(aes(size = log2(abs(median_value)+0.001), color = sig_color)) +  # Shape 21 is a filled circle
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Add a dashed line at x=0
  scale_color_identity() +  # Use actual color names specified in the sig_color column
  facet_grid(~ group, scales = "free_y") +  # Facets with free y-axis scale, arranged in one row
  scale_size_continuous(range = c(2, 8), guide = "none") +  # Adjust the size range as needed
  labs(x = "Value", y = "Feature") +
  theme_minimal() +
  theme(
    #panel.grid.major.x = element_blank(),  # Remove vertical grid lines
    panel.grid.minor.x = element_blank(),  # Remove minor vertical grid lines
    panel.grid.major.y = element_line(color = "grey", linetype = "dashed"),  # Keep horizontal grid lines
    axis.text.y = element_text(color = "black",lineheight = 1.5, size = 12),  # Uniform color for y-axis labels
    axis.ticks.y = element_line(color = "black"),  # Uniform color for y-axis ticks
    #axis.title.y = element_blank(),  # Removes y-axis label
    panel.spacing = unit(2, "lines"),  # Adjust spacing between panels
    strip.background = element_blank(),  # Optional: removes background from panel labels
    strip.text.x = element_text(face = "bold")  # Bold text for group names
  ) +
  scale_x_continuous(limits = c(-max(abs(merged_df$value)), max(abs(merged_df$value))))  # Ensuring zero is centered



print(tornado_plot)
