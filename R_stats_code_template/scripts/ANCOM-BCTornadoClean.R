##Paul's Dunnett's type test####
library(ANCOMBC)
library(tidyr)

#Import your objects or data that is pruned but not normalized! My object was "foals", so you will see that below 

#Paul's request was to compare groups to the oldest foals
#Thought process here is that the microbiome is growing to some stable point, so the oldest calves are the control group
# Specify the order of variables in a factor
# your first variable here will be used as the "reference" in the model, I wanted mine to be d 120 (oldest foals) 
#ALSO make sure everything is a factor first if you have number based data
sample_data(foals)$timepoint_ref <- sample_data(foals)$time_point #made a new column called timepoint_ref that had my timepoints together. 
sample_data(foals)$timepoint_ref <- factor(sample_data(foals)$timepoint_ref, levels = c("120","0","2","7","14","21","28","60","90")) #Where to put your first group/factor as the reference 

#Run ANCOM-BC with the reference group set up, this was from the ANCOM-BC tutorial site on setting up longitudinal data, 
#so make sure you have the right model set-up for your data, yours will probably look different if it is not longitudinal 
foal_ancom_output_refchange = ancombc2(data = foals, assay_name = "counts", tax_level = "Family",
                                       fix_formula = "timepoint_ref", rand_formula = "(1|Animal.ID)",
                                       p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                       group = "timepoint_ref", struc_zero = TRUE, neg_lb = TRUE,
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
##Extract the results, you will need this later for the tornado plot too - so remember the object name
#Note, you can extract different results depending on what you want - for example they have a res_pair that you can extract pairwise comps
redo_foal_res_dunn = foal_ancom_output_refchange$res_dunn
#make a datatable to have for writing and to see if your plots come out right, I like the .csv better
redo_foal_res_dunn[1:6, ] %>%
  datatable(caption = "ANCOM-BC2 Dunnett's Type of Test")
#jk maybe just do a csv
write.csv(redo_foal_res_dunn, "Foal_Dunn_datatable.csv")
##########Data Prep for HEAT MAP########################
#I like doing a heat map since they are easier to interpret and I can compare to the tornado plot if necessary
#Run the dunnet code, this is filtering out anything that shows a significant change  
df_foal_fig_dunn = redo_foal_res_dunn %>%
  dplyr::filter(diff_timepoint_ref0 == 1 |
                  diff_timepoint_ref2 ==1|
                  diff_timepoint_ref7==1|
                  diff_timepoint_ref14 == 1|
                  diff_timepoint_ref21 == 1|
                  diff_timepoint_ref28 == 1|
                  diff_timepoint_ref60 == 1|
                  diff_timepoint_ref90 == 1) %>%
  dplyr::mutate(lfc_0_to_120 = ifelse(diff_timepoint_ref0 == 1, 
                                      round(lfc_timepoint_ref0, 2), lfc_timepoint_ref0),
                lfc_2_to_120 = ifelse(diff_timepoint_ref2 == 1, 
                                      round(lfc_timepoint_ref2, 2), lfc_timepoint_ref2),
                lfc_7_to_120 = ifelse(diff_timepoint_ref7 == 1,
                                      round(lfc_timepoint_ref7, 2),lfc_timepoint_ref7),
                lfc_14_to_120 = ifelse(diff_timepoint_ref14 == 1,
                                       round(lfc_timepoint_ref14, 2),lfc_timepoint_ref14),
                lfc_21_to_120 = ifelse(diff_timepoint_ref21 == 1,
                                       round(lfc_timepoint_ref21, 2),lfc_timepoint_ref21),
                lfc_28_to_120 = ifelse(diff_timepoint_ref28 == 1,
                                       round(lfc_timepoint_ref28, 2),lfc_timepoint_ref28),
                lfc_60_to_120 = ifelse(diff_timepoint_ref60 == 1,
                                       round(lfc_timepoint_ref60, 2),lfc_timepoint_ref60),
                lfc_90_to_120 = ifelse(diff_timepoint_ref90 == 1,
                                       round(lfc_timepoint_ref90, 2),lfc_timepoint_ref90)) %>% 
  transmute(taxon, 
            `0 vs. 120` = round(lfc_0_to_120, 2),
            `2 vs. 120` = round(lfc_2_to_120, 2),
            `7 vs. 120` = round(lfc_7_to_120, 2),
            `14 vs. 120` = round(lfc_14_to_120, 2),
            `21 vs. 120` = round(lfc_21_to_120, 2),
            `28 vs. 120` = round(lfc_28_to_120, 2),
            `60 vs. 120` = round(lfc_60_to_120, 2),
            `90 vs. 120` = round(lfc_90_to_120, 2)) %>%
  tidyr::pivot_longer(cols = "0 vs. 120":"2 vs. 120":"7 vs. 120":"14 vs. 120":"21 vs. 120":"28 vs. 120":"60 vs. 120":"90 vs. 120", 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)

#Sensitivity testing, this is assigning colors to LogFC that do or do not pass sensitivity testing
df_sens_foal_fig_dunn = redo_foal_res_dunn %>%
  dplyr::filter(diff_timepoint_ref0 == 1 |
                  diff_timepoint_ref2 ==1|
                  diff_timepoint_ref7==1|
                  diff_timepoint_ref14 == 1|
                  diff_timepoint_ref21 == 1|
                  diff_timepoint_ref28 == 1|
                  diff_timepoint_ref60 == 1|
                  diff_timepoint_ref90 == 1) %>%
  dplyr::mutate(lfc_0_to_120 = ifelse(passed_ss_timepoint_ref0 == 1 & diff_timepoint_ref0 == 1, 
                                      "white", "black"),
                lfc_2_to_120 = ifelse(passed_ss_timepoint_ref2 == 1 & diff_timepoint_ref2 == 1, 
                                      "white", "black"),
                lfc_7_to_120 = ifelse(passed_ss_timepoint_ref7 == 1 & diff_timepoint_ref7 == 1, 
                                      "white", "black"),
                lfc_14_to_120 = ifelse(passed_ss_timepoint_ref14 == 1 & diff_timepoint_ref14 == 1, 
                                       "white", "black"),
                lfc_21_to_120 = ifelse(passed_ss_timepoint_ref21 == 1 & diff_timepoint_ref21 == 1, 
                                       "white", "black"),
                lfc_28_to_120 = ifelse(passed_ss_timepoint_ref28 == 1 & diff_timepoint_ref28 == 1, 
                                       "white", "black"),
                lfc_60_to_120 = ifelse(passed_ss_timepoint_ref60 == 1 & diff_timepoint_ref60 == 1, 
                                       "white", "black"),
                lfc_90_to_120 = ifelse(passed_ss_timepoint_ref90 == 1 & diff_timepoint_ref90 == 1, 
                                       "white", "black"),) %>%
  transmute(taxon, 
            `0 vs. 120` = lfc_0_to_120,
            `2 vs. 120` = lfc_2_to_120,
            `7 vs. 120` = lfc_7_to_120,
            `14 vs. 120` = lfc_14_to_120,
            `21 vs. 120` = lfc_21_to_120,
            `28 vs. 120` = lfc_28_to_120,
            `60 vs. 120` = lfc_60_to_120,
            `90 vs. 120` = lfc_90_to_120) %>% 
  
  tidyr::pivot_longer(cols = "0 vs. 120":"2 vs. 120":"7 vs. 120":"14 vs. 120":"21 vs. 120":"28 vs. 120":"60 vs. 120":"90 vs. 120", 
                      names_to = "group", values_to = "color") %>%
  dplyr::arrange(taxon)

#Plot the dunnet results in a figure 
#Join your sensitivty analysis with your ANCOM-BC
df_foal_fig_dunn1 = df_foal_fig_dunn%>%
  dplyr::left_join(df_sens_foal_fig_dunn, by = c("taxon", "group")) #add in sensitivity testing here

# Change order of X-value
df_foal_fig_dunn1$group <- factor(df_foal_fig_dunn1$group, levels = c("0 vs. 120","2 vs. 120","7 vs. 120","14 vs. 120","21 vs. 120","28 vs. 120","60 vs. 120","90 vs. 120"))

#Assign your map colors to values
lo = floor(min(df_foal_fig_dunn1$value))
up = ceiling(max(df_foal_fig_dunn1$value))
mid = 0
####HeatMap of Dunn Test####
foal_fig_dunn = df_foal_fig_dunn1 %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "black", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value, color = color), size = 4) +
  scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = "Log fold changes as foals aged compared to the most stable population (120 d)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
foal_fig_dunn

#############Tornado Plot Data Prep and Creation####################
# Format data into long form for plotting
#create a new object because you're scared to fuck up the other one 
volcano_resdunn_foal <- redo_foal_res_dunn

#This will do some filtering on your data to make it ready for your tornado plot
df_ed_long_volcano_dunn<- volcano_resdunn_foal %>%
  mutate(across(starts_with("lfc_timepoint_ref"), ~ifelse(get(gsub("lfc_", "diff_", cur_column())) == 1, .x, .x)), #FOR THE LOVE OF GOD PUT THIS AS .X OR IT JUST TAKE ALL THINGS TO 0, it won't show log fold changes 
         across(starts_with("lfc_timepoint_ref"), ~round(.x, 2), .names = "{.col}_rounded"),
         across(starts_with("diff_timepoint_ref"), #This whole string makes sure that it applies NA to the values that don't pass sensitivity testing, thank you Valeria
                ~ifelse(.x == 1 & get(gsub("diff_", "passed_ss_", cur_column())) == 1, "red3",
                        ifelse(.x == 1 | get(gsub("diff_", "passed_ss_", cur_column())) == 0, NA, "black")), 
                .names = "{.col}_color")
  )%>%
  pivot_longer(cols = contains("_rounded"), names_to = "group", values_to = "value", names_prefix = "lfc_",values_drop_na = FALSE) %>%
  pivot_longer(cols = contains("_color"), names_to = "color_group", values_to = "color", names_prefix = "diff_",values_drop_na = FALSE) %>%
  filter(gsub("_rounded", "", group) == gsub("_color", "", color_group)) %>%
  select(-color_group) %>%
  arrange(taxon)%>%
  mutate(value = ifelse(is.na(color), NA, value)) #here is where it actually makes sure the value AND color say "NA"

# Ensure 'group' column matches between LFC values and colors for direct comparison
df_ed_long_volcano_dunn$group <- sub("_rounded", "", df_ed_long_volcano_dunn$group)
#Just double check that your groups turned out correct
df_ed_long_volcano_dunn$group
unique(df_ed_long_volcano_dunn$group)

# Try tornado plot ####

# Filter the results to only include significant features
significant_features = df_ed_long_volcano_dunn$taxon

# Extract the names of the significant features
significant_feature_names = as.character(significant_features)

# Extract the bias_correct_log_table from the ancombc output
log_table = foal_ancom_output_refchange$bias_correct_log_table

# Subset the table to include only the significant features
log_table_significant = log_table[rownames(log_table) %in% significant_feature_names, ]

# Convert the row names into a column
log_table_significant$Feature = rownames(log_table_significant)

# Reshape the data into a format suitable for ggplot2
df_ed = reshape2::melt(log_table_significant, id.vars = "Feature")

# Calculating the average value for each Feature
median_values <- df_ed %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(median_value = median(value, na.rm = TRUE)) # na.rm = TRUE to handle any NA values, REMOVE NA OR ELSE IT WILL ALL BE NA!
median_values

df_ed$value


# Merge the df data frame with the sample data
ed_merged_df = merge(median_values, df_ed_long_volcano_dunn, by.x = "Feature", by.y = "taxon")

colnames(ed_merged_df)

ed_merged_df$Feature <- factor(ed_merged_df$Feature, levels = unique(ed_merged_df$Feature))

# Update the tornado_plot with new specifications
tornado_plot <- ggplot(ed_merged_df, aes(x = value, y = Feature)) +
  geom_point(aes(size =log2(abs(median_value)+15), color = color)) +# Shape 21 is a filled circle
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Add a dashed line at x=0
  scale_color_identity() +
  facet_grid(~ factor(group,levels = c("timepoint_ref0","timepoint_ref2","timepoint_ref7","timepoint_ref14","timepoint_ref21","timepoint_ref28","timepoint_ref60","timepoint_ref90")), scales = "free_y") +  # Facets with free y-axis scale, arranged in one row
  scale_size_continuous(range = c(2, 8), guide = "none") +  # Adjust the size range as needed
  labs(x = "LogFC Value", y = "Feature") +
  theme_minimal() +
  theme(
    #panel.grid.major.x = element_blank(),  # Remove vertical grid lines
    panel.grid.minor.x = element_blank(),  # Remove minor vertical grid lines
    panel.grid.major.y = element_line(color = "grey", linetype = "dashed"),  # Keep horizontal grid lines
    axis.text.y = element_text(color = "black",lineheight = 1.5, size = 12),  # Uniform color for y-axis labels
    axis.ticks.y = element_line(color = "black"),  # Uniform color for y-axis ticks
    #axis.title.y = element_blank(),  # Removes y-axis label
    panel.spacing = unit(1, "lines"),  # Adjust spacing between panels
    strip.background = element_blank(),  # Optional: removes background from panel labels
    strip.text.x = element_text(face = "bold"))+ # Bold text for group names
  scale_x_continuous(limits = c(-10, 15)) #I manually adjusted my x-axis values so the dots didn't get cut off, if you don't care, use following line of code
#scale_x_continuous(limits = c(-max(abs(ed_merged_df$value)), max(abs(ed_merged_df$value))))  # Ensuring zero is centered

print(tornado_plot)

