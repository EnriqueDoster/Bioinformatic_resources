# Load the required packages
library(phyloseq)
library(bnlearn)

# Make sure your phyloseq objects have been created: data_micro.ps and data_noSNP

# Merge the original phyloseq objects
merged_ps <- merge_phyloseq(data_micro.ps_subset, data_noSNP_subset)

fecal_merged_ps <- subset_samples(merged_ps, sample_type == "Fecal")

# Transform the merged phyloseq object into a data frame
merged_otu_table <- as.data.frame(otu_table(fecal_merged_ps))
merged_sample_data <- data.frame(sample_data(fecal_merged_ps))

merged_sample_data <-  merged_sample_data[, c("Combined_order", "treatment_order","blinded_treatment","pool_timepoint")]

# Transpose the OTU table and convert it to a data frame
transposed_otu_table <- as.data.frame(t(merged_otu_table))

# Set the row names of the transposed OTU table to match the sample identifiers
rownames(transposed_otu_table) <- colnames(merged_otu_table)

# Merge the OTU table and sample data by matching row names (sample identifiers)
merged_df <- cbind(transposed_otu_table, merged_sample_data)

merged_df <- na.omit(merged_df)


columns <- colnames(merged_df)
merged_df[, columns] <- lapply(columns, function(x) as.numeric(merged_df[[x]]))
merged_df.dedup <- dedup(merged_df, .95, debug = FALSE)

set.seed(42)
sink(file = "learning-log-hc.boot.strength.txt")

dag.hybrid.group <- hc(merged_df)
cl <- parallel::makePSOCKcluster(10,outfile="debug.txt")
boot.hc.hybrid.group = boot.strength(data=counts,R=1000,algorithm="hc", algorithm.args = list(cluster=cl))
save.image("temp_XIT_phylum_1000rep_network.RData")

avg.boot.hc.hybrid.group <- averaged.network(boot.hc.hybrid.group, threshold = .7)
avg.boot.hc.dag.hybrid.group <- cextend(avg.boot.hc.hybrid.group)
fitted.hybrid.group <- bn.fit(avg.boot.hc.dag.hybrid.group,counts)

sink()

stopCluster(cl)
save.image("XIT_phylum_1000rep_network.RData")

# avg.boot.hc.hybrid.group <- averaged.network(boot.hc.hybrid.group, threshold = .7)
# avg.boot.hc.dag.hybrid.group <- cextend(avg.boot.hc.hybrid.group)
# fitted.hybrid.group <- bn.fit(avg.boot.hc.dag.hybrid.group,counts)
# 
# #Visualize it:
# d.gph.hybrid.group <- graphviz.plot(avg.boot.hc.dag.hybrid.group)
# pdf(file="final_network_hybrid_group.pdf",width=20,height=20)
# plot(d.gph.hybrid.group, nodeAttrs=makeNodeAttrs(d.gph.hybrid.group, fontsize=34))
# dev.off()

#set.seed(1)
#cpquery(fitted.hybrid.group,event=(PREVCAT_ALL == "Top3rd"),evidence=(Group == "Treatment"),n=1000)


