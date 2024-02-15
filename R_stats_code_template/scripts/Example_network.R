library("ape")
### building trees ##
random_tree = rtree(ntaxa(data_noSNP), rooted=TRUE, tip.label=taxa_names(data_noSNP))
plot(random_tree)
data_noSNP=merge_phyloseq(data_noSNP,random_tree)


# Some plots
plot_net(data_noSNP , maxdist=0.4, color="sample_type")
plot_tree(data_noSNP, color="DOF_categories", label.tips="taxa_names", ladderize="left", plot.margin=0.3)
plot_tree(data_noSNP, color="Group", shape="Time", label.tips="taxa_names", ladderize="right", plot.margin=0.3)

ig = make_network(data_noSNP, type = "samples", distance = "bray", max.dist = 0.2)
plot_network(ig, data_noSNP, color = "DOF_categories", shape = "DDD_cat", line_weight = 0.4, label = NULL)
plot_network(ig, data_noSNP, color = "DDD_cat", shape = "DOF_categories", line_weight = 0.4, label = NULL,layout.method=layout.bipartite)

# layout.fruchterman.reingold
# layout.reingold.tilford
# layout.bipartite

## Top 10 features
data_noSNP_100 = prune_taxa(names(sort(taxa_sums(data_noSNP), TRUE))[1:6], data_noSNP)
jg = make_network(data_noSNP_100, "taxa", "jaccard", 0.3)
plot_network(jg, data_noSNP_100, "taxa",line_weight = 0.4, label = NULL)
