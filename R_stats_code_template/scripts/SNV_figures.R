remotes::install_github("cpauvert/psadd")
library(psadd)
plot_krona(SNV_data_noSNP,"SNV-krona", "specific_sample_type",trim=TRUE)

krona_combined_order <- plot_krona(SNV_data_noSNP,"SNV-krona-Combined_order", "Combined_order",trim=TRUE)


plot_sparsity(SNV_data_noSNP, title = NULL)