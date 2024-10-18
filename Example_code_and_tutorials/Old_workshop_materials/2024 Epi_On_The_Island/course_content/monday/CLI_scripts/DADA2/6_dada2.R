library(dada2); library(ggplot2); library(gridExtra)

# load environment
load("dada2.RData")

#------------------------------------------------------------------------
# Track sequence counts through the pipeline
#------------------------------------------------------------------------
getN <- function(x) sum(getUniques(x))
track <- cbind(trimmed_out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(asvtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.csv(track, "dada2_track_16S.csv")

# Write summary table
summary_tab <- data.frame(dada2_input=trimmed_out[,1],
                          filtered=trimmed_out[,2], dada_f=sapply(dadaFs, getN),
                          dada_r=sapply(dadaRs, getN), merged=sapply(mergers, getN),
                          nonchim=rowSums(asvtab.nochim),
                          final_perc_reads_retained=round(rowSums(asvtab.nochim)/trimmed_out[,1]*100, 1))

write.csv(summary_tab,"summary_table_DADA2.csv")

# save environment
save.image("dada2.RData")
