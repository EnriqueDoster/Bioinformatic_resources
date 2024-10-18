library(dada2); library(ggplot2); library(gridExtra)

# load environment
load("dada2.RData")


#------------------------------------------------------------------------
# Construct sequence table
#------------------------------------------------------------------------
seqtab <- makeSequenceTable(mergers)

# Further process ASVs (Amplicon Sequence Variants) by removing non-target-length sequences from your sequence table 
## This is analogous to cutting a band in-silico to get amplicons of the targeted length.
asvtab <- seqtab[,nchar(colnames(seqtab)) %in% 250:256] # Keep sequences of expected V4 length (optional for V3-V4 data)


#------------------------------------------------------------------------
# Remove chimeric sequences
#------------------------------------------------------------------------
asvtab.nochim <- removeBimeraDenovo(asvtab, method="consensus", multithread=TRUE, verbose=TRUE)


# save environment
save.image("dada2.RData")
