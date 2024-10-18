library(dada2); library(ggplot2); library(gridExtra)

# load environment
load("dada2.RData")

#------------------------------------------------------------------------
# Dereplicate sequences
#------------------------------------------------------------------------
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names


#------------------------------------------------------------------------
# Sample inference 
#------------------------------------------------------------------------
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)


#------------------------------------------------------------------------
# Merging paired-end reads
#------------------------------------------------------------------------
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

save.image("dada2.RData")
