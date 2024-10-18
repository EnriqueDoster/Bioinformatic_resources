
# load libraries
library(dada2); library(ggplot2); library(gridExtra)

# define path to trimmed reads
path.cut <- "../cutadapt/"

# defined primer removed forward and reverse fastq filenames format, check format:
cutFs <- sort(list.files(path.cut, pattern="_R1_001_trimmed.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001_trimmed.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have the format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))

#------------------------------------------------------------------------
# Inspect read quality profiles
#------------------------------------------------------------------------
p1 <- plotQualityProfile(cutFs[1:20], aggregate = TRUE) # this is combined plot can be done for all samples
p2 <- plotQualityProfile(cutRs[1:20], aggregate = TRUE) 
p3 <- grid.arrange(p1, p2, nrow = 1)

ggsave("plot_qscores_cutadapt.pdf", device= "pdf", p3, width = 7, height = 3)

save.image("dada2.RData")
