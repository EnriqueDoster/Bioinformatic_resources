# Load packages
library(ggplot2)
library(dplyr)
library(plyr)
#library(tidyr)
library(stringr)
library(scales)
library(vegan)
library(knitr)
library(readr)
library(kableExtra)
library(randomcoloR)
library(GUniFrac)
library(ggdendro)
library(cowplot)
library(stringr)
library(ggthemes)
library(tidyverse)
library(ggsci)
library(ggdendro)
# install.packages("ggthemes")

# Requires special installations
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("phyloseq")
#BiocManager::install("metagenomeSeq")

library(phyloseq) # BiocManager
#library(metagenomeSeq) # BiocManager

# install with devtools
#install.packages('devtools')
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
#library(btools) # github
library(pairwiseAdonis) # github

#install.packages("remotes")
#remotes::install_github("vmikk/metagMisc")
library(metagMisc) # gitHub

# Optional for figures with microViz
# BiocManager::install("microbiome")
# BiocManager::install("ComplexHeatmap") # optional
# install.packages("microViz",repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos")))

# Custom
# source some scripts
# source("Final_R_analysis/16S/changeSILVATaxaNames.R")
# source("Final_R_analysis/16S/R_analysis/MergeLowAbund.R")
# source("Final_R_analysis/16S/uw_unifrac.R")
# source("Final_R_analysis/16S/w_unifrac.R")
# 
# 


changeSILVAtaxa_w_species <- function(x) {
  # remove the D__ etc...
  tax.clean <- data.frame(row.names = row.names(x),
                          Domain = str_replace(x[,1], "d__",""),
                          Kingdom = str_replace(x[,2], "k__",""),
                          Phylum = str_replace(x[,3], "p__",""),
                          Class = str_replace(x[,4], "c__",""),
                          Order = str_replace(x[,5], "o__",""),
                          Family = str_replace(x[,6], "f__",""),
                          Genus = str_replace(x[,7], "g__",""),
                          Species = str_replace(x[,8], "s__",""),
                          stringsAsFactors = FALSE)
}
