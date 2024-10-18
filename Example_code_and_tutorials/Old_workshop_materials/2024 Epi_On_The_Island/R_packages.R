#-----------------------------------------------------------------------------------------------
# Install the latest version of R and RStudio on your local machine (Windows or Mac)
# Note that you must have R 4.2.0 or newer for to install the most current release 
# (https://posit.co/download/rstudio-desktop/)
#Install these packages and load libraries in your RStudio. We may not use all of these packages
#-----------------------------------------------------------------------------------------------



#---------------------------------------------
# Install microbiome analysis packages
#---------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("devtools")      # Install devtools for installing packages from GitHub

BiocManager::install("dada2")     # For 16S analysis follow https://benjjneb.github.io/dada2/index.html
                                  # https://github.com/benjjneb/dada2/releases v1.26
                                  
BiocManager::install("phyloseq")  # For handling microbiome data sets
install.packages("microbiome")    # For microbiome analysis (e.g. alphadiversity)
install.packages("vegan")         # Community ecology package
install.packages("qiime2R")       # To integrate QIIME2 outputs with R
BiocManager::install("Biostrings") # Manipulation of biological strings
BiocManager::install("ShortRead") # Reading, writing and manipulating short-read sequencing data in FASTQ format
install.packages("ape")           # Tools for analyzing phylogenetics distances
BiocManager::install("DECIPHER")  # Analysis/manipuation of sequence data
install.packages("decontam")      # For detecting and removing contaminant sequences


#---------------------------------------------
# Install packages for data manipulation 
#---------------------------------------------
install.packages("dplyr")    # For data manipulation
install.packages("tidyr")    # For data tidying
install.packages("tidyverse") # Data manipulation
install.packages("plyr")      # For data manipulation and transformation
install.packages("forcats")  # To manages factor variables
install.packages("matrixStats") # Functions for row and column operations on matrices
install.packages("metagenomeSeq") # For count data normalization (CSS) and statistical analysis 

#---------------------------------------------
# Install packages for statistical analysis
#---------------------------------------------
install.packages("Rmisc")      # To generate summary statistics (e.g. mean, SE, CI etc)
install.packages("lme4")       # For fitting linear mixed-effects models
install.packages("lmerTest")    # Extends lme4 to provide p-values and ANOVA tables
install.packages("car")        # For Type III ANOVA tests
install.packages("emmeans")    # For estimated marginal means (post-hoc analysis)
install.packages("nlme")       # For fitting non-linear and linear mixed-effects models
BiocManager::install("multcomp")  # For simultaneous tests and confidence intervals for general linear hypotheses
BiocManager::install("Maaslin2")  # For multivariable association analysis in microbiome
devtools::install_github("pedroembraga/pairwiseAdonis") # For pairwise adonis/PERMANOVA analysis (i.e. multivariate analysis)

#---------------------------------------------
# Install packages for data visualization
#---------------------------------------------
install.packages("ggthemes") # Themes and color palettes
install.packages("RColorBrewer") # Color palettes for data visualization
install.packages("viridis")  # Color palettes
install.packages("scales")  # Scaling functions for data visualization
install.packages("ggplot2")  # For data visualization
install.packages("pheatmap") # heatmaps with clustering and annotation options


#---------------------------------------------
# Load libraries for microbiome analysis
#---------------------------------------------
library(devtools)  
library(dada2); packageVersion("dada2") 
library(phyloseq)
library(microbiome)
library(vegan)
library(qiime2R)
library(Biostrings)
library(ShortRead)
library(ape)
library(DECIPHER)
library(decontam)
library(dplyr)
library(tidyr)
library(matrixStats)
library(tibble)
library(tidyverse)
library(plyr)
library(metagenomeSeq)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(nlme)
library(multcomp)
library(Maaslin2)
library(pairwiseAdonis)
library(ggthemes)
library(RColorBrewer)
library(viridis)
library(scales)
library(ggplot2)
library(pheatmap)
#---------------------------------------------



