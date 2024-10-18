library(dada2); library(ggplot2); library(gridExtra)

#------------------------------------------------------------------------
#Learn the Error Rates
#------------------------------------------------------------------------

# load environment
load("dada2.RData")

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Plot error profiles
forward.error.plot <- plotErrors(errF, nominalQ=F)
reverse.error.Plot <- plotErrors(errR, nominalQ=F)

error.plots <- grid.arrange(forward.error.plot, reverse.error.Plot, nrow = 1)

ggsave("error_plots.pdf", device = "pdf", error.plots, width = 7, height = 3)

save.image("dada2.RData")
