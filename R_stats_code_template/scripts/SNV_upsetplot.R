# Upset plots
#to install Microbiota
#BiocManager::install("MicrobiotaProcess")
library(UpSetR)
library("MicrobiotaProcess")


group_AMR.ps <- tax_glom(SNV_data_noSNP, taxrank = "group")
class_AMR.ps <- tax_glom(SNV_data_noSNP, taxrank = "class")


upsetda_edit <- get_upset(SNV_data_noSNP, factorNames="Combined_order") ## ASV



upset(upsetda_edit, sets=c("Start_CONV","Midpoint_CONV","End_CONV", "Start_RWA",
                           "Midpoint_RWA","End_RWA","Final_washout","Conventional Beef","RWA Beef"),
      
      sets.bar.color = c("#ff0000", "#cc0000", "#990000", "#0000ff", "#0000cc", "#000099", "#808080", "#ff6347", "#4169e1")
,text.scale = 2,
      
      order.by = "freq", empty.intersections = "on")

