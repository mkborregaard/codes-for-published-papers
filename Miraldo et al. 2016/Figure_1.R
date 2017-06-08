########################################
# Codes for the paper:
# An Anthropocene Map of Genetic Diversity
#
# Andreia Miraldo, Sen Li, Michael K. Borregaard, Alexander Floréz-Rodriguéz, Shyam Gopalakrishnan, Mirneza Risvanovic, Zhiheng Wang, Carsten Rahbek, Katharine A. Marske & David Nogués-Bravo
#
# Submitted to Science, 2016
# Code in this file by Michael K. Borregaard
########################################

library(nodiv) 
library(Cairo)

load("../terrestrial mammals data object.RData")


colors <- read.table("plotcols.txt", header = T)
plotcols <- rgb(colors, max = 255)

Fig1 <- read.csv("../Data/equalarea_numbers.csv", stringsAsFactors = FALSE)
ter <- add_sitestat(ter, Fig1)

for(group in c("total","amphibians","mammals")){
  Cairo(paste(group,".pdf", sep = ""), type = "pdf", width = 10, height = 7, units = "in")
  #plotval = log(sitestat(ter,"Fig1GD"))
  plotval <- rank(sitestat(ter, paste(group, "GDval", sep = "")),na.last = "keep")
  plotval[abs(plotval) == Inf] <- NA
  plot_sitestat(ter, plotval, col = plotcols, legend = F, gridlwd = 0.4, gridcol = "white", overlay_shape = TRUE, gridsites = which(!is.na(plotval)))
  dev.off()
}
