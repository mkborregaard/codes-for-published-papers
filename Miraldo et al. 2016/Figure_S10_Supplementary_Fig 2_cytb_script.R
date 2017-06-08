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


load("terrestrial mammals data object.RData")


#################


total_richness <- sitestat(ter, "mam_richness") + sitestat(ter, "amph_richness")
ter <- add_sitestat(ter, total_richness)

for(dataset in c("amph", "mam", "total")){
  
  dat <- read.csv(paste("Data/",dataset,"fig2.csv", sep = ""), stringsAsFactors = FALSE)
  bp_name <- paste(dataset,"_basepairs", sep = "")
  sp_name <- paste(dataset,"_gene_species", sep = "")
  IUCN_name <- paste(dataset,"_richness", sep = "")
  names(dat) <- c("site", bp_name, sp_name)
  #dat$site <- lat_long_to_long_lat(dat$site)
  
  ter <- add_sitestat(ter, dat)
  basepairs <- log(sitestat(ter, bp_name))
  sp_prop <- log(sitestat(ter, sp_name)/sitestat(ter, IUCN_name))
  
  Cairo(paste(dataset, "_FigS10.pdf", sep = ""), width = 10, height = 7, type = "pdf", unit = "in")
  two_color_map(ter, sp_prop, basepairs, 10, gridcol = "lightgrey", gridsites = which(!is.na(sp_prop)), overlay_shape = T, shapeborder = "white", shapefill = "grey", colour = "green_red")
  dev.off()
}

