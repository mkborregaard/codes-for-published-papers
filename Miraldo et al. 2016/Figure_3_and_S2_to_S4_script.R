########################################
# Codes for the paper:
# An Anthropocene Map of Genetic Diversity
#
# Andreia Miraldo, Sen Li, Michael K. Borregaard, Alexander Floréz-Rodriguéz, Shyam Gopalakrishnan, Mirneza Risvanovic, Zhiheng Wang, Carsten Rahbek, Katharine A. Marske & David Nogués-Bravo
#
# Submitted to Science, 2016
# Code in this file by Michael K. Borregaard
########################################

lat_long_to_long_lat <- function(site)
  unlist(lapply(strsplit(site, "_"), function(el) paste(el[2], el[1], sep = "_")))

#
# STEP 1: FIND THE RICHNESS OF SPECIES IN ANTHROMES AND BIOMES
#
# A function using a nodiv::distrib_data object to count the number of species in a region - the region is given by one of the sitestat attributes of the object
species_per_sitestat <- function(dist_data, site_stat){
    if(!inherits(dist_data, "distrib_data"))
        stop("First argument must be of type distrib_data")
    if(!is.character(site_stat))
        stop("site_stat must be a character vector")
    if(!site_stat %in% sitestat(dist_data))
        stop(paste("No such sitestat in the distrib_data object - possible sitestats are", paste(sitestat(dist_data), collapse = ", ")))

    vals <- sort(na.omit(unique(sitestat(dist_data, site_stat))))
    ret <- sapply(vals, function(val){
      ss <- nodiv:::subsample.distrib_data(dist_data, site = which(sitestat(dist_data, site_stat) == val))
      Nspecies(ss)
    })
    names(ret) <- vals
    ret
}

library(nodiv)

load("terrestrial_mammals.RData")
load("amph_new.RData")

and <- read.csv("Anthro1by1_world.csv", stringsAsFactors = FALSE)[,-(1:4)]
names(and)[1] <- "sites"
and$sites <- lat_long_to_long_lat(and$sites)
termam <- add_sitestat(termam, and)
amph_new <- add_sitestat(amph_new, and)

bio <- read.csv("BIOMES_1BY1.csv", stringsAsFactors = FALSE, na.strings = "#N/A")[,-(1:4)]
names(bio)[1] <- "sites"
bio$sites <- lat_long_to_long_lat(bio[,1])
termam <- add_sitestat(termam, bio)
amph_new <- add_sitestat(amph_new, bio)

latbands <- 10*floor(coords(termam)[,2] / 10)
names(latbands) <- sites(termam)
termam <- add_sitestat(termam, latbands)
latbands <- 10*floor(coords(amph_new)[,2] / 10)
names(latbands) <- sites(amph_new)
amph_new <- add_sitestat(amph_new, latbands)
rm(latbands)


#######################################################
# STEP 2 - CALCULATE SUMMARY VALUES

# Get the values from the genetic data

maketab <- function(name){
  tab <- read.csv(paste(name,"_numbers.csv", sep = ""))
  mam <- species_per_sitestat(termam, name)
  amp <- species_per_sitestat(amph_new, name)
  tab$mammalsdist_richness <- mam[match(tab$site, names(mam))]
  tab$amphibiansdist_richness <- amp[match(tab$site, names(amp))]
  tab$am_prop <- tab$amphibiansrichness / tab$amphibiansdist_richness 
  tab$ma_prop <- tab$mammalsrichness / tab$mammalsdist_richness 
  tab$tot_prop <- tab$totalrichness / (tab$mammalsdist_richness + tab$amphibiansdist_richness)
  tab
}
 

biomGD <- maketab("biomes")
biomGD <- biomGD[!biomGD$site %in% c(0, 99), ]
anthromeGD <- maketab("anthromes")
anthromeGD <- anthromeGD[!anthromeGD$site %in% c(0, 99), ]
latbandGD <- maketab("latbands")

# From here, they were taken to other plotting software. Similar plots can be generated with this function:

plotFig3 <- function(tab){
  
  # A function to add a bar
  addbar <- function(x, y, W, col, horiz = T){
    W <- W/2
    if(horiz)
      rect(xleft=0, ybottom=x-W, xright=y, ytop=x+W, col=col, border = NA)
    else
      rect(xleft=x-W, ybottom=0, xright=x+W, ytop=y, col=col, border = NA)
  }
  
  # A function to plot the right color
  create.cols <- function(vec, col, zlim)
  {
    if(missing(zlim)){
      zlim <- range(vec, na.rm = TRUE)
      if(zlim[1] * zlim[2] < 0)
        zlim <- c(-max(abs(vec), na.rm = TRUE), max(abs(vec), na.rm = TRUE))
    }
    
    if(length(col) == length(unique(vec)))
      return(col[match(vec, sort(unique(vec)))])
    
    vec = vec - zlim[1]
    vec = floor(vec * (length(col)-1)/(zlim[2]- zlim[1]))+1
    return(col[vec])
  }
  
  
  xs <- with(tab, barplot(totalGDval, border = NA, col = "grey", names.arg = site, horiz = T))
  addbar(xs, tab$totalbasepairs/320000000, 0.6, col = create.cols(tab$tot_prop, RColorBrewer::brewer.pal(9,"BuGn")[-1], zlim = c(0,0.25)))
  axis(3, at = 0:4*2/320, labels = 0:4*2)
}



plotFig3(biomGD[rev(order(biomGD$totalGDval)), ])
plotFig3(anthromeGD[rev(order(anthromeGD$site)), ])
plotFig3(latbandGD[order(latbandGD$site), ])


library(Cairo)

Cairo("total_biome.pdf", width = 6, height = 6, units = "in", type = "pdf")
plotFig3(biomGD[rev(order(biomGD$totalGDval)), ])
legend("topright", fill = RColorBrewer::brewer.pal(9,"BuGn")[-1], bty = "n", legend = paste(0:7*3, " - ", 1:8 * 3, "%", sep = ""))
dev.off()

Cairo("total_anthrome.pdf", width = 6, height = 6, units = "in", type = "pdf")
plotFig3(anthromeGD[rev(order(anthromeGD$totalGDval)), ])
legend("topright", fill = RColorBrewer::brewer.pal(9,"BuGn")[-1], bty = "n", legend = paste(0:7*3, " - ", 1:8 * 3, "%", sep = ""))
dev.off()

Cairo("total_latbands.pdf", width = 6, height = 6, units = "in", type = "pdf")
plotFig3(latbandGD[rev(order(latbandGD$totalGDval)), ])
legend("topright", fill = RColorBrewer::brewer.pal(9,"BuGn")[-1], bty = "n", legend = paste(0:7*3, " - ", 1:8 * 3, "%", sep = ""))
dev.off()

