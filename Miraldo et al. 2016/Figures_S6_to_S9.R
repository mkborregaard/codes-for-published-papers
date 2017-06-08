########################################
# Codes for the paper:
# An Anthropocene Map of Genetic Diversity
#
# Andreia Miraldo, Sen Li, Michael K. Borregaard, Alexander Floréz-Rodriguéz, Shyam Gopalakrishnan, Mirneza Risvanovic, Zhiheng Wang, Carsten Rahbek, Katharine A. Marske & David Nogués-Bravo
#
# Submitted to Science, 2016
# Code in this file by Michael K. Borregaard
########################################

#TODO update with the correct figure numbers

library(nodiv)
library(Cairo)

load("terrestrial mammals data object.RData")

colors <- read.table("Analysis_codes/plotcols.txt", header = T)
plotcols <- rgb(colors, max = 255)


# Read in all the different data files and put them into ter
Numbers <- read.csv("Data/equalarea_numbers.csv", stringsAsFactors = FALSE)
Fig1 <- Numbers[,c(1,2,4)]
names(Fig1)[2:3] <- c("Fig1GD", "totalbasepairs")
ter <- add_sitestat(ter, Fig1)

SFnew <- Numbers[,c(1,5,7,11,13)]
ter <- add_sitestat(ter, SFnew)

# Read in the data
SF1 <- read.csv("Data/Supp_Less10_numbers.csv", stringsAsFactors = FALSE)[1:3]
names(SF1)[2:3] <- c("Less10GD", "Less10Richness")
ter <- add_sitestat(ter, SF1)

SF4 <- read.csv("Data/Supp_Less5_numbers.csv", stringsAsFactors = FALSE)[1:3]
names(SF4)[2:3] <- c("Less5GD", "Less5Richness")
ter <- add_sitestat(ter, SF4)

SF5 <- read.csv("Data/Supp_Atleast5_numbers.csv", stringsAsFactors = FALSE)
names(SF5)[2] <- "AtLeast5GD"
ter <- add_sitestat(ter, SF5)

SF2 <- read.csv("Data/Supp_Atleast10_numbers.csv", stringsAsFactors = FALSE)
names(SF2)[2] <- "AtLeast10GD"
ter <- add_sitestat(ter, SF2)

SF9 <- read.csv("Data/Supp_WeightedGD_numbers.csv", stringsAsFactors = FALSE)
names(SF9)[2] <- "WeightedGD"
ter <- add_sitestat(ter, SF9)

SF3_10 <- read.csv("Data/Supp_Rarefy_numbers_10.csv", stringsAsFactors = FALSE)
names(SF3_10)[2:3] <- c("rare_mean_10", "rare_sd_10")
ter <- add_sitestat(ter, SF3_10)


SF3_5 <- read.csv("Data/Supp_Rarefy_numbers_5.csv", stringsAsFactors = FALSE)
names(SF3_5)[2:3] <- c("rare_mean_5", "rare_sd_5")
ter <- add_sitestat(ter, SF3_5)

SF9B <- read.csv("Data/weighted_gd_length.csv", stringsAsFactors = FALSE)[c("cell_michael", "GDbp")]  #This file created by Sen Li and not included in the code other places
names(SF9B) <- c("cell", "GD_Sen_weight")
ter <- add_sitestat(ter, SF9B)


#also read the matrices with individual runs from the rarefactions
AllRarefy10 <- read.csv("Data/All_rarefied_samples_10.csv", header = FALSE)  #Used for creating panels D in Figs S7 and S8 – code not implemented here
rare10 <- AllRarefy10[match(sites(ter), SF3_10$cell),]


AllRarefy5 <- read.csv("Data/All_rarefied_samples_5.csv", header = FALSE)  #Used for creating panels D in Figs S7 and S8 – code not implemented here
rare5 <- AllRarefy5[match(sites(ter), SF3_5$cell),]


# Plot all the maps with quantile color scales (Figs 1, )
for(name in c("Fig1GD","Less10GD", "AtLeast10GD", "AtLeast5GD", "rare_mean_10", "rare_sd_10", "rare_mean_5", "rare_sd_5", "Less5GD", "WeightedGD", "GD_Sen_weight")){
  Cairo(name, type = "pdf", width = 10.5, height = 7, units = "in")
  plotval <- sitestat(ter, name)
  plot_sitestat(ter, plotval, col = plotcols, legend = T, gridcol = "lightgrey", gridsites = which(!is.na(plotval)), overlay_shape = T, colscale = "quant")
  dev.off()
}


# Plot all the maps with equal interval color scale
for(name in c("Less10Richness", "Less5Richness")){
  Cairo(name, type = "pdf", width = 10.5, height = 7, units = "in")
  plotval <- (sitestat(ter,name))
  plotval[abs(plotval) == Inf] <- NA
  plot_sitestat(ter, plotval, col = plotcols, legend = T, gridcol = "lightgrey", gridsites = which(!is.na(plotval)), overlay_shape = T)
  dev.off()
}


# plot the scatterplots

# FIGURE S8d
Cairo("rarefy_10 vs Least10", type = "postscript", width = 8, height = 6, units = "in")
plot(log(sitestat(ter, "AtLeast10GD")), log(sitestat(ter, "rare_mean_10")), cex = sqrt(sitestat(ter, "totalbasepairs"))/300, xlab = "log GD for species with >= 10 sequences", ylab =  "log GD of 10 random sequences", ylim = c(-10, - 2))
apply(log(rare10), 2, points, x = log(sitestat(ter, "AtLeast10GD")), col = "grey",  cex = sqrt(sitestat(ter, "totalbasepairs"))/300)
abline(0,1, col = "red", lwd = 2)
points(log(sitestat(ter, "AtLeast10GD")), log(sitestat(ter, "rare_mean_10")), cex = sqrt(sitestat(ter, "totalbasepairs"))/300)
dev.off()


# FIGURE S9b
Cairo("Weighted vs GD", type = "postscript", width = 6.8, height = 7, units = "in")
plot(log(sitestat(ter, "Fig1GD")), log(sitestat(ter, "WeightedGD")), cex = sqrt(sitestat(ter, "totalbasepairs"))/300, xlab = "Log GD (Fig 1 main text)", ylab =  "Log weighted GD")
abline(0,1, col = "red")
dev.off()


# FIGURE S9b
Cairo("Sen_Weighted vs GD", type = "postscript", width = 6.8, height = 7, units = "in")
plot(log(sitestat(ter, "Fig1GD")), log(sitestat(ter, "GD_Sen_weight")), cex = sqrt(sitestat(ter, "totalbasepairs"))/300, xlab = "Log GD (Fig 1 main text)", ylab =  "Log base-pair-weighted GD")
abline(0,1, col = "red")
dev.off()



# FIGURE S9d
Cairo("rarefy_5 vs Least5", type = "postscript", width = 8, height = 6, units = "in")
plot(log(sitestat(ter, "AtLeast5GD")), log(sitestat(ter, "rare_mean_5")), cex = sqrt(sitestat(ter, "totalbasepairs"))/300, xlab = "log GD for species with >= 5 sequences", ylab =  "log GD of 5 random sequences", ylim = c(-10, - 2))
apply(log(rare5), 2, points, x = log(sitestat(ter, "AtLeast5GD")), col = "grey",  cex = sqrt(sitestat(ter, "totalbasepairs"))/300)
abline(0,1, col = "red", lwd = 2)
points(log(sitestat(ter, "AtLeast5GD")), log(sitestat(ter, "rare_mean_5")), cex = sqrt(sitestat(ter, "totalbasepairs"))/300)
dev.off()


# FIGURE S2
Cairo("mammals co1 vs cytb", type = "postscript", width = 6.8, height = 7, units = "in")
plot(log(SFnew$mammalsGDval), log(SFnew$macoiGDval), cex = 0.9, xlab = "log GD for mammals cytb (Fig 1B)", ylab =  "log GD for mammals co1 (Fig S1)")
abline(0,1, col = "red")
dev.off()



#Spatial regressions
library(SpatialPack)
zerocorrect <- function(var) (var * (sum(!is.na(var)) - 1) +  0.5)/sum(!is.na(var))

spacecor <- function(var1, var2){
  x <- log(zerocorrect(sitestat(ter, var1)))
  y <- log(zerocorrect(sitestat(ter, var2)))
  print(sum(!is.na(x + y)))
  modified.ttest(x, y, coords(ter))
}

spacecor("Fig1GD", "WeightedGD")        #Fig S9D
spacecor("AtLeast10GD", "rare_mean_10") #Fig S7D
spacecor("AtLeast10GD", "Fig1GD")       #Fig 1A vs FigS7A
spacecor("AtLeast5GD", "rare_mean_5")   #Fig S8D
spacecor("AtLeast5GD", "Fig1GD")        #Fig 1A vs FigS8A
spacecor("macoiGDval", "mammalsGDval")  #Fig1B vs FigS1
