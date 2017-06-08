########################################
# Codes for the paper:
# An Anthropocene Map of Genetic Diversity
#
# Andreia Miraldo, Sen Li, Michael K. Borregaard, Alexander Floréz-Rodriguéz, Shyam Gopalakrishnan, Mirneza Risvanovic, Zhiheng Wang, Carsten Rahbek, Katharine A. Marske & David Nogués-Bravo
#
# Submitted to Science, 2016
# Code in this file by Michael K. Borregaard
########################################


# NOTE:
# The codes in this file rests on the existence of text files detailing the geocoded sequences of any genes (not just cytb) for any species of mammal or amphibian. These text files are:
# * species_grid_amphs
# * species_grid_mammals
# * tot_length_grid_amphs
# * tot_length_grid_mammals
# * corrected_basepairs_gridcell.csv
#
#
# For information on how to generate all these files, please consult the supplementary methods
#


# WD is "/figure2_equalsizegrid"

library(nodiv)
load("../../terrestrial mammals data object.RData")

### Get all the data onto the right format - creating nodiv::distrib_data objects ######
########################################################################################

# read the genetic data species numbers (based on all genes, not just cytb like the other analyses)
sp_amphs <- read.table("species_grid_amphs", stringsAsFactors = F)
names(sp_amphs) <- c("site", "specs")
splist_amphs <- lapply(sp_amphs$specs, function(vec) unlist(strsplit(vec, ',')))

sp_mams <- read.table("species_grid_mammals", stringsAsFactors = F)
names(sp_mams) <- c("site", "specs")
splist_mams <- lapply(sp_mams$specs, function(vec) unlist(strsplit(vec, ',')))


# extract amphibian and mammal species - the format is a little weird
amph_gene_sp <- unique(unlist(splist_amphs))
mam_gene_sp <- unique(unlist(splist_mams))


# matching species of mammals
missing_mam <- mam_gene_sp[!(mam_gene_sp %in% gsub(" ", "_", species(ter)))]
length(missing_mam) #are theses some of the ones that need adding?

cd <- data.frame(sites = ter$coords$sites, data.frame(coords(ter)))
cor_amph <- cd[match(sp_amphs$site + 1, cd$sites),]
cor_mam <- cd[match(sp_mams$site + 1, cd$sites),]


#Nodiv data for amphibian genetic diversity at 4deg scales
gca <- t(sapply(splist_amphs, function(specs){
	ret <- numeric(length(amph_gene_sp))
	ind <- match(specs, amph_gene_sp)
	ret[ind] <- 1
	ret
}))

rownames(gca) <- cor_amph$sites
colnames(gca) <- amph_gene_sp

gca <- gca[!is.na(cor_amph$sites),]
cor_amph <- cor_amph[!is.na(cor_amph$sites),]

amph_gene <- distrib_data(gca, cor_amph)

#Nodiv data for mammalian genetic diversity at 4deg scales
gcm <- t(sapply(splist_mams, function(specs){
	ret <- numeric(length(mam_gene_sp))
	ind <- match(specs, mam_gene_sp)
	ret[ind] <- 1
	ret
}))

rownames(gcm) <- cor_mam$sites
colnames(gcm) <- mam_gene_sp

gcm <- gcm[!is.na(cor_mam$sites),]
cor_mam <- cor_mam[!is.na(cor_mam$sites),]

mam_gene <- distrib_data(gcm, cor_mam)

save(ter, mam_gene, amph_gene, file = "Fig2_nodiv_objects.RData")

##################################################
### Add the data to be plotted to the objects ####
##################################################

load("Fig2_nodiv_objects.RData") #to start from here
library(nodiv)

#Ensure that all is placed in the ter object and sorted correctly
mamgene_richness <- richness(mam_gene)
ampgene_richness <- richness(amph_gene)
ter <- add_sitestat(ter, ampgene_richness)
ter <- add_sitestat(ter, mamgene_richness)

# calculate total values, this detour is needed because of the presence of NA's in cells with 0 species
agr <- sitestat(ter, "ampgene_richness")
agr[is.na(agr)] <- 0
mgr <- sitestat(ter, "mamgene_richness")
mgr[is.na(mgr)] <- 0

allgene_richness <- agr + mgr
allgene_richness[allgene_richness == 0] <- NA

agr <- sitestat(ter, "amph_richness")
agr[is.na(agr)] <- 0
mgr <- sitestat(ter, "mam_richness")
mgr[is.na(mgr)] <- 0

alldist_richness <- agr + mgr
alldist_richness[alldist_richness == 0] <- NA
rm(agr, mgr)

ter <- add_sitestat(ter, allgene_richness)
ter <- add_sitestat(ter, alldist_richness)

# Calculate metrics as proportions of the total
prop_gene_mam <- sitestat(ter, "mamgene_richness") / sitestat(ter, "mam_richness")
prop_gene_amph <- sitestat(ter, "ampgene_richness") / sitestat(ter, "amph_richness")
prop_gene_total <- sitestat(ter, "allgene_richness") / sitestat(ter, "alldist_richness")


ter <- add_sitestat(ter, prop_gene_mam)
ter <- add_sitestat(ter, prop_gene_amph)
ter <- add_sitestat(ter, prop_gene_total)


#Add the numbers of base pairs
numbp_amphs <- read.table("tot_length_grid_amphs", stringsAsFactors = FALSE)
names(numbp_amphs) <- c("site", "bps_amph")
numbp_amphs$site <-  numbp_amphs$site +1 #Because the coordinates are messed up!

ter <- add_sitestat(ter, numbp_amphs)

#Add the numbers of base pairs for mammals
numbp_mams <- read.table("tot_length_grid_mammals", stringsAsFactors = FALSE)
names(numbp_mams) <- c("site", "bps_mam")
numbp_mams$site <-   numbp_mams$site+1 #Because the coordinates are messed up!

# TODO this presently does not work #why not?
ter <- add_sitestat(ter, numbp_mams)



 agr <- sitestat(ter, "bps_amphs")
 agr[is.na(agr)] <- 0
 mgr <- sitestat(ter, "bps_mams")
 mgr[is.na(mgr)] <- 0

 bps_total <- agr + mgr
 bps_total[bps_total == 0] <- NA

 ter <- add_sitestat(ter, bps_total)
 rm(agr, mgr)


save(ter, file = "ready for plotting2.RData")

#################################################
# Do the actual plotting ########################
#################################################

load("Analysis_codes/figure2_equalsizegrid/ready for plotting.RData") # To start from here


for(dataset in c("amph", "mam", "total")){

  bp_name <- paste("bps_", dataset, sep = "")
  sp_name <- paste("prop_gene_", dataset, sep = "")

  basepairs <- log10(sitestat(ter, bp_name))
  sp_prop <- log10(sitestat(ter, sp_name))

  Cairo(paste(dataset, "_Fig2.pdf", sep = ""), width = 10, height = 7, type = "pdf", unit = "in")
  two_color_map(ter, sp_prop, basepairs, 10, gridcol = "lightgrey", gridsites = which(!is.na(sp_prop)), overlay_shape = T, shapeborder = "white", shapefill = "grey", colour = "green_red")
  dev.off()
}



summaryTable <- save_twowaysplot(cap_prop(sitestat(ter, "prop_gene_all")), sitestat(ter, "bps_total"), "all")

summaryTable <- data.frame(summaryTable, save_twowaysplot(log10(sitestat(ter, "prop_gene_all")), log10(sitestat(ter, "bps_total")), "all_log")[,-(1:3)])

summaryTable <- data.frame(summaryTable, save_twowaysplot(cap_prop(sitestat(ter, "prop_gene_amph")), sitestat(ter, "bps_amphs"), "amph")[,-(1:3)])

summaryTable <- data.frame(summaryTable, save_twowaysplot(cap_prop(sitestat(ter, "prop_gene_mamm")), sitestat(ter, "bps_mams"), "mamm")[,-(1:3)])

write.table(summaryTable, file = "summaryTable.txt", sep = "\t", row.names = F)
