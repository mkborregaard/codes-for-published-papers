library(rgeos)    # for gIntersection(...), etc.
library(rgdal)    # for readOGR(...)


# FUNCTIONS FOR CALCULATING AREA
get.ar <- function(gr,sp){
  int <- gIntersection(sp,grid[grid$Id==gr,],byid=TRUE)
  if (length(int)==0) return(0)
  1 # gArea(int)
}

get.area <- function(spec, data){
  cat(paste("species ", i <<- i + 1, " of ", maxx,":\t",spec, "\n", sep = ""))
  sppol <- data[data$binomial == spec,]
  sppol <- gBuffer(sppol, width = 0)
  sppol <- gUnaryUnion(sppol)
  res <- numeric(length(Id))
  bb <- bbox(sppol) + rbind(c(-1000,1000), c(-1000, 1000))
  rels <- which(cd[,1] > bb[1,1] & cd[,1] < bb[1,2] & cd[,2] > bb[2,1] & cd[,2] < bb[2,2])
  maxx <- nrow(sppol)
  tmp <- sapply(Id[rels], get.ar, sp = sppol)
  res[rels] <- tmp
  res
}


# shapefile grid
grid <- readOGR("Equal area grid", "fishnet4")
proj4string(grid) <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs")
# They do not have IDs either, so I will assign one
grid$Id <- 1:nrow(grid@data)
Id <- grid$Id
cd <- coordinates(grid)

# read the shapefile
spec <- readOGR("Data/TERRESTRIAL_MAMMALS", "TERRESTRIAL_MAMMALS")

# select polygons according to the listed criteria
permitted <- which(spec@data$origin %in% c(1, 2, 5) & spec@data$presence %in% c(1, 2, 4, 5))
dataset <- spec[permitted, ]
specnames <- unique(dataset$binomial)
dataset <- spTransform(dataset, CRS(proj4string(grid)))
rangesizes <- gArea(dataset, byid = T)
ranges <- tapply(rangesizes, dataset$binomial, sum)
write.csv(as.data.frame(ranges), "mammal_ranges.csv")

i <- 0
maxx <- length(specnames)
# Area of each mountain in each grid cell
specareas <- sapply(specnames, get.area, data = dataset) #This here needs to be divided with the total area (allareas) - and then we will recalculater areas based on the full global area function.
colnames(specareas) <- gsub(" ", "_", specnames)



# read the shapefile
am <- readOGR("Data/AMPHIBIANS", "AMPHIBIANS")

# select polygons according to the listed criteria
permitted <- which(am@data$origin %in% c(1, 2, 5) & am@data$presence %in% c(1, 2, 4, 5))
dataam <- am[permitted, ]
amnames <- unique(dataam$binomial)
dataam<- spTransform(dataam, CRS(proj4string(grid)))

rangesizes <- gArea(dataam, byid = T)
ranges <- tapply(rangesizes, dataam$binomial, sum)
write.csv(as.data.frame(ranges), "amphibians_ranges.csv")

amareas <- sapply(amnames, get.area, data = dataam) #This here needs to be divided with the total area (allareas) - and then we will recalculater areas based on the full global area function.

colnames(amareas) <-gsub(" ", "_", amnames)
save(amareas, specareas, file = "New IUCN 4x4 equal area.RDAta")


## NEW Create nodiv objects (already at the right scales)
library(nodiv)
coords <- data.frame(cd)
names(coords) <- c("Lon", "Lat")
coords$sites <- Id
ter_Behr <- distrib_data(specareas, coords, type = "grid" )
amp_Behr <- distrib_data(amareas, coords, type = "grid" )

mam_richness <- richness(ter_Behr)
ter_Behr <- add_sitestat(ter_Behr, mam_richness)
amph_richness <- richness(amp_Behr)
ter_Behr <- add_sitestat(ter_Behr, amph_richness)

shape <- readOGR("Data/Shapefile","World_coast")
shape <- spTransform(shape, CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs"))
ter_Behr <- add_shape(ter_Behr, shape)
save(ter_Behr, amp_Behr, file = "nodiv objects 4 deg Behrmann.RData")
ter <- ter_Behr
save(ter, file = "terrestrial mammals data object.RData")
