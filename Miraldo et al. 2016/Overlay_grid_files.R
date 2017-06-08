library(rgdal)

grid <- readOGR("Equal area grid", "fishnet4")
proj4string(grid)
#unfortunately this shapefile does not have projection information. According to this website (https://gis.stackexchange.com/questions/138287/what-is-the-proj4-for-world-behrmann-54017) the correct proj4string to use is this:
proj4string(grid) <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs")
# They do not have IDs either, so I will assign one
grid$Id <- 1:nrow(grid@data)


# Go through all of the data files and put the coords files in their right place
for(folder in c("amphibians_cytb", "mammals_coi", "mammals_cytb")){
  outputfolder <- file.path("Data", "biome_anthrome_latband", folder, paste(substr(folder, nchar(folder) - 3, nchar(folder)), "equalarea_matlab", sep = "_"))
  dir.create(outputfolder)
  files <- list.files(file.path("Data", folder), pattern = ".coords")
  for(file in files){
    dat1 <- read.table(file.path("Data", folder, file))
    dat <- spTransform(SpatialPoints(dat1[,c(2,1)], CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")), CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs"))
    dat1$equ <- over(dat, grid)$Id
    write.table(dat1, file.path(outputfolder, file), sep = "\t", col.names = F, row.names = F)
  }
}

