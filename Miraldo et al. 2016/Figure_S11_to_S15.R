######## Figure S11 #######

# We provide the 2 files with range size data for all amphibians and mammals species that we retrieved from IUCN (amphibian_ranges.csv and mammals_ranges.csv)
# Data to plot Fig11 is given in files Fig11_data_amphs.csv and Fig11_data_mammals.csv
# the files have for each species the total number of sequences, range size (km2), and range size as in terms of percentile group (used for Fig.S12)

library(ggplot2)
Fig11_amphs<-read.csv("Fig11_data_amphs.csv")
Fig11_mam<-read.csv("Fig11_data_mammals.csv")

ggplot(Fig11_amphs, aes(x=log(RangeSize_km2), y=log(nseqs)))+geom_point(size=2, shape=1)+geom_smooth(method=lm, aes(group=1),colour="blue")+theme_linedraw()+ theme(panel.background = element_rect(fill="NA"))+theme(panel.grid.major=element_line(colour=NA)) + theme(panel.grid.minor=element_line(colour=NA))+ theme(axis.text=element_text(size=12), axis.title.y=element_text(vjust=1, size=12), axis.title.x=element_text(vjust=-1, size=12))+labs(x="log(Range size)", y="log(Number of sequences)")
ggplot(Fig11_mam, aes(x=log(RangeSize_km2), y=log(nseqs)))+geom_point(size=2, shape=1)+geom_smooth(method=lm, aes(group=1),colour="blue")+theme_linedraw()+ theme(panel.background = element_rect(fill="NA"))+theme(panel.grid.major=element_line(colour=NA)) + theme(panel.grid.minor=element_line(colour=NA))+ theme(axis.text=element_text(size=12), axis.title.y=element_text(vjust=1, size=12), axis.title.x=element_text(vjust=-1, size=12))+labs(x="log(Range size)", y="log(Number of sequences)")

#get the correlation values and significance
cor.test(log(Fig11_amphs$nseqs), log(Fig11_amphs$RangeSize_km2), method="pearson")
cor.test(log(Fig11_mam$nseqs), log(Fig11_mam$RangeSize_km2), method="pearson")

######## Figure S12 #######

Fig12_mam<-read.csv("Fig12_data_mam_persite.csv")
Fig12_amphs<-read.csv("Fig12_data_amphs_persite.csv")

#set variable RangeSize_percentile as factor
Fig12_mam$RangeSize_percentile<-as.factor(Fig12_mam$RangeSize_percentile)
Fig12_amphs$RangeSize_percentile<-as.factor(Fig12_amphs$RangeSize_percentile)

#subsetting the data to include only percentile 10 and 90 of range sizes
Fig12_mam_percentile<- Fig12_mam[Fig12_mam$RangeSize_percentile %in% c("perc_10", "perc_90"), ]
Fig12_amphs_percentile<- Fig12_amphs[Fig12_amphs$RangeSize_percentile %in% c("perc_10", "perc_90"), ]

#plotting the data as shown in FigS12 A and B
ggplot(Fig12_mam_percentile, aes(x=RangeSize_percentile, y=PiHat))+ theme_linedraw()+ theme(panel.background = element_rect(fill="NA"))+theme(panel.grid.major=element_line(colour=NA)) + theme(panel.grid.minor=element_line(colour=NA))+ stat_boxplot(geom="errorbar", colour="black", size=0.50)+ geom_boxplot (outlier.size=2, fill="gray70", alpha=1, outlier.colour="black") + stat_summary(fun.y=mean, geom="point", shape=16, color="red", size=3) + theme(axis.text=element_text(size=12), axis.title.y=element_text(vjust=1, size=0), axis.title.x=element_text(vjust=-1, size=0))+ scale_x_discrete(breaks=c("perc_10", "perc_90"),labels=c("Small range species", "Large range species"))
ggplot(Fig12_amphs_percentile, aes(x=RangeSize_percentile, y=PiHat))+ theme_linedraw()+ theme(panel.background = element_rect(fill="NA"))+theme(panel.grid.major=element_line(colour=NA)) + theme(panel.grid.minor=element_line(colour=NA))+ stat_boxplot(geom="errorbar", colour="black", size=0.50)+ geom_boxplot (outlier.size=2, fill="gray70", alpha=1, outlier.colour="black") + stat_summary(fun.y=mean, geom="point", shape=16, color="red", size=3) + theme(axis.text=element_text(size=12), axis.title.y=element_text(vjust=1, size=0), axis.title.x=element_text(vjust=-1, size=0))+ scale_x_discrete(breaks=c("perc_10", "perc_90"),labels=c("Small range species", "Large range species"))

#Kruskal-Wallis test
kruskal.test(PiHat~RangeSize_percentile, data=Fig12_mam_percentile)
kruskal.test(PiHat~RangeSize_percentile, data=Fig12_amphs_percentile)

######## Figure S13 #######
library(rgdal)
library(nodiv)
library(ggplot2)
library(betareg)

load("terrestrial mammals data object.RData") # this object is used to store all variables for plotting, not just the mammals

Fig1 <- read.csv("Data/equalarea_numbers.csv", stringsAsFactors = FALSE)
ter <- add_sitestat(ter, Fig1)

proj4string(ter$coords) <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs")  #TODO I need to edit nodiv so that the projections are always maintained
latlongpoints <- data.frame(coordinates(spTransform(ter$coords, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "))))
names(latlongpoints) <- c("Lon", "Lat")


library(Cairo)
lat <- latlongpoints$Lat
for(var in c("total", "mammals", "amphibians", "macoi")){
  GD <- sitestat(ter, paste(var, "GDval", sep = ""))
  mod <- lm(GD ~ lat + I(lat^2))
  #print(summary(mod))


  mm <- data.frame(lat, GD)
  mm <- mm[!is.na(mm$GD), ]
  Cairo(paste(var, "parabolic", sep = "_"), type = "postscript", width = 8, height = 6, units = "in")
  p <- ggplot(mm, aes(mm$lat, mm$GD)) + geom_point() + geom_smooth(col = "red", method = "lm", formula = y ~ x + I(x^2) ) + theme_classic() + xlab("Latitude") + ylab("GD of grid cell") + ylim(0, 0.1)
  print(p)
  dev.off()

  mm$abslat <- abs(mm$lat)
  Cairo(paste(var, "absolute", sep = "_"), type = "postscript", width = 8, height = 6, units = "in")
  p <- ggplot(mm, aes(abslat, GD)) + geom_point() + geom_smooth(col = "red", method = "lm", formula = y ~ x) + theme_classic() + xlab("Absolute latitude") + ylab("GD of grid cell") #+ ylim(0, 0.1)
  print(p)
  dev.off()
  mod2 <- lm(GD ~ abslat, data = mm)
  #print(summary(mod2))


  # these plots used for the supplementary materials
  mm$betaGD <-(mm$GD * (nrow(mm) - 1) +  0.5)/nrow(mm)
  Cairo(paste(var, "absolute beta", sep = "_"), type = "postscript", width = 8, height = 6, units = "in")
  mod3 <- betareg(betaGD ~ abslat, data = mm)
  p <- ggplot(mm, aes(abslat, GD)) + geom_point() + geom_line(col = "red", aes(y = predict(mod3), x = mm$abslat)) + theme_classic() + xlab("Absolute latitude") + ylab("GD of grid cell")
  print(p)
  dev.off()
  cat(toupper(paste("\n\n", var, " beta regression, absolute latitude\n", sep = "")))
  print(summary(mod3))
}

######## Figure S14 #######

data_amphs <- read.csv ("FigS14_data_amphs.csv")
data_mam <- read.csv ("FigS14_data_mam.csv")
data_mam_coi<- read.csv ("FigS14_data_mam_coi.csv")

# Transform PiHat values using the formula (y · (n − 1) + 0.5)/n, where n is the sample size (Smithson and Verkuilen 2006)
data_amphs$PiHat_transformed <- (data_amphs$PiHat * (length(data_amphs$PiHat) - 1) + 0.5) / length(data_amphs$PiHat)
data_mam$PiHat_transformed <- (data_mam$PiHat * (length(data_mam$PiHat) - 1) + 0.5) / length(data_mam$PiHat)
data_mam_coi$PiHat_transformed <- (data_mam_coi$PiHat * (length(data_mam_coi$PiHat) - 1) + 0.5) / length(data_mam_coi$PiHat)
data_combined$PiHat_transformed <- (data_combined$PiHat * (length(data_combined$PiHat) - 1) + 0.5) / length(data_combined$PiHat)

##Kruskal Wallis
kt.amph=kruskal.test(PiHat_transformed~factor(Anthrome), data=data_amphs)
kt.amph
kt.mam=kruskal.test(PiHat_transformed~factor(Anthrome), data=data_mam)
kt.mam
kt.mamcoi=kruskal.test(PiHat_transformed~factor(Anthrome), data=data_mam_coi)
kt.mamcoi

#Jonckheere.test
library(clinfun)
GD_amphs<-data_amphs$PiHat_transformed
n<-c(17,52,114,135,173,24)
grp<-as.ordered(factor(rep(1:length(n),n)))
jonckheere.test(GD_amphs,grp,alternative="increasing", nperm=1000)

GD_mamcoi<-data_mam_coi$PiHat_transformed
n1<-c(21,127,300,242,367,252)
grp<-as.ordered(factor(rep(1:length(n1),n1)))
jonckheere.test(GD_mamcoi,grp,alternative="increasing", nperm=1000)

GD_mam<-data_mam$PiHat_transformed
n1<-c(44,109,248,310,293,117)
grp<-as.ordered(factor(rep(1:length(n1),n1)))
jonckheere.test(GD_mam,grp,alternative="increasing", nperm=1000)


######## Figure 15 #######
library(ggplot2)

#plotting number of cytb sequences per species (combined) per grid cell. Note that the full range graph y limit goes up to ~400 sequences. It was cut at 25 to aid visualization, removing therefore all ouliers
cytb_seqs<-read.csv("FigS15_GrainSize_data_cytb.csv")
cytb_seqs$GrainSize<-as.factor(cytb_seqs$GrainSize)
cytb_seqs[, 3][cytb_seqs[, 3] == 1] <- NA #removing species with only 1 sequence, as these are excluded from analyses

graph<-ggplot(cytb_seqs, aes(x=GrainSize, y=num_seqs)) + theme_linedraw()+ theme(panel.background = element_rect(fill="NA")) +theme(panel.grid.major=element_line(colour=NA)) + theme(panel.grid.minor=element_line(colour=NA))
FigS15<-graph + coord_cartesian(ylim=c(0,25)) + stat_boxplot(geom="errorbar", colour="black", size=0.50) + geom_boxplot (coef=0, outlier.size=0, fill="#cccccc", alpha=1, outlier.colour=NA) + labs(x="Grid cell grain (km)", y="Number of sequences per species in each grid cell") + theme(axis.text=element_text(size=8), axis.title.y=element_text(vjust=1, size=9), axis.title.x=element_text(vjust=-1, size=9))
FigS15
