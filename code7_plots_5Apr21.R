#export LD_LIBRARY_PATH=/usr/local/lib
#R CMD BATCH filename.r
#setwd("~/Downloads/maxent_analysis")

#set working directory
setwd("/home/bharti/D/PostDoc data/centipede_spatial")

library(raster)
library(rgdal)
library(rgeos)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(stringr)

#setting projection
wgs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#define extent of area of interest
wg.ext<-extent(68,91,8,24)
wg.ext.sp<-as(wg.ext, "SpatialPolygons")
crs(wg.ext.sp)<-wgs

#loading coastline shape file
coastline<-readOGR(dsn='dataframes/sdm_layers/ne_10m_coastline', layer='ne_10m_coastline')

#cropping to the extent defined above
coastline<-crop(coastline, wg.ext)
coast.map<-fortify(coastline)

#loading western ghats shapefile
wg<-readOGR(dsn='dataframes/sdm_layers/wg_vector', layer='wg_boundary_mod_23Jun20')
wg.map<-fortify(wg)

#loading elevation file - I used the GEBCO file for the simulations, but since it also has bathymetry information using the elevation file from WORLDCLIM just for mapping purposes
elev<-raster("dataframes/sdm_layers/wc2.1_30s_elev.tif")
elev<-crop(elev, coastline)
elev.wg<-mask(elev, wg)

###Plot all sampling locations, distinguish points by genus###
#load the sampling locations that were saved to disk to be used in the supplementary material
loc<-read.csv("results/results_figures/supp_table_locations.csv", header=TRUE)
loc.sp<-loc
coordinates(loc.sp)<-~Longitude+Latitude
crs(loc.sp)<-wgs

#turns out ggplot2 cannot plot a raster as such, I need to convert it into a dataframe and then plot it using geom_raster or geom_tile
elev.map<-rasterToPoints(elev)
elev.map<-as.data.frame(elev.map)
colnames(elev.map)[3]<-"value"

elev.wg.map<-rasterToPoints(elev.wg)
elev.wg.map<-as.data.frame(elev.wg.map)
colnames(elev.wg.map)[3]<-"value"

#checking if all locations fall within the extent
which(over(loc.sp, wg.ext.sp)!=1)

#subsetting locations within western ghats alone
loc.sp.wg<-loc.sp[!is.na(over(loc.sp, wg)[,1]),]

#plot data points
my.col<-scale_color_manual(name="Genus", labels=c("Digitipes", "Ethmostigmus", "Rhysida"), values=c("black", "blue", "red"))
my.shp<-scale_shape_manual(name="Genus", labels=c("Digitipes", "Ethmostigmus", "Rhysida"), values=c(1,1,1))
my.fill<-scale_fill_gradientn(name="Elevation (m)", colours=terrain.colors(10))

pdf(file=paste0("results/results_figures/sampling_locations_wg_outline.pdf"), height=5, width=8, useDingbats=FALSE)
ggplot() + geom_raster(data=elev.map, aes(x,y,fill=value)) + geom_point(data=loc, aes(x=Longitude, y=Latitude, color=Genus, shape=Genus)) + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + xlab("Longitude") + ylab("Latitude") + theme_bw() + theme(legend.text=element_text(face=c(rep("italic", 3)))) + my.col + my.shp + my.fill + coord_fixed()
dev.off()

my.scale<-scale_y_continuous(breaks=seq(8, 24, 2))

pdf(file=paste0("results/results_figures/sampling_locations_wg.pdf"), height=8, width=5, useDingbats=FALSE)
ggplot() + geom_raster(data=elev.wg.map, aes(x,y,fill=value)) + geom_point(data=as.data.frame(loc.sp.wg), aes(x=Longitude, y=Latitude, color=Genus, shape=Genus)) + xlab("Longitude") + ylab("Latitude") + theme_bw() + theme(legend.text=element_text(face=c(rep("italic", 3)))) + my.col + my.shp + my.fill + my.scale + coord_fixed()
dev.off()

###Plot the bias map###
#loading bias map and overlay bias locations
bias<-raster("dataframes/biasMap.tif")
bias.map<-as.data.frame(rasterToPoints(bias))
colnames(bias.map)[3]<-"value"

#loading bias background locations
bias.loc<-read.csv("dataframes/biasBg.csv")
colnames(bias.loc)

#plot bias map
my.fill<-scale_fill_gradientn(name="Bias weight", colours=c("white", "red"))

pdf(file=paste0("results/results_figures/bias_map.pdf"),  height=5, width=8, useDingbats=FALSE)
ggplot() + geom_raster(data=bias.map, aes(x,y,fill=value)) + geom_path(data=coast.map, aes(x=long, y=lat, group=group)) + xlab("Longitude") + ylab("Latitude") + theme_bw() + my.fill + coord_fixed()
dev.off()

###Plot variables that are important to the model###
#load the file derived from Table 1
var.imp<-read.csv("results/revisions/results_figures/results_summary_fig3_5Apr21.csv", header=TRUE)
head(var.imp)

#keep only select columns, also choosing only biasBg
var.imp<-var.imp[var.imp$model=="Transferability", c(1,2,11:14)]
head(var.imp)
colnames(var.imp)[c(3:6)]<-c("imp1", "imp1.val", "imp2", "imp2.val")

#look at the frequency of different important variables
vars<-c(as.character(var.imp$imp1), as.character(var.imp$imp2))
vars<-as.data.frame(table(vars))
vars<-vars[order(vars$Freq, decreasing=TRUE),]
vars<-vars[-1,]

#creating a new column which indicates the rank of importance
var.imp1<-var.imp[,-c(5,6)]
var.imp1$imp<-rep("1", times=nrow(var.imp1))
colnames(var.imp1)[3:4]<-c("var", "val")

var.imp2<-var.imp[,-c(3,4)]
var.imp2$imp<-rep("2", times=nrow(var.imp2))
colnames(var.imp2)[3:4]<-c("var", "val")

var.imp1<-rbind(var.imp1, var.imp2)
var.imp1<-var.imp1[complete.cases(var.imp1),]
var.imp1$var<-droplevels(var.imp1$var)

#changing the order of levels in the variable
levels(var.imp1$var)
var.imp1$var<-factor(var.imp1$var, levels=as.character(vars$vars))

#create a new column which has the genus and species names together, first removing extra spaces at the end of the Genus names - I don't know why they are there
var.imp1$Genus<-str_trim(var.imp1$Genus)
var.imp1$sp<-paste0(var.imp1$Genus, " ", var.imp1$Species)

#changing the order of species in the variable
var.imp1$sp<-factor(var.imp1$sp)
levels(var.imp1$sp)
var.imp1$sp<-factor(var.imp1$sp, levels=rev(levels(var.imp1$sp)))

#plot
my.size<-scale_size(name="Percentage \nImportance", limits=c(1,100), range=c(1, 6))

pdf(file=paste0("results/revisions/results_figures/variable_importance.pdf"), height=8, width=11, useDingbats=FALSE)
ggplot(data=var.imp1, aes(x=var, y=sp)) + geom_point(aes(size=val), col="#e31a1c", alpha=0.7) + xlab("Environmental variables") + ylab("Species")  + my.size + theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1), axis.text.y=element_text(face="italic"))
dev.off()

####