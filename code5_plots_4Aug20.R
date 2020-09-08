###This code plots Figure 1 and some other figures that have been placed in the supplementary material

#set working directory - this should lead to the location where the /results, /dataframes and /R scripts folders are present
setwd("/home/bharti/D/PostDoc data/spatial_data")

library(raster)
library(rgdal)
library(rgeos)
library(reshape2)
library(ggplot2)
library(ggpubr)

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

###Figure 1: Plot all sampling locations, distinguish points by genus###
#load the sampling locations that were saved to disk to be used in the supplementary material
loc<-read.csv("dataframes/supp_table_locations.csv", header=TRUE)
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
which(over(loc.sp, wg.ext.sp)!=1) #155 does not fall within the extent - checking the row details from the Excel sheet - the Bondla, Goa point for Rhysida longipes was incorrect - the Longitude was given as 1.1053 instead of 74.1053 #this is fixed from the source onwards now

#subsetting locations within western ghats alone
loc.sp.wg<-loc.sp[!is.na(over(loc.sp, wg)[,1]),]

#plot data points
my.col<-scale_color_manual(name="Genus", labels=c("Digitipes", "Ethmostigmus", "Rhysida"), values=c("black", "blue", "red"))
my.shp<-scale_shape_manual(name="Genus", labels=c("Digitipes", "Ethmostigmus", "Rhysida"), values=c(1,1,1))
my.fill<-scale_fill_gradientn(name="Elevation (m)", colours=terrain.colors(10))

pdf(file=paste0("results/results_figures/fig1_sampling_locations.pdf"),  height=5, width=8, useDingbats=FALSE)
ggplot() + geom_raster(data=elev.map, aes(x,y,fill=value)) + geom_point(data=loc, aes(x=Longitude, y=Latitude, color=Genus, shape=Genus)) + xlab("Longitude") + ylab("Latitude") + theme_bw() + theme(legend.text=element_text(face=c(rep("italic", 3)))) + my.col + my.shp + my.fill + coord_fixed()
dev.off()

pdf(file=paste0("results/results_figures/fig1_sampling_locations_wg_outline.pdf"), height=5, width=8, useDingbats=FALSE)
ggplot() + geom_raster(data=elev.map, aes(x,y,fill=value)) + geom_point(data=loc, aes(x=Longitude, y=Latitude, color=Genus, shape=Genus)) + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + xlab("Longitude") + ylab("Latitude") + theme_bw() + theme(legend.text=element_text(face=c(rep("italic", 3)))) + my.col + my.shp + my.fill + coord_fixed()
dev.off()

my.scale<-scale_y_continuous(breaks=seq(8, 24, 2))

pdf(file=paste0("results/results_figures/fig1_sampling_locations_wg.pdf"), height=8, width=5, useDingbats=FALSE)
ggplot() + geom_raster(data=elev.wg.map, aes(x,y,fill=value)) + geom_point(data=as.data.frame(loc.sp.wg), aes(x=Longitude, y=Latitude, color=Genus, shape=Genus)) + xlab("Longitude") + ylab("Latitude") + theme_bw() + theme(legend.text=element_text(face=c(rep("italic", 3)))) + my.col + my.shp + my.fill + my.scale + coord_fixed()
dev.off()

###Figure 2: Plot the bias map###
#loading bias map and overlay bias locations
bias<-raster("dataframes/biasMap.tif")
bias.map<-as.data.frame(rasterToPoints(bias))
colnames(bias.map)[3]<-"value"

#loading bias background locations
bias.loc<-read.csv("dataframes/biasBg.csv")
colnames(bias.loc)

#plot bias map
my.fill<-scale_fill_gradientn(name="Bias weight", colours=c("white", "red"))

pdf(file=paste0("results/results_figures/fig2_bias_map.pdf"),  height=5, width=8, useDingbats=FALSE)
ggplot() + geom_raster(data=bias.map, aes(x,y,fill=value)) + geom_path(data=coast.map, aes(x=long, y=lat, group=group)) + xlab("Longitude") + ylab("Latitude") + theme_bw() + my.fill + coord_fixed()
dev.off()

###Figure 3: Plot variables that are important to the model###
#load the file derived from Table 1
var.imp<-read.csv("results/results_figures/results_summary_fig3_10Jun20.csv", header=TRUE)
head(var.imp)

#keep only select columns, also choosing only biasBg
var.imp<-var.imp[var.imp$model=="randomBg",c(1,2,8:11)]
head(var.imp)
colnames(var.imp)[c(3:6)]<-c("imp1", "imp1.val", "imp2", "imp2.val")

#look at the frequency of different important variables
vars<-c(as.character(var.imp$imp1), as.character(var.imp$imp2))
vars<-as.data.frame(table(vars))
vars<-vars[order(vars$Freq, decreasing=TRUE),]
vars<-vars[-9,]

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

pdf(file=paste0("results/results_figures/fig3_variable_importance.pdf"), height=5, width=8, useDingbats=FALSE)
ggplot(data=var.imp1, aes(x=var, y=sp)) + geom_point(aes(size=val), col="#e31a1c", alpha=0.7) + xlab("Environmental variables") + ylab("Species")  + my.size + theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1), axis.text.y=element_text(face="italic"))
dev.off()

###Figure 4: Plot endemism and diversity measures###
div.names<-c("SR", "WE", "PD", "PE")
div<-list()
map.wg<-list()
log.map<-list()
bg.names<-c("Random", "Bias corrected")
h<-2

#load the data files 
for(i in 1:length(div.names)){
div[[i]]<-raster(paste0("results/", div.names[i], ".tif"))
ras<-div[[i]]
ras<-as.data.frame(rasterToPoints(ras))
colnames(ras)<-c("Longitude", "Latitude", "div")

if(div.names[i] %in% c("SR", "PD")){
  
#normal scale
my.fill1<-scale_fill_gradientn(name=div.names[i], colours=rev(c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')))

my.scale<-scale_y_continuous(breaks=seq(8, 24, 2))

#plot for western ghats alone
map.wg[[i]]<-ggplot(data=ras, aes(x=Longitude, y=Latitude)) + geom_tile(aes(fill=div)) + my.fill1 + ggtitle(paste0(bg.names[h], " - ", div.names[i])) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.scale + coord_fixed()

} else{
  
#log scale
my.fill2<-scale_fill_gradientn(name=paste0("log(",div.names[i], ")"), colours=rev(c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')))

#plot for western ghats alone
map.wg[[i]]<-ggplot(data=ras, aes(x=Longitude, y=Latitude)) + geom_tile(aes(fill=log(div))) + my.fill2 + ggtitle(paste0(bg.names[h], " - log(", div.names[i], ")")) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.scale + coord_fixed()
}
}

#plotting the graphs individually
pdf(file=paste0("results/results_figures/fig4a_SR.pdf"), height=8, width=5, useDingbats=FALSE)
map.wg[[1]]
dev.off()

jpeg(file=paste0("results/results_figures/fig4a_SR.jpeg"), height=8, width=5, res=300, units="in")
map.wg[[1]]
dev.off()

pdf(file=paste0("results/results_figures/fig4b_PD.pdf"), height=8, width=5, useDingbats=FALSE)
map.wg[[3]]
dev.off()

jpeg(file=paste0("results/results_figures/fig4b_PD.jpeg"), height=8, width=5, res=300, units="in")
map.wg[[3]]
dev.off()

pdf(file=paste0("results/results_figures/fig4c_WE.pdf"), height=8, width=5, useDingbats=FALSE)
map.wg[[2]]
dev.off()

jpeg(file=paste0("results/results_figures/fig4c_WE.jpeg"), height=8, width=5, res=300, units="in")
map.wg[[2]]
dev.off()

pdf(file=paste0("results/results_figures/fig4b_PE.pdf"), height=8, width=5, useDingbats=FALSE)
map.wg[[4]]
dev.off()

jpeg(file=paste0("results/results_figures/fig4d_PE.jpeg"), height=8, width=5, res=300, units="in")
map.wg[[4]]
dev.off()

#arranging the above graphs
pdf(file=paste0("results/results_figures/fig4_diversity_measures.pdf"),  height=15, width=12, useDingbats=FALSE)
ggarrange(map.wg[[1]], map.wg[[3]], map.wg[[2]], map.wg[[4]], ncol=2, nrow=2)
dev.off()

jpeg(file=paste0("results/results_figures/fig4_diversity_measures.jpeg"),  height=15, width=12, res=300, units="in")
ggarrange(map.wg[[1]], map.wg[[3]], map.wg[[2]], map.wg[[4]], ncol=2, nrow=2)
dev.off()

###     

###LOOSE CODE###
# div[[i]]<-read.csv(paste0("results/", div.names[i], ".csv"), header=TRUE)
# ras<-div[[i]]

# #creating another raster following the WG shape file
# temp<-ras
# coordinates(temp)<-~Longitude + Latitude
# crs(temp)<-wgs
# ras.wg<-temp[complete.cases(over(temp, wg)),]
# ras.wg<-as.data.frame(ras.wg)

#ras.wg<-rasterFromXYZ(ras)
#crs(ras.wg)<-wgs
#ras.wg<-mask(ras.wg, wg)
#ras.wg2<-as.data.frame(rasterToPoints(ras.wg))
#colnames(ras.wg2)<-c("Longitude", "Latitude", "div")
