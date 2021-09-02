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
library(dismo)
library(reshape2)
library(stringr)
library(phylobase)
library(ape)
library(dplyr)
library(vegan)
library(picante)
library(dendextend)
library(abind)
library(betapart)

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

############################
####LOAD OCCURRENCE DATA####

#loading the gps.pts dataframe
gps<-read.csv("dataframes/gps_9Jun20.csv")

#creating a new column with Species.Name split into Genus and Species columns
new.cols<-colsplit(gps$Species.Name," ",c("Genus","Species"))
gps<-cbind(new.cols, gps[,c(3,4)])

#another way of splitting strings
#parts<-strsplit(as.character(gps$Species.Name), "[\\ ]")
#new.sp<-lapply(parts, "[", 2)
#new.sp<-unlist(new.sp)

#running the models for a subset of genera
genus.sub<-c("Digitipes", "Rhysida", "Ethmostigmus")

#subsetting only genera in genus.sub
gps<-gps[gps$Genus %in% genus.sub,]
gps<-droplevels.data.frame(gps)

#looking at the number of sampling points in each species
table(gps$Species)

#creating a new species column where multiple putative subspecies in Digitipes are collapsed
gps$Species[grep("barnabasi", gps$Species)]<-"barnabasi"
gps$Species[grep("coonoorensis", gps$Species)]<-"coonoorensis"
gps$Species[grep("jonesii", gps$Species)]<-"jonesii"

#looking at the number of sampling points in each species
sp.freq<-as.data.frame(table(gps$Species))
sp.rm<-sp.freq$Var1[which(sp.freq$Freq<3)]

#removing species where there are less than three presence locations
gps<-gps[!(gps$Species %in% sp.rm),]
gps$Species<-factor(gps$Species)

#drop factor levels
gps<-droplevels.data.frame(gps)

#creating a sp object with species names
sp<-levels(gps$Species)
sp2<-substr(sp, 1, 5)

####################################
####VIEW 1 km THRESHOLDED RASTER####
#choosing the model
mods<-"bc2"

#choosing the model selection method
mod.sel<-"rank"

#plot thresholded maps for bc2-rank models
for(i in 1:length(sp)){
#provide the species specific folder name
fold<-paste0("results/revisions/", unique(gps$Genus[gps$Species==sp[i]]), "_", sp[i], "/processed")

#loading the model - the object called model.sel would be loaded into the environment
load(paste0(fold, "/", sp2[i], "_", mods, "_", mod.sel, "_model"))

#load the occurrence locations
loc.sp<-read.csv(paste0("results/revisions/", unique(gps$Genus[gps$Species==sp[i]]), "_", sp[i], "/", sp2[i], "_", mods, "_occ.csv"), header=TRUE)

#in some of the models the model.sel object is saved as model.rank
env<-ls()
if("model.rank" %in% env){
model.sel<-get("model.rank", envir=.GlobalEnv)
rm(model.rank)
}

#find the threshold corresponding to MSS
mss<-model.sel@results[grep("Maximum.training.sensitivity.plus.specificity.Cloglog.threshold", rownames(model.sel@results))]

#load the raster
predict.sel<-raster(paste0(fold, "/", sp2[i], "_", mods, "_", mod.sel, "_map.tif"))

#cropping the prediction raster to the extent of Western Ghats by first creating a template
#using mask directly on sdm.100 does not work as it removes cells with small overlaps with the polygon - which is unusual but I wasn't able to find a solution to this.
sdm<-raster::mask(predict.sel, wg)

#converting occurrence locations into spatial points object
coordinates(loc.sp)<-~LON+LAT
crs(loc.sp)<-wgs

#subsetting locations within Western Ghats alone
loc.sp<-loc.sp[!is.na(over(loc.sp, wg)[,1]),]

#convert the raster layer into a linear vector of raster values
sdm.vec<-as.data.frame(rasterToPoints(sdm))

#rename the last column
colnames(sdm.vec)[3]<-"prob"

#converting values less than the threshold to 0 and greater than or equal to threshold to 1
sdm.vec$prob<-ifelse(sdm.vec$prob<mss, 0, 1)

#use the above dataframe for plotting
plot.sdm<-sdm.vec
plot.sdm$prob<-as.character(plot.sdm$prob)
plot.sdm$prob<-factor(plot.sdm$prob, levels=c("0", "1"))

#plotting the thresholded raster and saving it to disk
#normal scale
my.fill<-scale_fill_manual(name=NULL, labels=c("Absence", "Presence"), values=c("grey", "#31a354"), drop=F)
my.scale<-scale_y_continuous(breaks=seq(8, 24, 2))

#plot for western ghats alone
plot<-ggplot(data=plot.sdm, aes(x=x, y=y)) + geom_raster(aes(fill=prob)) + geom_point(data=as.data.frame(loc.sp), aes(x=LON, y=LAT), shape=1) + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + xlab("Longitude") + ylab("Latitude") + ggtitle(paste0(sp[[i]], " ", mods, " ",  mod.sel)) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_fixed() + my.fill + my.scale

pdf(file=paste0(fold, "/", sp2[i], "_", mods, "_", mod.sel, "_map_1km.pdf"), height=8, width=5, useDingbats=FALSE)
print(plot)
dev.off()
} 

####################################
####CREATE BIODIVERSE INPUT FILE####

#choosing the model
mods<-"bc2"

#choosing the model selection method
mod.sel<-"rank"

#create a list where results from different species will be stored
biodiv<-list()

#opening the folder corresponding to each species and loading the evaluate results from the species model
for(i in 1:length(sp)){

#provide the species specific folder name
fold<-paste0("results/revisions/", unique(gps$Genus[gps$Species==sp[i]]), "_", sp[i], "/processed")

#loading the model - the object called model.sel would be loaded into the environment
load(paste0(fold, "/", sp2[i], "_", mods, "_", mod.sel, "_model"))

#load the occurrence locations
loc.sp<-read.csv(paste0("results/revisions/", unique(gps$Genus[gps$Species==sp[i]]), "_", sp[i], "/", sp2[i], "_", mods, "_occ.csv"), header=TRUE)

#in some of the models the model.sel object is saved as model.rank
env<-ls()
if("model.rank" %in% env){
model.sel<-get("model.rank", envir=.GlobalEnv)
rm(model.rank)
}

#find the threshold corresponding to MSS
mss<-model.sel@results[grep("Maximum.training.sensitivity.plus.specificity.Cloglog.threshold", rownames(model.sel@results))]

#load the raster
predict.sel<-raster(paste0(fold, "/", sp2[i], "_", mods, "_", mod.sel, "_map.tif"))

#cropping the prediction raster to the extent of Western Ghats by first creating a template
#using mask directly on sdm.100 does not work as it removes cells with small overlaps with the polygon - which is unusual but I wasn't able to find a solution to this.
temp<-mask(predict.sel, wg)
temp<-raster::aggregate(temp, fact=100, fun=max, na.rm=TRUE)

#aggregate raster cells to a higher resolution
#0.0083 is 1 km and 0.83 is 100 km
sdm.100<-raster::aggregate(predict.sel, fact=100, fun=max, na.rm=TRUE)
sdm.100<-mask(sdm.100, temp)

#loading the sampling points corresponding to the ith species
coordinates(loc.sp)<-~LON+LAT
crs(loc.sp)<-wgs

#subsetting locations within western ghats alone
loc.sp<-loc.sp[!is.na(over(loc.sp, wg)[,1]),]

#since there is wide variation in the predictions for sada and sp. 1, using only the presence locations and turning all other predictions to 0
if(sp[i] %in% c("sada", "sp. 1")){
sdm.100[!is.na(values(sdm.100))]<-0
}

#find the cell positions of loc.sp and assign those as presence locations
#this is being done since some of the presence locations have been assigned low values of habitat suitability by the model, falling below the threshold and are removed from the thresholded predictions.
pres.cell<-cellFromXY(sdm.100, loc.sp)
sdm.100[pres.cell]<-1

#convert the raster layer into a linear vector of raster values
sdm.vec<-as.data.frame(rasterToPoints(sdm.100))

#rename the last column
colnames(sdm.vec)[3]<-"prob"

#converting values less than the threshold to 0 and greater than or equal to threshold to 1
sdm.vec$prob<-ifelse(sdm.vec$prob<mss, 0, 1)

#use the above dataframe for plotting
plot.sdm<-sdm.vec
plot.sdm$prob<-as.character(plot.sdm$prob)
plot.sdm$prob<-factor(plot.sdm$prob, levels=c("0", "1"))

#plotting the thresholded raster and saving it to disk
#normal scale
my.fill<-scale_fill_manual(name=NULL, labels=c("Absence", "Presence"), values=c("grey", "#31a354"), drop=F)
my.scale<-scale_y_continuous(breaks=seq(8, 24, 2))

#plot for western ghats alone
plot<-ggplot(data=plot.sdm, aes(x=x, y=y)) + geom_raster(aes(fill=prob)) + geom_point(data=as.data.frame(loc.sp), aes(x=LON, y=LAT), shape=1) + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + xlab("Longitude") + ylab("Latitude") + ggtitle(paste0(sp[[i]], " ", mods, " ",  mod.sel)) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_fixed() + my.fill + my.scale
    
pdf(file=paste0(fold, "/", sp2[i], "_", mods, "_", mod.sel, "_map_100km.pdf"), height=8, width=5, useDingbats=FALSE)
print(plot)
dev.off()

#adding the species name to the presence-absence column
colnames(sdm.vec)[3]<-paste0(unique(gps$Genus[gps$Species==sp[i]]), "_", sp[i])

#saving this dataframe to the list
biodiv[[i]]<-sdm.vec  

#remove modEval
rm(model.sel, mss, predict.sel)
}

#There were no presence locations for R. konda in the Western Ghats according to the threshold defined for the Maxent model when the SDM predictions were cropped before the aggregation. However, now since the cropping is done after the aggregation, there is one cell with presence in the Maxent prediction. 

#combining results from different species
#also checking if the x and y coordinates are the same across species
#starting with saving the x and y coordinates of the first
biodiv.mat<-biodiv[[1]][,c(1,2)]

for(i in 1:length(biodiv)){
if(identical(biodiv.mat[,1:2], biodiv[[i]][,1:2])){
biodiv.mat<-cbind(biodiv.mat, biodiv[[i]][,3])
colnames(biodiv.mat)[2+i]<-colnames(biodiv[[i]])[3]
}else{
print(paste0("The coordinates don't match for the ", i, "th species."))
}
}

#change presences for E. coonooranus and E. sahyadrensis
biodiv.mat$Ethmostigmus_coonooranus[biodiv.mat$y>16]<-0

#For E. agasthyamalaensis, removing predictions above 10.5 degrees
biodiv.mat$Ethmostigmus_sahyadrensis[biodiv.mat$y<14]<-0

#give id names to the rows
biodiv.mat$id<-seq(1:nrow(biodiv.mat))
colnames(biodiv.mat)
biodiv.mat<-biodiv.mat[,c(22,1:21)]

#change the name of Rhysida sp 1
colnames(biodiv.mat)[20]<-"Rhysida_sp_1"

#writing this to disk as a .csv file
write.csv(biodiv.mat, file="dataframes/biodiverse/biodiverse_16Mar21.csv", row.names=FALSE)

####################################
####CREATE BIODIVERSE INPUT TREE####
#creating the input tree file for biodiverse
#loading the tree file
cen.tree<-phylo4(read.tree("dataframes/biodiverse/old_files/threegenera.tre"))

#inputing the tip label remap code
remap<-read.csv("dataframes/tiplabels.csv", header=FALSE)
colnames(remap)<-c("tip_name", "sp_name")

#removing extra quotes from tip labels
remap$tip_name<-str_replace_all(remap$tip_name, "\"", "")

#giving the final name
remap$final_name<-c(rep("Digitipes_jonesii", 2), rep("Digitipes_coonoorensis", 2), "Digitipes_jangii", "Digitipes_nudus", rep("Digitipes_barnabasi", 3), "Rhysida_pazhuthara", "Rhysida_crassispina", "Rhysida_konda", "Rhysida_longipes", "Rhysida_sp_2", "Rhysida_lewisi", "Rhysida_sp_1", "Rhysida_trispinosa", "Rhysida_sada", "Rhysida_aspinosa", "Rhysida_ikhalama", "Rhysida_immarginata", "Ethmostigmus_agasthyamalaiensis", "Ethmostigmus_coonooranus", "Ethmostigmus_sahyadrensis", "Ethmostigmus_praveeni", "Ethmostigmus_tristis")

#write remap to disk
write.csv(remap, "dataframes/biodiverse/tiplabels_biodiverse.csv", row.names=FALSE)

#convert from factor to character
remap[,1]<-as.character(remap[,1])
remap[,2]<-as.character(remap[,2])
remap[,3]<-as.character(remap[,3])

#rename the tip labels with sp_name so that tips can be more easily identified and removed
rename.tips<-function(subtree){
for(i in 1:nTips(subtree)) {
tip_label<-labels(subtree)[i]
tip_label<-str_replace_all(tip_label, "\"", "")
new_label<-remap$sp_name[remap$tip_name==tip_label]
labels(subtree)[i]<-new_label
}
return(subtree)
}

#renaming tips on the tree
cen.tree<-rename.tips(cen.tree)

#remove within species edges - keep the longest branch from an intra-specific clade
cen.tree2<-subset(cen.tree, tips.exclude=c("Digitipes jonesii 2", "Digitipes coonoorensis - 1", "Digitipes barnabasi 2", "Digitipes barnabasi 3", "Rhysida sp 2"))
#Jahnavi has already removed Rhysida ikhalama and Rhysida immarginata from the tree when she gave me the tree file

#checking to see if the tree looks okay
plot(cen.tree)

#rename the tip labels with final_name
rename.tips.fin<-function(subtree){
for(i in 1:nTips(subtree)) {
tip_label<-labels(subtree)[i]
new_label<-remap$final_name[remap$sp_name==tip_label]
labels(subtree)[i]<-new_label
}
return(subtree)
}

#renaming the tip labels again so that they match the names from spatial data
cen.tree2<-rename.tips.fin(cen.tree2)
plot(cen.tree2)

#convert tree object to phylo format and save to disk
cen.tree2<-as(cen.tree2, "phylo")

#write tree to disk
write.tree(cen.tree2, file = "dataframes/biodiverse/threegenera_biodiverse.tre")

####