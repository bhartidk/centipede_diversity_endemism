#export LD_LIBRARY_PATH=/usr/local/lib
#R CMD BATCH filename.r
setwd("~/Downloads/maxent_analysis")

#set working directory
#setwd("/home/bharti/D/PostDoc data/centipede_spatial")

#load necessary libraries
library(raster)
library(sp)
library(rgdal)
library(maps)
library(rgeos)
library(dismo)
library(dplyr)
library(Hmisc)
library(ggplot2)
library(devtools)
library(digest)
library(rJava)
library(geosphere)
library(stringr)
library(ncdf4)
library(sf)
library(reshape2)
library(ENMeval)
library(parallel)
library(devtools)
library(MASS)

#setting projection
wgs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#define extent of area of interest
#wg.ext<-extent(72.5,88,8,22)
wg.ext<-extent(68,91,8,24)

#loading coastline shape file
coastline<-readOGR(dsn='dataframes/sdm_layers/ne_10m_coastline', layer='ne_10m_coastline')
pi.coast<-crop(coastline, wg.ext)

#load Western Ghats shapefile
wg<-readOGR(dsn='dataframes/sdm_layers/wg_vector', layer='wg_boundary')

##PREPARING INPUT DATA##
#load relevant input files from WorldClim - loading all files and loading select files which I had discussed with JJ which make sense for centipedes

#WORLDCLIM 1km layers - downloaded from worldclim.org on 26th March 2020
wc<-list.files(path="dataframes/sdm_layers/wc2.1_30s_bio", pattern="*.tif")

#The file is input, cropped to the required extent and put in the stack. 
wc_layers<-stack()
for(i in 1:length(wc)){
wc.raster<-raster(paste0("dataframes/sdm_layers/wc2.1_30s_bio/", wc[i]))
wc.raster<-crop(wc.raster, wg.ext)
wc_layers<-stack(wc_layers, wc.raster)
}

#Giving the layers meaningful names
names(wc_layers)
names(wc_layers)<-c("ann_mean_temp", "mean_temp_warmQ", "mean_temp_coldQ", "ann_prec", "prec_wet_month", "prec_dry_month", "prec_season", "prec_wetQ", "prec_dryQ", "prec_warmQ", "prec_coldQ", "mean_diur_range", "isotherm", "temp_season", "max_temp", "min_temp", "ann_range_temp", "mean_temp_wetQ", "mean_temp_dryQ")

#Setting min and max
wc_layers<-setMinMax(wc_layers)

#input soil layers
soil.stype<-raster("dataframes/sdm_layers/soil_stype.tif", RAT=TRUE)

#input elevation layer (from gebco)
elev<-raster("dataframes/sdm_layers/wc2.1_30s_elev.tif")
elev<-crop(elev, wg.ext)
names(elev)<-"elevation"

#creating a new stack called predictors
preds<-stack(wc_layers, elev, soil.stype)

#creating a substack of rasters which uses only the primary predictors from WC (as in BC2 of Low et al., 2020)
#BC2 includes mean temperature (bio1), maximum temperature of the warmest month (bio5), minimum temperature of the coldest month (bio6), annual precipitation (bio12), precipitation of the wettest month (bio13) and precipitation of the driest month (bio14) - a total of 6 variables
preds.bc2<-preds[[c(1,15,16,4,5,6,20,21)]]
names(preds.bc2)

#creating another substack of rasters where the predictors are chosen based on the ecology of centipedes
preds.eco<-preds[[(c(1,4,7,9,16,20,21))]]

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

#creating a list where each element is a species-specific dataframe
ind.occ<-list()
for(i in 1:length(levels(gps$Species))){
ind.occ[[i]]<-subset(gps, gps$Species %in% sp[i], select=c(4,3))
}

#input the background data
biasBg<-read.csv("dataframes/biasBg.csv")

##MAXENT MODELS WITH CENTIPEDE DATA##
#running a loop perform the model exercise over species and multiple background locations
for(i in 1:length(ind.occ)){
#create a directory for the species in question
dir.create(paste0("results/revisions/", unique(gps$Genus[gps$Species==sp[i]]), "_", sp[i]))
fold<-paste0("results/revisions/", unique(gps$Genus[gps$Species==sp[i]]), "_", sp[i])

#taking in occurrence values by species
sp.occ<-as.data.frame(ind.occ[i])

#removing multiple sampling locations in a grid cell
source('R_scripts/Eliminate Points in Same Cell of a Raster.r')
sp.occ<-elimCellDups(sp.occ, preds[[1]], longLatFields=c('Longitude', 'Latitude'), priority=NULL)
nrow(sp.occ) 

#if the number of occurrences is between 10-79 I would like to test models with different feature classes (LQH, LQ, QH, L, H). Not using LP - because L is an instance of H. While for number of occurrences less than 10, by default Maxent only uses LC and I would like to do that too. However, hinge features can be described as instances of linear and threshold features, so using that too. Also, based on https://groups.google.com/g/maxent/c/vAB5G3Eo1ho, looks like I don't need to add C (for categorical) to the feature class argument -
fc.comb<-c("LQH", "LQ", "QH", "L", "Q", "H")

#creating a string of regularization multiplier values to use in the model
reg.mult<-seq(0.5,5,by=0.5)

#selecting the columns which have the coordinates of the background sites
spBg.occ<-biasBg[,c(1,2)]

#running and comparing the different models
if(nrow(sp.occ)>19){
eval.preds<-ENMevaluate(occ=sp.occ, env=preds, bg.coords=spBg.occ, method='block', categoricals="soil_stype", RMvalues=reg.mult, fc=fc.comb, algorithm='maxent.jar', overlap=TRUE, bin.output=TRUE, clamp=TRUE, rasterPreds=TRUE, parallel=TRUE, numCores=7)

eval.preds.bc2<-ENMevaluate(occ=sp.occ, env=preds.bc2, bg.coords=spBg.occ, method='block', categoricals="soil_stype", RMvalues=reg.mult, fc=fc.comb, algorithm='maxent.jar', overlap=TRUE, bin.output=TRUE, clamp=TRUE, rasterPreds=TRUE, parallel=TRUE, numCores=7)

eval.preds.eco<-ENMevaluate(occ=sp.occ, env=preds.eco, bg.coords=spBg.occ, method='block', categoricals="soil_stype", RMvalues=reg.mult, fc=fc.comb, algorithm='maxent.jar', overlap=TRUE, bin.output=TRUE, clamp=TRUE, rasterPreds=TRUE, parallel=TRUE, numCores=7)

}else{

eval.preds<-ENMevaluate(occ=sp.occ, env=preds, bg.coords=spBg.occ, method='jackknife', categoricals="soil_stype", RMvalues=reg.mult, fc=fc.comb, algorithm='maxent.jar', overlap=TRUE, bin.output=TRUE, clamp=TRUE, rasterPreds=TRUE, parallel=TRUE, numCores=7)

eval.preds.bc2<-ENMevaluate(occ=sp.occ, env=preds.bc2, bg.coords=spBg.occ, method='jackknife', categoricals="soil_stype", RMvalues=reg.mult, fc=fc.comb, algorithm='maxent.jar', overlap=TRUE, bin.output=TRUE, clamp=TRUE, rasterPreds=TRUE, parallel=TRUE, numCores=7)

eval.preds.eco<-ENMevaluate(occ=sp.occ, env=preds.eco, bg.coords=spBg.occ, method='jackknife', categoricals="soil_stype", RMvalues=reg.mult, fc=fc.comb, algorithm='maxent.jar', overlap=TRUE, bin.output=TRUE, clamp=TRUE, rasterPreds=TRUE, parallel=TRUE, numCores=7)
}

##saving each of the runs to disk in case they need to be called later on
save(eval.preds, file=paste0(fold,"/", sp2[i], "_preds"))
save(eval.preds.bc2, file=paste0(fold,"/", sp2[i], "_bc2"))
save(eval.preds.eco, file=paste0(fold,"/", sp2[i], "_eco"))

##saving the input data to disk
#model1 - preds
write.csv(cbind(eval.preds@occ.pts, eval.preds@occ.grp), file=paste0(fold,"/", sp2[i], "_preds_occ.csv"), row.names=FALSE)
write.csv(cbind(eval.preds@bg.pts, eval.preds@bg.grp), file=paste0(fold,"/", sp2[i], "_preds_bg.csv"), row.names=FALSE)

#model2 - bc2
write.csv(cbind(eval.preds.bc2@occ.pts, eval.preds.bc2@occ.grp), file=paste0(fold,"/", sp2[i], "_bc2_occ.csv"), row.names=FALSE)
write.csv(cbind(eval.preds.bc2@bg.pts, eval.preds.bc2@bg.grp), file=paste0(fold,"/", sp2[i], "_bc2_bg.csv"), row.names=FALSE)

#model3 - eco
write.csv(cbind(eval.preds.eco@occ.pts, eval.preds.eco@occ.grp), file=paste0(fold,"/", sp2[i], "_eco_occ.csv"), row.names=FALSE)
write.csv(cbind(eval.preds.eco@bg.pts, eval.preds.eco@bg.grp), file=paste0(fold,"/", sp2[i], "_eco_bg.csv"), row.names=FALSE)

##saving the evaluation results to disk
#model1 - preds
write.csv(eval.preds@results, file=paste0(fold,"/", sp2[i], "_preds_results.csv"), row.names=FALSE)

#model2 - bc2
write.csv(eval.preds.bc2@results, file=paste0(fold,"/", sp2[i],"_bc2_results.csv"), row.names=FALSE)

#model3 - eco
write.csv(eval.preds.eco@results, file=paste0(fold,"/", sp2[i],"_eco_results.csv"), row.names=FALSE)
}

####