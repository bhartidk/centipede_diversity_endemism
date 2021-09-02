#export LD_LIBRARY_PATH=/usr/local/lib
#R CMD BATCH filename.r
#setwd("~/Downloads/maxent_analysis")

#set working directory
setwd("/home/bharti/D/PostDoc data/centipede_spatial")

#load necessary libraries
#library(sdmpredictors)
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
#library(rmaxent)
library(MASS)
#library(ecospat)
#library(combinat)
#library(ROCR)

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
wg.map<-fortify(wg)

##LOAD PREDICTOR LAYERS##
#load relevant input files from WorldClim - loading all files and loading select files which I had discussed with JJ which make sense for centipedes

#WORLDCLIM 1km layers - downloaded from worldclim.org on 26th March 2020
#wc<-list.files(path="dataframes/sdm_layers/wc2.1_30s_bio", pattern="*.tif")
wc<-list.files(path="/media/bharti/LITTLEONE", pattern="*.tif")

#The file is input, cropped to the required extent and put in the stack. 
wc_layers<-stack()
for(i in 1:length(wc)){
wc.raster<-raster(paste0("/media/bharti/LITTLEONE/", wc[i]))
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

##LOAD OCCURRENCE DATA##
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

##LOAD AND PROCESS ENMEvaluate RESULTS FILES FROM DISK##
#creating a list that stores all the results of each species as a separate object
enm.res<-list()
enm.res.sel<-list()
mod.sel<-c("rank", "aic")
mod.th<-c("mss", "mtp")
biodiverse<-data.frame()

#name the three kinds of models that were run for each species
mods<-c("preds", "bc2", "eco")

for(i in 1:length(sp)){
#provide the species specific folder name
dir.create(paste0("results/revisions/", unique(gps$Genus[gps$Species==sp[i]]), "_", sp[i], "/processed"))
fold<-paste0("results/revisions/", unique(gps$Genus[gps$Species==sp[i]]), "_", sp[i], "/processed")

#initializing elements of the list as a dataframe
enm.res[[i]]<-data.frame()
enm.res.sel[[i]]<-data.frame()
temp<-data.frame()

for(j in 1:length(mods)){
#load the ENMevaluation results file
load(paste0("results/revisions/", unique(gps$Genus[gps$Species==sp[i]]), "_", sp[i], "/", sp2[i], "_", mods[j]))
  
#load the ENMevaluation object as a new variable
if(j==1){
eval.model<-get(paste0("eval.", mods[j]), envir=.GlobalEnv)
}else{
eval.model<-get(paste0("eval.preds.", mods[j]), envir=.GlobalEnv)  
}

#save the results as a new object along with species and model type information
temp<-eval.model@results
temp<-cbind(species=rep(sp[i], times=nrow(temp)), model=rep(mods[j], times=nrow(temp)), temp)

#change the order of levels for temp$features and convert temp$rm into a factor
temp$features<-factor(temp$features, levels=c("L", "H", "Q", "LQ", "QH", "LQH"))
temp$rm<-as.factor(as.character(temp$rm))
temp$rm<-factor(temp$rm, levels=rev(seq(0.5, 5, by=0.5)))

#pool all the information into a single dataframe
enm.res[[i]]<-rbind(enm.res[[i]], temp)

for(k in 1:length(mod.sel)){
if(mod.sel[k]=="rank"){
#find the best model giving preference to minimum OR-MTP (minimum training presence omission rate) followed by minimum AUCdiff (difference between training and test AUC) followed by highest AUCtest
temp.sel<-temp %>% 
filter(avg.test.AUC > 0.6) %>%
filter(avg.test.orMTP==min(avg.test.orMTP)) %>%
filter(avg.diff.AUC==min(avg.diff.AUC)) %>%
filter(avg.test.AUC==max(avg.test.AUC)) %>%
as.data.frame()
#if there are no models with AUC-test greater than 0.5, assign NA values
if(nrow(temp.sel)==0){
next
}
}else if(mod.sel[k]=="aic"){
#based only on information theoretic approach, choose the best model based on the AICc values alone
temp.sel<-temp %>% 
filter(delta.AICc==0) %>%
as.data.frame()

#if there are no models that fit the criteria, move to next loop
if(nrow(temp.sel)==0){
next
}
}

#if there are multiple models which are filtered out, choosing the simplest one with the highest regularization value
if(nrow(temp.sel)>1){
temp.sel$fs<-as.numeric(temp.sel$features)
temp.sel$rms<-as.numeric(temp.sel$rm)
temp.sel<-temp.sel %>% 
filter(fs==min(fs)) %>%
filter(rms==min(rms))
temp.sel<-temp.sel[ ,colnames(temp)]
}
  
#find the id of the maxent models corresponding to each of the best models selected above
id.sel<-which(temp$settings %in% temp.sel$settings)

#saving the threshold values for each of the best models
#models ranked best in terms of OR-MTP, AUCdiff and AUCtest
model.sel<-eval.model@models[[id.sel]]
mss<-model.sel@results[grep("Maximum.training.sensitivity.plus.specificity.Cloglog.threshold", rownames(model.sel@results))]
mtp<-model.sel@results[grep("Minimum.training.presence.Cloglog.threshold", rownames(model.sel@results))]
tp10<-model.sel@results[grep("X10.percentile.training.presence.Cloglog.threshold", rownames(model.sel@results))]

#writing the model results to disk for future reference
write.csv(model.sel@results, paste0(fold, "/", sp2[i], "_", mods[j], "_", mod.sel[k], "_results.csv"))

#saving the model to disk - it says the object doesnot exist
save(model.sel, file=paste0(fold, "/", sp2[i], "_", mods[j], "_", mod.sel[k],  "_model"))

#saving the prediction raster to disk
#note that the predict function in dismo by default produces a complementary log-log output prediction map (look at the "A brief tutorial on Maxent, by Phillips, 2017"), but specifying it in the arguments anyway just so that I remember
predict.sel<-dismo::predict(model.sel, preds, args="outputformat=cloglog")
writeRaster(predict.sel, filename=paste0(fold, "/", sp2[i], "_", mods[j], "_", mod.sel[k], "_map"), format="GTiff", overwrite=TRUE)

#plot predicted map
jpeg(file=paste0(fold, "/", sp2[i], "_", mods[j], "_", mod.sel[k], "_map.jpeg"),  height=8, width=11, units='in', res=300)
plot(predict.sel, main=paste0(sp[i], "_", mods[j], "_", mod.sel[k]))
points(eval.model@occ.pts)
dev.off()

#plot variable percent contribution 
jpeg(file=paste0(fold, "/", sp2[i], "_", mods[j], "_", mod.sel[k], "_varCon.jpeg"), height=8, width = 11, units = 'in', res=300)
per.con<-var.importance(model.sel)[,c(1,2)]
per.con<-per.con[order(per.con$percent.contribution, decreasing=TRUE),]
per.con[,1]<-factor(per.con[,1], levels=per.con[order(per.con$percent.contribution, decreasing=TRUE), 1])
p<-ggplot(per.con, aes(x=variable, y=percent.contribution)) + geom_point() + ylim(0,100) + theme(axis.text.x=element_text(angle=90, hjust=0))
print(p)
dev.off()

#plot variable permutation importance
jpeg(file=paste0(fold, "/", sp2[i], "_", mods[j], "_", mod.sel[k], "_perImp.jpeg"), height=8, width = 11, units = 'in', res=300)
per.imp<-var.importance(model.sel)[,c(1,3)]
per.imp<-per.imp[order(per.imp$permutation.importance, decreasing=TRUE),]
per.imp[,1]<-factor(per.imp[,1], levels=per.imp[order(per.imp$permutation.importance, decreasing=TRUE), 1])
p<-ggplot(per.imp, aes(x=variable, y=permutation.importance)) + geom_point() + ylim(0,100) + theme(axis.text.x=element_text(angle=90, hjust=0))
print(p)
dev.off()

#saving the response curve
jpeg(file=paste0(fold, "/", sp2[i], "_", mods[j], "_", mod.sel[k], "_response.jpeg"),  height=8, width=11, units='in', res=300)
modRes<-response(model.sel, range='p')
dev.off()

#saving the thresholds as a vector
mod.th.values<-c(mss, mtp)

#pool all the information into a single dataframe
enm.res.sel[[i]]<-rbind(enm.res.sel[[i]], cbind(temp.sel, selection=mod.sel[k], mss, mtp, tp10))

##create and save the thresholded map##
for(l in 1:length(mod.th)){
#cropping the prediction raster to the extent of Western Ghats by first creating a template
#using mask directly on sdm.100 does not work as it removes cells with small overlaps with the polygon - which is unusual but I wasn't able to find a solution to this.
map<-raster::mask(predict.sel, wg)
map<-raster::aggregate(map, fact=100, fun=max, na.rm=TRUE)

#aggregate raster cells to a higher resolution
#0.0083 is 1 km and 0.83 is 100 km
sdm.100<-raster::aggregate(predict.sel, fact=100, fun=max, na.rm=TRUE)
sdm.100<-raster::mask(sdm.100, map)

#loading the sampling points corresponding to the ith species
loc.sp<-eval.model@occ.pts
coordinates(loc.sp)<-~LON+LAT
crs(loc.sp)<-wgs

#subsetting locations within western ghats alone
loc.sp<-loc.sp[!is.na(over(loc.sp, wg)[,1]),]

#find the cell positions of loc.sp and assign those as presence locations
#this is being done since some of the presence locations have been assigned low values of habitat suitability by the model, falling below the threshold and are removed from the thresholded predictions.
pres.cell<-cellFromXY(sdm.100, loc.sp) 
sdm.100[pres.cell]<-1

#convert the raster layer into a linear vector of raster values
sdm.vec<-as.data.frame(rasterToPoints(sdm.100))

#rename the last column
colnames(sdm.vec)[3]<-"prob"

#converting values less than the threshold to 0 and greater than or equal to threshold to 1
sdm.vec$prob<-ifelse(sdm.vec$prob<mod.th.values[l], 0, 1)

#use the above dataframe for plotting
plot.sdm<-sdm.vec
plot.sdm$prob<-as.character(plot.sdm$prob)
plot.sdm$prob<-factor(plot.sdm$prob, levels=c("0", "1"))

#plotting the thresholded raster and saving it to disk
#normal scale
my.fill<-scale_fill_manual(name=NULL, labels=c("Absence", "Presence"), values=c("grey", "#31a354"), drop=F)
my.scale<-scale_y_continuous(breaks=seq(8, 24, 2))

#plot for western ghats alone
plot<-ggplot(data=plot.sdm, aes(x=x, y=y)) + geom_raster(aes(fill=prob)) + geom_point(data=as.data.frame(loc.sp), aes(x=LON, y=LAT), shape=1) + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + xlab("Longitude") + ylab("Latitude") + ggtitle(paste0(sp[[i]], " ", mods[j], " ",  mod.sel[k])) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_fixed() + my.fill + my.scale

pdf(file=paste0(fold, "/", sp2[i], "_", mods[j], "_", mod.sel[k], "_", mod.th[l], "_map.pdf"), height=8, width=5, useDingbats=FALSE)
print(plot)
dev.off()

#add species, model and model selection criteria information to sdm vec
biodiverse<-as.data.frame(cbind(sdm.vec, species=rep(sp[i], times=nrow(sdm.vec)), model=rep(mods[j], times=nrow(sdm.vec)), model.selection=rep(mod.sel[k], times=nrow(sdm.vec))))

#save the biodiverse data to disk
write.csv(biodiverse, file=paste0(fold, "/", sp2[i], "_", mods[j], "_", mod.sel[k], "_", mod.th[l], "_biodiverse.csv"), row.names=FALSE)
}
}
}

#converting enm.res and enm.res.rank into a dataframe and saving it to disk
write.csv(as.data.frame(enm.res[[i]]), file=paste0(fold, "/", sp2[i], "_all_results.csv"), row.names=FALSE)

write.csv(enm.res.sel[[i]], file=paste0(fold, "/", sp2[i], "_best_models.csv"), row.names=FALSE)

#clear memory
gc()
}

####