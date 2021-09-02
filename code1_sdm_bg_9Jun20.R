#export LD_LIBRARY_PATH=/usr/local/lib
#R CMD BATCH filename.r
setwd("~/Downloads/maxent_analysis")

#set working directory
#setwd("/home/bharti/D/PostDoc data/spatial_data")

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
#library(ENMeval)
library(parallel)
library(devtools)
#library(rmaxent)
library(MASS)
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
names(elev)<-"elevation"
elev<-crop(elev, wg.ext)

pred<-stack(wc_layers, elev, soil.stype)

#load the spatial data file
#gps<-read.csv("GPS_data_20Mar20.csv")
gps<-read.csv("dataframes/Scolopendridae_5May20.csv")
levels(gps$Species.Name)

#look at the first few lines of the file
head(gps)

#looking at the structure of the dataframe
str(gps)

#drop voucher number
gps<-gps[,-c(1,2,4,5,6,8)]
head(gps)

#change column names for Latitude and longitude
colnames(gps)[c(3,4)]<-c("Latitude", "Longitude")

#changing the class of some of the columns
#gps$Voucher.No.<-as.character(gps$Voucher.No.)
gps$Location<-as.character(gps$Location)
gps$Latitude<-as.numeric(as.character(gps$Latitude))
gps$Longitude<-as.numeric(as.character(gps$Longitude))
gps$Species.Name<-as.character(gps$Species.Name)

#choosing only the rows which are from India - retain only rows which have "India" in the Location column
gps<-gps[grep("India", gps$Location), ]

#removing rows with NA values in Longitude and Latitude
gps<-gps[!(is.na(gps$Latitude) | is.na(gps$Longitude) | is.na(gps$Species.Name) | gps$Species.Name==""),]

#creating a new column with the genus names
gps$Genus<-rep(NA, times=nrow(gps))

gps$Genus[grep("Digitipes", gps$Species.Name)]<-"Digitipes"
gps$Genus[grep("Rhysida", gps$Species.Name)]<-"Rhysida"
gps$Genus[grep("Otostigmus", gps$Species.Name)]<-"Otostigmus"
gps$Genus[grep("Ethmostigmus", gps$Species.Name)]<-"Ethmostigmus"
gps$Genus[grep("Cormocephalus", gps$Species.Name)]<-"Cormocephalus"
gps$Genus[grep("Scolopendra", gps$Species.Name)]<-"Scolopendra"
gps$Genus[grep("Asanada", gps$Species.Name)]<-"Asanada"

#converting species and genus columns into factor   
gps$Genus<-factor(gps$Genus)
gps$Species.Name<-factor(gps$Species.Name)

#remove duplicate rows - having duplicate points for the same species doesn't help
nrow(gps)
gps<-unique(gps)
nrow(gps)

#remove species and genera which are no longer present in the dataframe
gps$Species.Name<-droplevels(gps$Species.Name)
gps$Genus<-droplevels(gps$Genus)

#writing the gps.pts object to disk
write.csv(gps, file=paste0("dataframes/gps_9Jun20.csv"), row.names=FALSE)

#removing the Scolopendra and Asanada group locations - because it was collected differently and retaining only the longitude and latitude columns
gps.pts<-unique(gps[gps$Genus!="Scolopendra" & gps$Genus!="Asanada" ,c(4,3)])
nrow(gps.pts) #the number of points is 186

#retaining only one point from each cell - these would essentially be environmental duplicates - using the code I has used for my PhD chapter
#Removing multiple records from the same cell - Adam Smith code
#this can also be mentioned as an argument within the maxent function but it is perhaps better to do it here
#inputing the function
source('R scripts/Eliminate Points in Same Cell of a Raster.r')

gps.pts<-elimCellDups(gps.pts, pred[[1]], longLatFields=c('Longitude', 'Latitude'), priority=NULL)
nrow(gps.pts) #there are now 133 points

#converting gps.pts into a SpatialPoints dataframe and setting the projection - this is to use the function extract
gps.sp.pts<-gps.pts
coordinates(gps.sp.pts)<-~Longitude+Latitude
projection(gps.sp.pts)<-wgs

#plot points to check that there aren't any outside India
plot(gps.sp.pts)
plot(pi.coast, add=TRUE, col="blue")

#there is a point way up north - can maybe remove it because it is outside the model extent
bb<-as(wg.ext, 'SpatialPolygons')
crs(bb)<-wgs

#plotting to see if we have got the right points
plot(gps.sp.pts)
plot(pi.coast, add=TRUE, col="blue")
plot(gps.sp.pts[is.na(over(gps.sp.pts, bb))], add=TRUE, col="red")

#remove this point
gps.sp.pts<-gps.sp.pts[-which(is.na(over(gps.sp.pts, bb)))]

#replot to check if this point has been removed
plot(gps.sp.pts)
plot(pi.coast, add=TRUE, col="blue")

#extract values from the continuous layers - since I am extracting these values for the bias layer, calling this bias.ext
bias.ext<-as.data.frame(extract(pred, gps.sp.pts))

#combining bias.ext with the coordinates
bias.ext<-cbind(as.data.frame(gps.sp.pts), bias.ext)
nrow(bias.ext) #132

#checking if there are any NA values here, and keeping only the ones with non-NA values
bias.ext<-bias.ext[complete.cases(bias.ext),]
nrow(bias.ext) #131

#extracting random background points - drawing more points than needed because earlier there were some points drawn where the categorical values had NAs.
#I am going to save these random background points and use it later as one means of building a Maxent model
set.seed(454)
randomSites<-randomPoints(pred, 11000)
colnames(randomSites)<-c("Longitude", "Latitude")
randomSites<-as.data.frame(randomSites)

#extracting environmental values for random data-points from pred
randomBg<-as.data.frame(extract(pred, randomSites))

#attaching the random pts coordinates with the extracted values
randomBg<-cbind(randomSites, randomBg)

#looking at the number of complete cases in rand.ext
nrow(randomBg[complete.cases(randomBg),])

#keeping 10,000 points from the complete cases in rand.ext
randomBg<-randomBg[complete.cases(randomBg),]
randomBg<-randomBg[1:10000,]

##CREATING INPUT DATA##
#Maxent model using random background points
trainData<-rbind(bias.ext, randomBg)
presentBg<-c(rep(1, times=nrow(bias.ext)), rep(0, times=nrow(randomBg)))

biasIn<-cbind(trainData, presentBg)
ncol(biasIn)
nrow(biasIn)
colnames(biasIn)

#checking if there are any NAs in the data
nrow(biasIn[!complete.cases(biasIn),])

#saving to disk
write.csv(biasIn, file=paste0("dataframes/biasIn.csv"), row.names=FALSE)

#saving just the random background points to disk
write.csv(randomBg, file=paste0("dataframes/randomBg.csv"), row.names=FALSE)

##RUNNING THE MAXENT MODEL##
#the meaning of iterations in the Maxent output - from Phillips et al. 2006 - "modeling species geographic distributions, we implemented an efficient algorithm together with a choice of feature types that are well suited to the task. Our implementation uses a sequential-update algorithm that iteratively picks a weight λj and adjusts it so as to minimize the resulting regularized log loss. The algorithm is deterministic, and is guaranteed to converge to the Maxent probability distribution. The algorithm stops when a user-specified number of iterations has been performed, or when the change in log loss in an iteration falls below a user-specified value (convergence), whichever happens first"

#the path variable saves the temp files of Maxent in the folder I am pointing to. The args argument asks to make predictions for the background points as well (need to look into this a little more)
biasModel<-maxent(x=biasIn[,3:(ncol(biasIn)-1)], p=biasIn[,ncol(biasIn)], factors=c('soil_stype', 'geol_lith'), args=c('writebackgroundpredictions=true', 'randomtestpoints=20'), path='/home/bharti/D/PostDoc data/spatial_data/dataframes/maxent_temp')

#saving model results to disk
write.csv(biasModel@results, paste0("dataframes/biasModel.csv"))

#saving the model to disk
save(biasModel, file="dataframes/biasModel")

#predicting the model to the extent of peninsular India
biasMap<-predict(biasModel, pred)

#saving the prediction raster to disk
writeRaster(biasMap, filename="dataframes/biasMap", format="GTiff", overwrite=TRUE)

#extract 10,000 random points from the bias raster
set.seed(30985)
biasSites<-randomPoints(biasMap, 10000, prob=TRUE)
biasEnv<-as.data.frame(extract(pred, biasSites))
biasBg<-cbind(biasSites, biasEnv)
colnames(biasBg)[1:2]<-c("Longitude", "Latitude")
nrow(biasBg[complete.cases(biasBg),]) #checking if there are NA extracts

#writing the bias background to disk
write.csv(biasBg, file="dataframes/biasBg.csv", row.names=FALSE)

####