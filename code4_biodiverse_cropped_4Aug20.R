###This code creates the input data for biodiverse (thresholding raw Maxent output and aggregating cells to a coarser resolution), uses the output of biodiverse to plot maps of diversity, calculates MPD and MNTD, and contains analysis of betadiversity patterns and generates associated plots. These different sections are clearly marked.

#set working directory - this should lead to the location where the /results, /dataframes and /R scripts folders are present
setwd("/home/bharti/D/PostDoc data/spatial_data")

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

####################################
####CREATE BIODIVERSE INPUT FILE####

#load the sampling locations that were saved to disk to be used in the supplementary material
loc<-read.csv("dataframes/supp_table_locations.csv", header=TRUE)

#create a vector of folder names where the sdm results have already been stored
loc$Species.Name2<-paste0(loc$Genus, "_", loc$Species)
sp<-unique(loc$Species.Name2)

#removing Rhysida_sp. 1 from this list because the distribution model is not significant
sp<-sp[-19]

#using the specific background we are interested in
bg.names<-"biasBg"

#creating names to store map results
sp2<-substr(unique(loc$Species), 1, 5)
sp2<-sp2[-19]

#saving unique genus and species names
sp.names<-colsplit(sp, "_", c("Genus", "Species"))

#create a list where results from different species will be stored
biodiv<-list()

#opening the folder corresponding to each species and loading the evaluate results from the species model
for(i in 1:length(sp)){

#assigning folder name
fold<-paste0("results/", sp[i], "/")

#find the filename of Maxent evaluate output for the ith species
eval.file<-list.files(path=fold, pattern=paste0("*", bg.names, "_eval"))

#load the above output - the object is saved with the name modEval
load(paste0(fold, eval.file))

#find the threshold which corresponds to maximum sum of sensitivity and specificity
mod.th<-threshold(modEval, stat='spec_sens')
print(mod.th)

#load the Maxent predictions of the ith species
sdm.file<-list.files(path=fold, pattern=paste0("*", bg.names, "_map.tif"))
sdm<-raster(paste0(fold, sdm.file))

#cropping the prediction raster to the extent of Western Ghats by first creating a template
#using mask directly on sdm.100 does not work as it removes cells with small overlaps with the polygon - which is unusual but I wasn't able to find a solution to this.
temp.1<-mask(sdm, wg)
temp<-raster::aggregate(temp.1, fact=100, fun=max, na.rm=TRUE)

#aggregate raster cells to a higher resolution
#0.0083 is 1 km and 0.83 is 100 km
sdm.100<-raster::aggregate(sdm, fact=100, fun=max, na.rm=TRUE)
sdm.100<-mask(sdm.100, temp)
sdm.1<-mask(sdm, temp.1)

#loading the sampling points corresponding to the ith species
loc.sp<-loc[loc$Species.Name2==sp[i],]
coordinates(loc.sp)<-~Longitude+Latitude
crs(loc.sp)<-wgs

#subsetting locations within western ghats alone
loc.sp<-loc.sp[!is.na(over(loc.sp, wg)[,1]),]

#find the cell positions of loc.sp and assign those as presence locations
#this is being done since some of the presence locations have been assigned low values of habitat suitability by the model, falling below the threshold and are removed from the thresholded predictions.
pres.cell<-cellFromXY(sdm.100, loc.sp) 
sdm.100[pres.cell]<-1

#convert the raster layer into a linear vector of raster values
sdm.vec.1km<-as.data.frame(rasterToPoints(sdm.1))
sdm.vec<-as.data.frame(rasterToPoints(sdm.100))

#rename the last column
colnames(sdm.vec.1km)[3]<-"prob"
colnames(sdm.vec)[3]<-"prob"

#converting values less than the threshold to 0 and greater than or equal to threshold to 1
sdm.vec.1km$prob<-ifelse(sdm.vec.1km$prob<mod.th, 0, 1)
sdm.vec$prob<-ifelse(sdm.vec$prob<mod.th, 0, 1)

#use the above dataframe for plotting
plot.sdm.1km<-sdm.vec.1km
plot.sdm.1km$prob<-as.factor(as.character(plot.sdm.1km$prob))

plot.sdm<-sdm.vec
plot.sdm$prob<-as.factor(as.character(plot.sdm$prob))

#removing the absence locations and saving the dataframe to disk
if(length(sdm.vec$prob[sdm.vec$prob==1])==0){
print(paste0("There are no presence locations for ", sp[i]))
}else{
#plotting the thresholded raster and saving it to disk
#normal scale
my.fill<-scale_fill_manual(name=NULL, labels=c("Absence", "Presence"), values=c("grey", "#31a354"))
my.scale<-scale_y_continuous(breaks=seq(8, 24, 2))

#plot for western ghats alone
plot.1km<-ggplot(data=plot.sdm.1km, aes(x=x, y=y)) + geom_raster(aes(fill=prob)) + xlab("Longitude") + ylab("Latitude") + ggtitle(sp[[i]]) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_fixed() + my.fill + my.scale

pdf(file=paste0(fold, sp2[i], "_", bg.names, "_map_th_1km.pdf"), height=8, width=5, useDingbats=FALSE)
print(plot.1km)
dev.off()

#plot for western ghats alone
plot<-ggplot(data=plot.sdm, aes(x=x, y=y)) + geom_raster(aes(fill=prob)) + geom_point(data=as.data.frame(loc.sp), aes(x=Longitude, y=Latitude), shape=1) + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + xlab("Longitude") + ylab("Latitude") + ggtitle(sp[[i]]) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_fixed() + my.fill + my.scale

pdf(file=paste0(fold, sp2[i], "_", bg.names, "_map_th.pdf"), height=8, width=5, useDingbats=FALSE)
print(plot)
dev.off()
}

#adding the species name to the presence-absence column
colnames(sdm.vec)[3]<-sp[i]

#saving this dataframe to the list
biodiv[[i]]<-sdm.vec  

#remove modEval
rm(modEval)
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

#Add locations for Rhysida sp. 1 where Maxent predictions were not used because the model was not significant
#creating a template where Rhysida sp. 1 points will be loaded
Rhysida_sp_1<-rasterFromXYZ(biodiv.mat[,1:3])
crs(Rhysida_sp_1)<-wgs
names(Rhysida_sp_1)<-"Rhysida_sp_1"

#turn the non-NA values of the template to 0
Rhysida_sp_1[!is.na(values(Rhysida_sp_1))]<-0

#select the locations of Rhysida sp. 1 which fall within WG
loc.sub<-loc[loc$Species.Name2=="Rhysida_sp. 1", 4:5]
coordinates(loc.sub)<-~Longitude+Latitude
crs(loc.sub)<-wgs

loc.sub<-loc.sub[!is.na(over(loc.sub, wg)[,1]),]

#find the cell positions of Rhysida sp. 1
pres.cell<-cellFromXY(Rhysida_sp_1, loc.sub) 
Rhysida_sp_1[pres.cell]<-1

#plot the new raster to see if the value has been assigned correctly
plot(Rhysida_sp_1)
points(loc.sub)

#converting this new raster into a dataframe
Rhysida_sp_1<-as.data.frame(rasterToPoints(Rhysida_sp_1))

#check if the lat longs are identical to biodiverse.mat
identical(Rhysida_sp_1[,1:2], biodiv.mat[,1:2]) #this turns out false but these locations are identical
setdiff(Rhysida_sp_1[,1:2], biodiv.mat[,1:2])
cbind(Rhysida_sp_1[,1], biodiv.mat[,1]) #they are identical
cbind(Rhysida_sp_1[,2], biodiv.mat[,2]) #they are identical

#adding the Rhysida_sp_1 data to biodiv.mat
biodiv.mat$Rhysida_sp_1<-Rhysida_sp_1[,3]
colnames(biodiv.mat)[21]<-"Rhysida_sp. 1"
head(biodiv.mat)

#change presences for D. coonoorensis, E. agasthyamalaensis, E. coonooranus, R. sada
#Decided not to meddle with E. coonooranus and R. trispinosa

#For D. coonoorensis, removing predictions above 14 degrees
biodiv.mat$Digitipes_coonoorensis[biodiv.mat$y>14]<-0

#For E. agasthyamalaensis, removing predictions above 10.5 degrees
biodiv.mat$Ethmostigmus_agasthyamalaiensis[biodiv.mat$y>10.5]<-0

#For E. coonooranus, removing predictions above 16.5
biodiv.mat$Ethmostigmus_coonooranus[biodiv.mat$y>16.5]<-0

#For R. sada removing predictions below 17
biodiv.mat$Rhysida_sada[biodiv.mat$y<17]<-0

#plotting the new thresholded graphs for these species
plot.sub<-biodiv.mat[,c(1,2,4,8,10,19,21)]
sp.names.sub<-sp[c(2,6,8,17)]
sp.names.sub<-c(sp.names.sub, "Rhysida_sp. 1")
sp2.sub<-c(sp2[c(2,6,8,17)], "sp. 1")

#normal scale
my.fill<-scale_fill_manual(name=NULL, labels=c("Absence", "Presence"), values=c("grey", "#31a354"))
my.scale<-scale_y_continuous(breaks=seq(8, 24, 2))

for(i in 1:length(sp.names.sub)){
#loading the sampling points corresponding to the ith species
loc.sp<-loc[loc$Species.Name2==sp.names.sub[i],]
coordinates(loc.sp)<-~Longitude+Latitude
crs(loc.sp)<-wgs

#subsetting locations within western ghats alone
loc.sp<-loc.sp[!is.na(over(loc.sp, wg)[,1]),]

plot.sdm<-plot.sub[,c("x", "y", sp.names.sub[i])]
colnames(plot.sdm)[3]<-"prob"

plot.sdm$prob<-as.factor(as.character(plot.sdm$prob))

#plot for western ghats alone
plot<-ggplot(data=plot.sdm, aes(x=x, y=y)) + geom_raster(aes(fill=prob)) + geom_point(data=as.data.frame(loc.sp), aes(x=Longitude, y=Latitude), shape=1) + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + xlab("Longitude") + ylab("Latitude") + ggtitle(sp2.sub[[i]]) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_fixed() + my.fill + my.scale

pdf(file=paste0("results/", sp.names.sub[i], "/", sp2.sub[i], "_", bg.names, "_map_th_edit.pdf"), height=8, width=5, useDingbats=FALSE)
print(plot)
dev.off()

rm(plot.sdm)
}

#give id names to the rows
biodiv.mat$id<-seq(1:nrow(biodiv.mat))
colnames(biodiv.mat)
biodiv.mat<-biodiv.mat[,c(22,1:21)]

#changing the colname of Rhysida sp. 1
colnames(biodiv.mat)[22]<-"Rhysida_sp_1"

#writing this to disk as a .csv file
write.csv(biodiv.mat, file="dataframes/biodiverse/biodiverse_edit_6Jul20.csv", row.names=FALSE)

##
#creating the input tree file for biodiverse
#loading the tree file
cen.tree<-phylo4(read.tree("dataframes/biodiverse/threegenera.tre"))

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

#######################################
####ANALYZE BIODIVERSE OUTPUT FILES####
#load biodiverse input file
bd.ip<-read.csv("dataframes/biodiverse/biodiverse_edit_6Jul20.csv", header=TRUE)

#load biodiverse results
bd<-read.csv("results/biodiverse/threegenera_results_6Jul20.csv", header=TRUE)
bd.r<-read.csv("results/biodiverse/threegenera_struc_randomization_results_6Jul20.csv", header=TRUE)

#keeping the columns of interest from bd
colnames(bd)
identical(bd$ENDW_RICHNESS, bd$RICHNESS_ALL) #identical - keeping only richness
identical(bd$RICHNESS_ALL, bd$RICHNESS_SET1) #identical - keeping only richness
identical(bd$ENDW_CWE, bd$ENDW_WE) #not the same, ENDW_CWE is weighted by richness - but I don't need this, keeping only ENDW_WE
identical(bd$ENDW_WE, bd$ENDW_SINGLE) #these are the same, retaining only ENDW_WE
#_P is scaled by tree length, _per_taxon is proportion contribution by each taxon 0 I'm choosing not to consider these for now - but I think the scaled indices are what have been reported in the papers we have discussed so far. _RPD_DIFF2 tells us how much more or less PD is there than expected, _RPD_NULL - gives the value of the denominator used for RPD calculations - choosing to ignore these as well.

#based on the above reference, only keeping columns that are of interest to us
bd<-bd[,c(2,19,6,7,8,11,12,13,15,16,18)]

#renaming columns
colnames(bd)
colnames(bd)<-c("id", "sr", "we", "pd", "pd.p", "pe", "pe.p", "rpd", "rpd.denom", "rpe", "rpe.denom")

#keeping the column names of interest from bd.r
colnames(bd.r) #what we are interested in is the percentiles - ranking of original value against randomizations
bd.r<-bd.r[,c(2,25,26,29,30,31,33,34,36)]
colnames(bd.r)<-c("id", "p.pd", "p.pd.p", "p.pe", "p.pe.p", "p.rpd", "p.rpd.denom", "p.rpe", "p.rpe.denom") 

#left_join bd.ip and bd
bd.all<-left_join(bd.ip, bd, by="id")
bd.all<-left_join(bd.all, bd.r, by="id")

#eyeballing the distribution of various values
par(mfrow=c(2,2))
hist((bd.all$p.pe)*100, xlab="Percentile", main="Histogram of PE percentiles")
hist((bd.all$p.rpe)*100, xlab="Percentile", main="Histogram of RPE percentiles")
hist((bd.all$p.pd)*100, xlab="Percentile", main="Histogram of PD percentiles")
hist((bd.all$p.rpd)*100, xlab="Percentile", main="Histogram of RPD percentiles")

#plot sr, pd, we, pe
#normal scale
my.fill1<-scale_fill_gradientn(name=NULL, colours=rev(c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')))

  
plot.raw<-list()

#plot diversity measures for western ghats alone
plot.raw[[1]]<-ggplot(data=bd.all, aes(x=x, y=y)) + geom_raster(aes(fill=sr)) + my.fill1 + ggtitle("SR") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.scale + coord_fixed()

plot.raw[[2]]<-ggplot(data=bd.all, aes(x=x, y=y)) + geom_raster(aes(fill=we)) + my.fill1 + ggtitle("WE") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.scale + coord_fixed()

plot.raw[[3]]<-ggplot(data=bd.all, aes(x=x, y=y)) + geom_raster(aes(fill=pd.p)) + my.fill1 + ggtitle("PD") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.scale + coord_fixed()

plot.raw[[4]]<-ggplot(data=bd.all, aes(x=x, y=y)) + geom_raster(aes(fill=pe.p)) + my.fill1 + ggtitle("PE") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.scale + coord_fixed()

plot.raw[[5]]<-ggplot(data=bd.all, aes(x=x, y=y)) + geom_raster(aes(fill=rpd)) + my.fill1 + ggtitle("RPD") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.scale + coord_fixed()

plot.raw[[6]]<-ggplot(data=bd.all, aes(x=x, y=y)) + geom_raster(aes(fill=rpe)) + my.fill1 + ggtitle("RPE") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.scale + coord_fixed()

#plot percentiles of different diversity measures as such
my.fill<-scale_fill_gradientn(name=NULL, colours=rev(c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')), breaks=seq(0,1, by=0.2), limits=c(0,1))

my.scale<-scale_y_continuous(breaks=seq(8, 24, 2))

plot.q<-list()

plot.q[[1]]<-ggplot(data=bd.all, aes(x=x, y=y)) + geom_raster(aes(fill=p.pd.p)) + ggtitle("PD-percentile") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill + my.scale + coord_fixed()

plot.q[[2]]<-ggplot(data=bd.all, aes(x=x, y=y)) + geom_raster(aes(fill=p.rpd)) + ggtitle("RPD-percentile") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill + my.scale + coord_fixed()

plot.q[[3]]<-ggplot(data=bd.all, aes(x=x, y=y)) + geom_raster(aes(fill=p.rpd.denom)) + ggtitle("Null PD-percentile") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill + my.scale + coord_fixed()

#plot diversity measures for western ghats alone
plot.q[[4]]<-ggplot(data=bd.all, aes(x=x, y=y)) + geom_raster(aes(fill=p.pe.p)) + ggtitle("PE-percentile") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill + my.scale + coord_fixed()

plot.q[[5]]<-ggplot(data=bd.all, aes(x=x, y=y)) + geom_raster(aes(fill=p.rpe)) + ggtitle("RPE-percentile") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill + my.scale + coord_fixed()

plot.q[[6]]<-ggplot(data=bd.all, aes(x=x, y=y)) + geom_raster(aes(fill=p.rpe.denom)) + ggtitle("Null PE-percentile") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill + my.scale + coord_fixed()

#sub-selecting only those sites where there is significantly high PD - falls within the highest 2.5 percentile
#plot pe, significant pe, paleoendemism and neoendemism
bd.all$pd.cat<-cut(bd.all$p.pd.p, breaks=c(0, 0.95, 0.99, 1), labels=c("Not Significant", ">0.95", ">0.99"))
bd.all$rpd.cat<-cut(bd.all$p.rpd, breaks=c(0, 0.01, 0.025, 0.975, 0.99, 1), labels=c("<0.01", "<0.025", "Not Significant", ">0.975", ">0.99"))

bd.all$pe.cat<-cut(bd.all$p.pe.p, breaks=c(0, 0.95, 0.99, 1), labels=c("Not Significant", ">0.95", ">0.99"))
bd.all$rpe.cat<-cut(bd.all$p.rpe, breaks=c(0, 0.01, 0.025, 0.975, 0.99, 1), labels=c("<0.01", "<0.025", "Not Significant", ">0.975", ">0.99"))

bd.all$rpd.denom.cat<-cut(bd.all$p.rpd.denom, breaks=c(0, 0.95, 0.99, 1), labels=c("Not Significant", ">0.95", ">0.99"))
bd.all$rpe.denom.cat<-cut(bd.all$p.rpe.denom, breaks=c(0, 0.95, 0.99, 1), labels=c("Not Significant", ">0.95", ">0.99"))

#plot the different significance levels
#normal scale
my.fill<-scale_fill_manual(name=NULL, labels=c("<0.01", "<0.025", "Not\nSignificant", ">0.975", ">0.99"), values=c('#a50f15','#de2d26','#fee5d9','#2b8cbe','#045a8d'), drop=FALSE)

my.fill.pde<-scale_fill_manual(name=NULL, labels=c("Not\nSignificant", ">0.95", ">0.99"), values=c('#fee5d9','#2b8cbe','#045a8d'), drop=FALSE)

my.scale<-scale_y_continuous(breaks=seq(8, 24, 2))

plot.cat<-list()

#plot diversity measures for western ghats alone
#plot pd
plot.cat[[1]]<-ggplot(data=bd.all[complete.cases(bd.all),], aes(x=x, y=y)) + geom_raster(aes(fill=pd.cat)) + ggtitle("Phylogenetic diversity") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill.pde + my.scale + coord_fixed()

#plot rpd
plot.cat[[2]]<-ggplot(data=bd.all[complete.cases(bd.all),], aes(x=x, y=y)) + geom_raster(aes(fill=rpd.cat)) + ggtitle("Relative Phylogenetic diversity") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill + my.scale + coord_fixed()

#plot rpd.denom
plot.cat[[3]]<-ggplot(data=bd.all[complete.cases(bd.all),], aes(x=x, y=y)) + geom_raster(aes(fill=rpd.denom.cat)) + ggtitle("Null Phylogenetic diversity") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill.pde + my.scale + coord_fixed()

#plot pe
plot.cat[[4]]<-ggplot(data=bd.all[complete.cases(bd.all),], aes(x=x, y=y)) + geom_raster(aes(fill=pe.cat)) + ggtitle("Phylogenetic endemism") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill.pde + my.scale + coord_fixed()

#plot rpe
plot.cat[[5]]<-ggplot(data=bd.all[complete.cases(bd.all),], aes(x=x, y=y)) + geom_raster(aes(fill=rpe.cat)) + ggtitle("Relative Phylogenetic endemism") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill + my.scale + coord_fixed()

#plot rpe.denom
plot.cat[[6]]<-ggplot(data=bd.all[complete.cases(bd.all),], aes(x=x, y=y)) + geom_raster(aes(fill=rpe.denom.cat)) + ggtitle("Null Phylogenetic endemism") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill.pde + my.scale + coord_fixed()

#arranging the above graphs
pdf(file=paste0("results/results_figures/diversity_canape_significance_biodiverse.pdf"), height=8.5, width=11, useDingbats=FALSE)
ggarrange(plot.cat[[1]], plot.cat[[3]], plot.cat[[2]], plot.cat[[4]], plot.cat[[6]], plot.cat[[5]], ncol=3, nrow=2)
dev.off()

#creating a new column called endemism classification
bd.all$canape<-rep("not significant", times=nrow(bd.all))
bd.all$canape[(bd.all$p.pe>0.95 | bd.all$p.rpe.denom>0.95) & bd.all$p.rpe>0.975]<-"palaeoendemism"
bd.all$canape[(bd.all$p.pe>0.95 | bd.all$p.rpe.denom>0.95) & bd.all$p.rpe<0.025]<-"neoendemism"
bd.all$canape[(bd.all$p.pe>0.95 | bd.all$p.rpe.denom>0.95) & (bd.all$p.rpe>0.025 | bd.all$p.rpe<0.975)]<-"mixed endemism"

bd.all$canape<-as.factor(bd.all$canape)
levels(bd.all$canape)<-c(levels(bd.all$canape), "neoendemism", "palaeoendemism")

bd.all$canape<-factor(bd.all$canape, levels=c("neoendemism", "palaeoendemism", "mixed endemism", "not significant"))

#plot CANAPE results
my.fill<-scale_fill_manual(name=NULL, values=c('#a50f15','#045a8d','#9e9ac8', '#fee5d9'), drop=FALSE)

my.scale<-scale_y_continuous(breaks=seq(8, 24, 2))

pdf(file=paste0("results/results_figures/canape_classification_biodiverse.pdf"), height=11, width=8.5, useDingbats=FALSE)
ggplot(data=bd.all[complete.cases(bd.all),], aes(x=x, y=y)) + geom_raster(aes(fill=canape)) + ggtitle("CANAPE") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill + my.scale + coord_fixed()
dev.off()

################
##MPD and MNTD##
#load biodiverse input file
bd.ip<-read.csv("dataframes/biodiverse/biodiverse_edit_6Jul20.csv", header=TRUE)
bd.mat<-as.matrix(bd.ip[,4:22])
rownames(bd.mat)<-bd.ip$id

#load the phylogenetic tree file
#loading the tree file
cen.tree<-read.tree("dataframes/biodiverse/threegenera_biodiverse.tre")

#prune the tree to have the same labels as in the sitesxspecies matrix
cen.tree<-prune.sample(bd.mat, cen.tree)

#calculating MPD
#it requires a distance matrix - which is obtained by using the cophenetic function from the stats package
#I think cophenetic just calculates the branch lengths between species - takes a tree file (a dendrogram) as the input - the documentation says it takes a heirarchical clustering R object as input
cen.dist<-cophenetic(cen.tree)
cen.mpd<-mpd(bd.mat, cen.dist, abundance.weighted=FALSE)

#calculating standard effect size around cen.mpd
set.seed(53823)
cen.ses.mpd<-ses.mpd(bd.mat, cen.dist, null.model="taxa.labels", abundance.weighted=FALSE, runs=999, iterations=10000)

#adding the id column
cen.ses.mpd$id<-as.numeric(row.names(cen.ses.mpd))

#calculating MNTD
cen.mntd<-mntd(bd.mat, cen.dist, abundance.weighted=FALSE)
set.seed(53823)
cen.ses.mntd<-ses.mntd(bd.mat, cen.dist, null.model="taxa.labels", abundance.weighted=FALSE, runs=999, iterations=10000)
cen.ses.mntd$id<-as.numeric(row.names(cen.ses.mntd))

#calculating Faith's PD - checked that the results match with results from biodiverse when include.root=TRUE
set.seed(53823)
cen.ses.pd<-ses.pd(bd.mat, cen.tree, null.model="taxa.labels", runs=999, iterations=10000, include.root=TRUE)
cen.ses.pd$id<-as.numeric(row.names(cen.ses.pd))

#left_join the mpd and mntd dataframes
#checked and found that ntaxa for both ses.mpd and ses.mntd matched
cen.ses<-left_join(cen.ses.mpd, cen.ses.mntd, by=c('id', 'ntaxa', 'runs'))
cen.ses<-left_join(cen.ses, cen.ses.pd, by=c('id', 'ntaxa', 'runs'))
cen.ses<-left_join(cen.ses, bd.ip[,c(1:3)], by='id')

#find the range of finite values
pd.range<-range(cen.ses$pd.obs.z[!is.na(cen.ses$pd.obs.z)])
mpd.range<-range(cen.ses$mpd.obs.z[!is.na(cen.ses$mpd.obs.z)])
mntd.range<-range(cen.ses$mntd.obs.z[!is.na(cen.ses$mntd.obs.z)])

#the column mpd.obs.z is the standard effect size - plotting that on the map
my.fill<-scale_fill_gradientn(name=NULL, colours=rev(c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')))

#my.fill2<-scale_fill_gradientn(name=NULL, colours=rev(c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')), breaks=seq(-2, 2, by=0.5), limits=mpd.range)

#my.fill3<-scale_fill_gradientn(name=NULL, colours=rev(c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')), breaks=seq(-2, 2.5, by=0.5), limits=mntd.range)

my.scale<-scale_y_continuous(breaks=seq(8, 24, 2))

plot.ses<-list()

plot.ses[[1]]<-ggplot(data=cen.ses[complete.cases(cen.ses),], aes(x=x, y=y)) + geom_raster(aes(fill=pd.obs.z)) + ggtitle("SES-PD") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill + my.scale + coord_fixed()

plot.ses[[2]]<-ggplot(data=cen.ses[complete.cases(cen.ses),], aes(x=x, y=y)) + geom_raster(aes(fill=mpd.obs.z)) + ggtitle("SES-MPD") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill + my.scale + coord_fixed()

plot.ses[[3]]<-ggplot(data=cen.ses[complete.cases(cen.ses),], aes(x=x, y=y)) + geom_raster(aes(fill=mntd.obs.z)) + ggtitle("SES-MNTD") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill + my.scale + coord_fixed()

#arranging the above graphs
pdf(file=paste0("results/results_figures/SR_biodiverse.pdf"),  height=8.5, width=11, useDingbats=FALSE)
print(plot.raw[[1]])
dev.off()

pdf(file=paste0("results/results_figures/PD_MPD_MNTD.pdf"),  height=8.5, width=11, useDingbats=FALSE)
ggarrange(plot.raw[[3]], plot.raw[[5]], plot.q[[2]], plot.ses[[1]], plot.ses[[2]], plot.ses[[3]], ncol=3, nrow=2)
dev.off()

pdf(file=paste0("results/results_figures/WE_biodiverse.pdf"),  height=8.5, width=11, useDingbats=FALSE)
print(plot.raw[[2]])
dev.off()

pdf(file=paste0("results/results_figures/PE_biodiverse.pdf"),  height=8.5, width=11, useDingbats=FALSE)
ggarrange(plot.raw[[4]], plot.q[[4]], plot.raw[[6]], plot.q[[5]], ncol=2, nrow=2)
dev.off()

##################
##Beta diversity##
#I was not able to load the cluster analysis results from Biodiverse properly, so redoing the analysis here
#load biodiverse input file
bd.ip<-read.csv("dataframes/biodiverse/biodiverse_edit_6Jul20.csv", header=TRUE)
bd.mat<-as.matrix(bd.ip[,4:22])
rownames(bd.mat)<-bd.ip$id

#load the phylogenetic tree file
#loading the tree file
cen.tree<-read.tree("dataframes/biodiverse/threegenera_biodiverse.tre")

##Using Leprieur's R code##
#call the functions from outside
source("R scripts/PBD_Leprieur_2012.R")

#calculate components of beta phylogenetic diversity
set.seed(73635)
bd.comp<-beta.pd.decompo(com=bd.mat, tree=cen.tree, type="PhyloSor", random=999)
head(bd.comp)
#checked that PhyloSor calculated here and from the picante package give the same results

#split the first column into site 1 and site 2
new.cols<-colsplit(row.names(bd.comp$betadiv),"-",c("site1","site2"))
bd.res<-cbind(new.cols, bd.comp$betadiv)

#classifying the sites as SWG, Nilgiris, CWG, NWG based on the id in the site1 and site2 columns
sw<-c(38:45)
nilg<-c(30:37)
cw<-c(18:29)
nw<-c(1:17)

bd.res$site1.cat[bd.res$site1 %in% sw]<-"SWG"
bd.res$site1.cat[bd.res$site1 %in% nilg]<-"Nilgiris"
bd.res$site1.cat[bd.res$site1 %in% cw]<-"CWG"
bd.res$site1.cat[bd.res$site1 %in% nw]<-"NWG"

bd.res$site2.cat[bd.res$site2 %in% sw]<-"SWG"
bd.res$site2.cat[bd.res$site2 %in% nilg]<-"Nilgiris"
bd.res$site2.cat[bd.res$site2 %in% cw]<-"CWG"
bd.res$site2.cat[bd.res$site2 %in% nw]<-"NWG"

#create a new column of the comparison
bd.res$site.pair<-paste0(bd.res$site1.cat, "-", bd.res$site2.cat)

#reordering columns  
bd.res<-bd.res[,c(1,2,9:11,3:8)]

#converting site.pair into a factor
bd.res$site.pair<-factor(bd.res$site.pair, levels=c("SWG-SWG", "Nilgiris-Nilgiris", "CWG-CWG", "NWG-NWG", "Nilgiris-SWG", "CWG-SWG", "NWG-SWG", "CWG-Nilgiris", "NWG-Nilgiris", "NWG-CWG"))

#creating a new column which tells us if phylogenetic turn-over is significantly different from expected
bd.res$phylosorturn.sig<-rep("ns", times=nrow(bd.res))
bd.res$phylosorturn.sig[bd.res$SES_PhyloSor_turn>1.96]<-">97.5%"
bd.res$phylosorturn.sig[bd.res$SES_PhyloSor_turn<(-1.96)]<-"<2.5%"

#creating a new column which tells us if phylogenetic turn-over is significantly different from expected
bd.res$phylosorpd.sig<-rep("ns", times=nrow(bd.res))
bd.res$phylosorpd.sig[bd.res$SES_PhyloSor_PD>1.96]<-">97.5%"
bd.res$phylosorpd.sig[bd.res$SES_PhyloSor_PD<(-1.96)]<-"<2.5%"

#looking at the distribution of number of significant differences where turnover is greater or lower than random expectation
PBDturn.sig.high<-bd.res[bd.res$phylosorturn.sig==">97.5%", c("site.pair", "phylosorturn.sig")]
PBDturn.sig.low<-bd.res[bd.res$phylosorturn.sig=="<2.5%", c("site.pair", "phylosorturn.sig")]
PBDturn.sig.all<-rbind(PBDturn.sig.high, PBDturn.sig.low)
PBDturn.sig.all$phylosorturn.sig<-factor(PBDturn.sig.all$phylosorturn.sig, levels=c(">97.5%", "<2.5%"))

my.fill<-scale_fill_manual(name="Significance", labels=c(">0.975", "<0.025"), values=c('#de2d26', '#2b8cbe'), drop=FALSE)
my.scale<-scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, 2))

plot.sig<-list()

plot.sig[[1]]<-ggplot(data=PBDturn.sig.all, aes(x=site.pair, fill=phylosorturn.sig)) + geom_bar(position=position_dodge(preserve="single")) + xlab("Comparison") + ylab("No. of sig. comparisons") + theme_bw() + theme(text=element_text(size=16), axis.text.x=element_text(angle=45, hjust=1)) + my.fill + my.scale + scale_x_discrete(drop=FALSE)

#looking at the distribution of number of significant differences where beta diversity PD is greater or lower than random expectation
PBDpd.sig.high<-bd.res[bd.res$phylosorpd.sig==">97.5%", c("site.pair", "phylosorpd.sig")]
PBDpd.sig.low<-bd.res[bd.res$phylosorpd.sig=="<2.5%", c("site.pair", "phylosorpd.sig")]
PBDpd.sig.all<-rbind(PBDpd.sig.high, PBDpd.sig.low)
PBDpd.sig.all$phylosorpd.sig<-factor(PBDpd.sig.all$phylosorpd.sig, levels=c(">97.5%", "<2.5%"))

plot.sig[[2]]<-ggplot(data=PBDpd.sig.all, aes(x=site.pair, fill=phylosorpd.sig)) + geom_bar(position=position_dodge(preserve="single")) + xlab("Comparison") + ylab("No. of sig. comparisons") + theme_bw() + theme(text=element_text(size=16), axis.text.x=element_text(angle=45, hjust=1)) + my.fill + my.scale + scale_x_discrete(drop=FALSE)

pdf(file=paste0("results/results_figures/PBD_significance_boxplot.pdf"), height=8, width=8, useDingbats=FALSE)
ggarrange(plot.sig[[2]], plot.sig[[1]], ncol=1, nrow=2)
dev.off()

#getting some basic statistics out from the results
phylosor<-sum(bd.res$PhyloSor)
phylosor.turn<-sum(bd.res$PhyloSor_turn)
phylosor.pd<-sum(bd.res$PhyloSor_PD)

phylosor.turn/phylosor
phylosor.pd/phylosor

#plotting boxplot of phylosor.ses for different comparisons
my.scale<-scale_y_continuous(limits=c(-3, 3), breaks=seq(-2, 2, 1))

plot.bd<-list()

plot.bd[[1]]<-ggplot(data=bd.res[is.finite(bd.res$SES_PhyloSor),], aes(x=site.pair, y=SES_PhyloSor)) + geom_boxplot() + xlab("Comparison") + ylab("PhyloSor-SES") + geom_hline(aes(yintercept=0)) + theme_bw() +  theme(text=element_text(size=16), axis.text.x=element_text(angle=45, hjust=1)) + my.scale

plot.bd[[2]]<-ggplot(data=bd.res[is.finite(bd.res$SES_PhyloSor_turn),], aes(x=site.pair, y=SES_PhyloSor_turn)) + geom_boxplot() + xlab("Comparison") + ylab("PhyloSor(Turn)-SES") + geom_hline(aes(yintercept=0)) + theme_bw() +  theme(text=element_text(size=16), axis.text.x=element_text(angle=45, hjust=1)) + my.scale

plot.bd[[3]]<-ggplot(data=bd.res[is.finite(bd.res$SES_PhyloSor_PD),], aes(x=site.pair, y=SES_PhyloSor_PD)) + geom_boxplot() + xlab("Comparison") + ylab("PhyloSor(PD)-SES") + geom_hline(aes(yintercept=0)) + theme_bw() +  theme(text=element_text(size=16), axis.text.x=element_text(angle=45, hjust=1)) + my.scale

#plot the above figures
pdf(file=paste0("results/results_figures/PBD_components_boxplot.pdf"), height=8, width=8, useDingbats=FALSE)
ggarrange(plot.bd[[3]], plot.bd[[2]], ncol=1, nrow=2)
dev.off()

##SORENSEN INDEX##
#calculate beta-diversity-Sorensen dissimilarity index using vegan
sor.dist<-vegdist(bd.mat, method="bray", binary=TRUE)
sor<-as.matrix(sor.dist)
sor[upper.tri(sor, diag=TRUE)]<-NA
sor<-melt(sor)
colnames(sor)<-c("site1", "site2", "sor")

#calculate Phylogenetic Sorensen's index
phylo.sor<-phylosor(bd.mat, cen.tree)

#phylosor calculates the fraction of PD (branch lengths) shared between two sites
phylo.sor.dist<-1-phylo.sor

#carrying out clustering analysis with sor.dist and phylo.sor.dist
sor.tree<-hclust(sor.dist, method="average")
sor.dend<-as.dendrogram(sor.tree)
sor.phylo<-as.phylo(sor.tree)

phylo.sor.tree<-hclust(phylo.sor.dist, method="average")
phylo.sor.dend<-as.dendrogram(phylo.sor.tree)
phylo.sor.phylo<-as.phylo(phylo.sor.tree)

#plot the dedrograms from beta-diversity and phylogenetic beta-diversity side by side
#matching the same names on the lat longs of dataframe
sor.clust<-cutree(sor.tree, k=4) 
phylosor.clust<-cutree(phylo.sor.tree, k=4)
clust.df<-data.frame(id=bd.ip$id, x=bd.ip$x, y=bd.ip$y, sor.clust=sor.clust, phylosor.clust=phylosor.clust)

pdf(file=paste0("results/results_figures/sorensen_upgma.pdf"),  height=8.5, width=11, useDingbats=FALSE)
par(mfrow=c(1,2))
sor.dend %>% color_branches(col=rev(c("#d7191c", "#fdae61", "#a6d96a", "#1a9641"))[clust.df[labels(sor.dend),"sor.clust"]]) %>% set("labels_cex", 0.6) %>% plot(main="Sorensen-UPGMA")

phylo.sor.dend %>% color_branches(col=rev(c("#d7191c", "#fdae61", "#a6d96a", "#1a9641"))[clust.df[labels(phylo.sor.dend),"phylosor.clust"]]) %>% set("labels_cex", 0.6) %>% plot(main="PhyloSorensen-UPGMA")
dev.off()

#plot cluster analysis results
clust.df$sor.clust<-factor(clust.df$sor.clust)
clust.df$phylosor.clust<-factor(clust.df$phylosor.clust)

my.fill<-scale_fill_manual(name=NULL, values=rev(c("#d7191c", "#fdae61", "#a6d96a", "#1a9641")), drop=FALSE)

my.scale<-scale_y_continuous(breaks=seq(8, 24, 2))

plot.beta=list()

plot.beta[[1]]<-ggplot(data=clust.df, aes(x=x, y=y)) + geom_raster(aes(fill=sor.clust)) + ggtitle("Sorensen dissimilarity") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill + my.scale + coord_fixed()

plot.beta[[2]]<-ggplot(data=clust.df, aes(x=x, y=y)) + geom_raster(aes(fill=phylosor.clust)) + ggtitle("PhyloSorensen dissimilarity") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill + my.scale + coord_fixed()

pdf(file=paste0("results/results_figures/sorensen_map.pdf"), height=8.5, width=11, useDingbats=FALSE)
ggarrange(plot.beta[[1]], plot.beta[[2]], ncol=2, nrow=1)
dev.off()

#annotate raster cell numbers 
pdf(file=paste0("results/results_figures/wg_raster_cell_no.pdf"), height=8.5, width=11, useDingbats=FALSE)
plot.beta[[1]] + annotate("text", x=clust.df$x, y=clust.df$y, label=clust.df$id)
dev.off()

##SIMPSON INDEX##
beta<-beta.pair(bd.mat, index.family="sorensen")
sim<-beta$beta.sim
sim.tree<-as.dendrogram(hclust(sim, method="average"))
sim.dend<-as.dendrogram(sim.tree)

phylobeta<-phylo.beta.pair(bd.mat, cen.tree, index.family="sorensen")
phylosim<-phylobeta$phylo.beta.sim
phylosim.tree<-as.dendrogram(hclust(phylosim, method="average"))
phylosim.dend<-as.dendrogram(phylosim.tree)

#matching the same names on the lat longs of dataframe
sim.clust<-cutree(sim.tree, k=4) 
phylosim.clust<-cutree(phylosim.tree, k=4)
clust.df<-data.frame(id=bd.ip$id, x=bd.ip$x, y=bd.ip$y, sim.clust=sim.clust, phylosim.clust=phylosim.clust)

pdf(file=paste0("results/results_figures/simpson_upgma.pdf"),  height=8.5, width=11, useDingbats=FALSE)
par(mfrow=c(1,2))
sim.dend %>% color_branches(col=rev(c("#d7191c", "#fdae61", "#a6d96a", "#1a9641"))[clust.df[labels(sim.dend),"sim.clust"]]) %>% set("labels_cex", 0.6) %>% plot(main="Simpson-UPGMA")

phylosim.dend %>% color_branches(col=rev(c("#d7191c", "#fdae61", "#a6d96a", "#1a9641"))[clust.df[labels(phylosim.dend),"phylosim.clust"]]) %>% set("labels_cex", 0.6) %>% plot(main="PhyloSimpson-UPGMA")
dev.off()

#plot cluster analysis results
clust.df$sim.clust<-factor(clust.df$sim.clust)
clust.df$phylosim.clust<-factor(clust.df$phylosim.clust)

my.fill<-scale_fill_manual(name=NULL, values=rev(c("#d7191c", "#fdae61", "#a6d96a", "#1a9641")), drop=FALSE)

my.scale<-scale_y_continuous(breaks=seq(8, 24, 2))

plot.beta.sim<-list()

plot.beta.sim[[1]]<-ggplot(data=clust.df, aes(x=x, y=y)) + geom_raster(aes(fill=sim.clust)) + ggtitle("Simpson dissimilarity") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill + my.scale + coord_fixed()

plot.beta.sim[[2]]<-ggplot(data=clust.df, aes(x=x, y=y)) + geom_raster(aes(fill=phylosim.clust)) + ggtitle("PhyloSimpson dissimilarity") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill + my.scale + coord_fixed()

pdf(file=paste0("results/results_figures/simpson_map.pdf"), height=8.5, width=11, useDingbats=FALSE)
ggarrange(plot.beta.sim[[1]], plot.beta.sim[[2]], ncol=2, nrow=1)
dev.off()

################################
##plotting beta vs. phylo-beta##
#calculate Phylogenetic Sorensen's index for randomized tree
set.seed(3246)
phylo.sor.rand<-phylosor.rnd(samp=bd.mat, tree=cen.tree, cstSor=TRUE, null.model="taxa.labels", runs=999)

#creating a function to calculate SES-phylosor
ses.sor<-function(obs, rand){
#convert distance object in each list into a matrix - the dist object only has the lower triangle listed - converting it into a matrix creates a symmetric matrix
rand<-lapply(rand, as.matrix)

#align all the random matrices as a matrices
rand<-abind(rand, along=3)

#convert it into a dissimilarity matrix
rand<-1-rand

#calculate mean matrix across all random matrices
sor.mean<-apply(rand, c(1,2), mean)

#calculate sd matrix across all random matrices
sor.sd<-apply(rand, c(1,2), sd)

#find the SES of observed values
sor.obs<-as.matrix(obs)
sor.ses<-(sor.obs-sor.mean)/sor.sd

#convert the upper triangle (including the diagonal) into NA so that it can be removed
sor.ses[upper.tri(sor.ses, diag=TRUE)]<-NA

#converting the above into a dataframe
sor.ses<-melt(sor.ses)
colnames(sor.ses)<-c("site1", "site2", "phylosor.ses")

#remove rows with NA values - checked that the number of rows retained is as expected (990)
sor.ses<-sor.ses[!is.na(sor.ses$phylosor.ses),]

return(sor.ses)
}

#running the ses.sor function with our data
phylo.sor.ses<-ses.sor(phylo.sor, phylo.sor.rand)

#eyeballing the plot of phylo.sor.ses
plot(phylo.sor.ses$phylosor.ses)
abline(h=0, col="red")

#left_join phylosor.sor.ses with sor
phylo.sor.ses<-left_join(phylo.sor.ses, sor, by=c("site1", "site2"))

#classifying the sites as SWG, Nilgiris, CWG, NWG based on the id in the site1 and site2 columns
sw<-c(38:45)
nilg<-c(30:37)
cw<-c(18:29)
nw<-c(1:17)

phylo.sor.ses$site1.cat[phylo.sor.ses$site1 %in% sw]<-"SWG"
phylo.sor.ses$site1.cat[phylo.sor.ses$site1 %in% nilg]<-"Nilgiris"
phylo.sor.ses$site1.cat[phylo.sor.ses$site1 %in% cw]<-"CWG"
phylo.sor.ses$site1.cat[phylo.sor.ses$site1 %in% nw]<-"NWG"

phylo.sor.ses$site2.cat[phylo.sor.ses$site2 %in% sw]<-"SWG"
phylo.sor.ses$site2.cat[phylo.sor.ses$site2 %in% nilg]<-"Nilgiris"
phylo.sor.ses$site2.cat[phylo.sor.ses$site2 %in% cw]<-"CWG"
phylo.sor.ses$site2.cat[phylo.sor.ses$site2 %in% nw]<-"NWG"

#create a new column of the comparison
phylo.sor.ses$site.pair<-paste0(phylo.sor.ses$site1.cat, "-", phylo.sor.ses$site2.cat)
                                
#reordering columns  
phylo.sor.ses<-phylo.sor.ses[,c(1,2,5:7,4,3)]

#converting site.pair into a factor
phylo.sor.ses$site.pair<-factor(phylo.sor.ses$site.pair, levels=c("SWG-SWG", "Nilgiris-Nilgiris", "CWG-CWG", "NWG-NWG", "SWG-Nilgiris", "SWG-CWG", "SWG-NWG", "Nilgiris-CWG", "Nilgiris-NWG", "CWG-NWG"))

#plotting boxplot of phylosor.ses for different comparisons
my.scale<-scale_y_continuous(limits=c(-40, 65), breaks=seq(-60, 60, 20))

pdf(file=paste0("results/results_figures/phylosor_ses_boxplot.pdf"), height=8.5, width=11, useDingbats=FALSE)
ggplot(data=phylo.sor.ses[is.finite(phylo.sor.ses$phylosor.ses),], aes(x=site.pair, y=phylosor.ses)) + geom_boxplot() + xlab("Comparison") + ylab("PhyloSorensen-SES") + geom_hline(aes(yintercept=0)) + theme_bw() + my.scale
dev.off()

#plotting phylosor.ses against sor
my.col<-scale_color_manual(name="Species", values=c('#beaed4','#7fc97f','#fdc086', '#bf5b17', '#984ea3', '#377eb8', '#4daf4a', '#ff7f00', '#e41a1c', '#a65628'), drop=FALSE)
my.shape<-scale_shape_manual(name="Comparisons", values=c(17,17,17,17,19,19,19,19,19,19), drop=FALSE)
my.alpha<-scale_alpha_discrete(name="Comparisons", values=c(1,1,1,1,rep(0.5, 6)))
my.scale<-scale_y_continuous(limits=c(-60, 65), breaks=seq(-60, 60, 20))

pdf(file=paste0("results/results_figures/fig4_phylobeta_vs_beta.pdf"), height=8.5, width=11.5, useDingbats=FALSE)
ggplot(data=phylo.sor.ses[is.finite(phylo.sor.ses$phylosor.ses),], aes(x=sor, y=phylosor.ses)) + geom_point(aes(col=site.pair, shape=site.pair), size=1.5, alpha=0.6) + xlab("Sorensen's dissimilarity index") + ylab("PhyloSor-SES") + geom_hline(aes(yintercept=0)) + my.col + my.shape + my.scale
dev.off()

###



####LOOSE CODE####
#comparing the observed values of beta-diversity with null model - calculating the standardized effect sizes using a function borrowed from https://pedrohbraga.github.io/PhyloCompMethods-in-R-workshop/PhyloCompMethodsMaterial.html#community-phylogenetic-patterns
ses.phylosor<-function(obs, rand){
#the as.data.frame function turns each of the 999 randomized matrices into a column of a dataframe, the transpose function changes it into rows. So there would be 999 rows and 990 columns ((45*45)-45)/2 (considering the lower or upper triangle in the matrix)
rand<-t(as.data.frame(lapply(rand, as.vector)))
#I added this to convert the phylogenetic similarity matrix into a dissimilarity matrix
rand<-1-rand
phySor.obs<-as.numeric(obs)
#find column mean and sd for each cell across randomizations
phySor.mean<-apply(rand, MARGIN=2, FUN=mean, na.rm=TRUE)
phySor.sd<-apply(rand, MARGIN=2, FUN=sd, na.rm = TRUE)
phySor.ses<-(phySor.obs-phySor.mean)/phySor.sd
#rank the observed and randomized values - to get some distribution of values
phySor.obs.rank<-apply(X=rbind(phySor.obs, rand), MARGIN=2, FUN=rank)[1, ]
phySor.obs.rank<-ifelse(is.na(phySor.mean), NA, phySor.obs.rank)
#calculate the p value based on ranks (similar to the percentile analysis)
data.frame(phySor.obs, phySor.mean, phySor.sd, phySor.obs.rank, phySor.ses, phySor.obs.p=phySor.obs.rank/(dim(rand)[1] + 1))
}

##USING RECLUSTER##
library(recluster)

#load biodiverse input file
bd.ip<-read.csv("dataframes/biodiverse/biodiverse_edit_6Jul20.csv", header=TRUE)
bd.mat<-as.matrix(bd.ip[,4:22])
rownames(bd.mat)<-bd.ip$id

#load the phylogenetic tree file
#loading the tree file
cen.tree<-read.tree("dataframes/biodiverse/threegenera_biodiverse.tre")

#calculating beta-diversity 
#first create a consensus UPGMA tree
#tr=number of trees for consensus
#p=proportion for clade to be represented in the consensus tree - not sure if this is only for phylogenetic beta diversity. What is the proportion that they are referring to?, look up dist options from designdist in vegan package - the one I use here is for Sorensen's, method options from hclust for clustering methods - average=UPGMA, blenghts - A logical asking if non-negative least squares branch lengths should be computed.
set.seed(3443)
sor.cons<-recluster.cons(bd.mat, phylo=NULL, tr=100, p=1, dist="(A+B-2*J)/(A+B)", method="average", blenghts=TRUE, select=FALSE)

#compute node bootstraps - for each bootstrap iteration, the row order is resampled several times and a consensus rule is applied. How is this different from recluster.cons? That function also resamples sites
#level referes to the ratio of species numbers to be included in the analysis and number of species in matrix
#the result gives percentage bootstrap values for each node
set.seed(5432)
sor.boot<-recluster.boot(tree=sor.cons$cons, mat=bd.mat, tr=100, p=1, dist="(A+B-2*J)/(A+B)", method="average", boot=1000, level=1)

#plot the consensus tree
recluster.plot(sor.cons$cons,sor.boot,direction="downwards")

#calculating phylogenetic beta diversity
#the recluster package already converts phylosor into a dissimilarity matrix - can use that directly
set.seed(3443)
phylo.sor.cons<-recluster.cons(bd.mat, phylo=cen.tree, tr=100, p=1, dist="phylosor", method="average")

#compute node bootstraps
set.seed(5432)
phylo.sor.boot<-recluster.boot(tree=phylo.sor.cons$cons, mat=bd.mat, tr=100, p=1, method="average", boot=1000, level=1)

recluster.plot(phylo.sor.cons$cons,phylo.sor.boot,direction="downwards")

#removing the lat and lon columns from the ip file
head(bd.ip)
clus.ip<-bd.ip[,4:ncol(bd.ip)]
clus.op<-vegdist(clus.ip, method="bray", binary=TRUE)
clus.tree<-hclust(clus.op, method="average")
clus.dend<-as.dendrogram(clus.tree)
clus.phylo<-as.phylo(clus.tree)

library(dendextend)
pdf(file=paste0("results/results_figures/fig9_beta_diversity_dendrogram.pdf"), height=5, width=8, useDingbats=FALSE)
clus.dend %>% set("branches_k_color", value=c("#d73027", "#fc8d59", "#fee08b", "#d9ef8b", "#91cf60", "#1a9850"), k=6) %>%plot(main="Sorensen's dissimilarity - UPGMA tree")
dev.off()

pdf(file=paste0("results/results_figures/fig9_beta_diversity_dendrogram_guide.pdf"), height=8, width=5, useDingbats=FALSE)
plot.dend<-as.ggdend(clus.dend %>% set("branches_k_color", value=c("#d73027", "#fc8d59", "#fee08b", "#d9ef8b", "#91cf60", "#1a9850"), k=6))
ggplot(plot.dend, horiz=TRUE, labels=FALSE)
dev.off()

plot.dend<-as.ggdend(clus.dend %>% set("branches_k_color", value=c("#d73027", "#fc8d59", "#fee08b", "#d9ef8b", "#91cf60", "#1a9850"), k=6))
ggplot(plot.dend, horiz=TRUE, labels=FALSE)

#matching the same names on the lat longs of dataframe
clust<-cutree(clus.tree, k=6)                    
clust.df<-data.frame(id=bd.ip$id, x=bd.ip$x, y=bd.ip$y, cluster=factor(clust))
clust.df$cluster2<-clust.df$cluster
levels(clust.df$cluster2)<-c("3", "2", "1", "6", "4", "5")
clust.df$cluster2<-factor(clust.df$cluster2, levels=c("1", "2", "3", "4", "5", "6"))

#the way the clusters are cut in the dendextend package are different from how they are cut using the cutree function - the colours wouldn't match - I need to match it at the level of the ggplot function 
#plot cluster analysis results
#my.fill<-scale_fill_manual(name=NULL, values=c('#fee08b', '#fc8d59', '#d73027', '#1a9850', '#d9ef8b', '#91cf60'), drop=FALSE)
my.fill<-scale_fill_manual(name=NULL, values=c("#d73027", "#fc8d59", "#fee08b", "#d9ef8b", "#91cf60", "#1a9850"), drop=FALSE)

my.scale<-scale_y_continuous(breaks=seq(8, 24, 2))

pdf(file=paste0("results/results_figures/fig8_beta_diversity.pdf"), height=8, width=5, useDingbats=FALSE)
plot.beta<-ggplot(data=clust.df, aes(x=x, y=y)) + geom_raster(aes(fill=cluster2)) + ggtitle("Sorensen's dissimilarity") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill + my.scale + coord_fixed()
dev.off()

#arranging the above graphs
pdf(file=paste0("results/results_figures/fig8_beta_diversity_cluster.pdf"),  height=15, width=12, useDingbats=FALSE)
ggarrange(plot.fin, ncol=2, nrow=1)
dev.off()

##plot three genera phylogenetic tree
#creating the input tree file for biodiverse
#loading the tree file
cen.tree<-phylo4(read.tree("dataframes/biodiverse/threegenera.tre"))

#inputing the tip label remap code
remap<-read.csv("dataframes/tiplabels.csv", header=FALSE)
colnames(remap)<-c("tip_name", "sp_name")

#removing extra quotes from tip labels
remap$tip_name<-str_replace_all(remap$tip_name, "\"", "")

#giving the final name
remap$final_name<-c(rep("D. jonesii", 2), rep("D. coonoorensis", 2), "D. jangii", "D. nudus", rep("D. barnabasi", 3), "R. pazhuthara", "R. crassispina", "R. konda", "R. longipes", "R. sp. 2", "R. lewisi", "R. sp. 1", "R. trispinosa", "R. sada", "R. aspinosa", "R. ikhalama", "R. immarginata", "E. agasthyamalaiensis", "E. coonooranus", "E. sahyadrensis", "E. praveeni", "E. tristis")

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

#plot tree
my.col<-c(rep("blue", times=5), rep("red", times=9), rep("black", times=5))

pdf(file=paste0("results/results_figures/fig1_phylogenetic_tree.pdf"), height=5, width=5, useDingbats=FALSE)
plot(cen.tree2,no.margin=TRUE,edge.width=2, tip.color=my.col)
add.scale.bar()
dev.off()

test<-biodiv[biodiv$Genus=="Digitipes",]
coordinates(test)<-~x+y
crs(test)<-wgs
plot(test)

#convert tree object to phylo format and save to disk
#cen.tree2<-as(cen.tree2, "phylo")
labels(cen.tree2)
my.col<-c(rep("blue", times=5), rep("red", times=9), rep("black", times=5))
ggtree(cen.tree2) + theme_tree2() + geom_text(aes(label=label), hjust=-0.1, size=3)), color=my.col)

ggplot(cen.tree3, aes(x, y)) + geom_tree() + theme_tree()

########
#plotting beta-SES by index
my.col<-scale_color_manual(name="Comparisons", values=c('#045a8d','#2b8cbe','#bdc9e1', '#d53e4f', '#f46d43', '#fee08b'), drop=FALSE)
my.shape<-scale_shape_manual(name="Comparisons", values=c(4,4,4,19,19,19))
ggplot(data=phylo.sor.ses[is.finite(phylo.sor.ses$phyloSor.ses),], aes(y=phyloSor.ses, x=as.numeric(rownames(phylo.sor.ses[is.finite(phylo.sor.ses$phyloSor.ses),])))) + geom_point(aes(col=site.pair, shape=site.pair)) + ggtitle("PhyloSor-SES") + my.col + my.shape


