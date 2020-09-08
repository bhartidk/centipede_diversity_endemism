###This code takes the Maxent habitat suitability predictions at 1km resolution (not thresholded, proportional to probability of presence), stacks them across species within a genus and uses a genus level phylogenetic tree to calculate diversity and endemism indices. These indices are finally summed across all genera. The final output is raster maps of diversity and endemism measures at 1km resolution.

#set working directory - this should lead to the location where the folders /results, /dataframes and /R scripts folders are present
#setwd("/home/bharti/D/PostDoc data/spatial_data")

library(reshape2)
library(raster)
library(phylobase)
library(ape)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(abind)
library(rgdal)

#setting projection
wgs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#loading the western ghats shape file
wg<-readOGR(dsn='dataframes/sdm_layers/wg_vector', layer='wg_boundary_mod_23Jun20')

#loading the gps.pts dataframe
gps<-read.csv("dataframes/gps_9Jun20.csv")

#creating a new column with Species.Name split into Genus and Species columns
new.cols<-colsplit(gps$Species.Name," ",c("Genus","Species"))
gps<-cbind(new.cols, gps[,c(1,3,4)])

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

#change Species.Names to match tip labels
gps$Species.Name<-as.character(gps$Species.Name)
gps$Species.Name[gps$Species.Name=="Digitipes barnabasi"]<-"Digitipes barnabasi 1"
gps$Species.Name[gps$Species.Name=="Digitipes barnabasi - 2"]<-"Digitipes barnabasi 2"
gps$Species.Name[gps$Species.Name=="Digitipes barnabasi - 3"]<-"Digitipes barnabasi 3"
gps$Species.Name[gps$Species.Name=="Rhysida sp. 1"]<-"Rhysida sp 1"

#creating a sp object with species names
gps$Species.Name2<-paste0(gps$Genus, "_", gps$Species)

#match Species.Name2 with Species.Name
remap.sp<-unique(gps[,c(3,6)])
colnames(remap.sp)<-c("sp_name", "final_name")

#create a vector of folder names where the sdm results have already been stored
sp<-unique(gps$Species.Name2)

#loading the calc_div_end_22May20.R which is a combination of calc_PE.R script and phylogenetic endemism.r scripts from the github page of Dan Rosauer - https://github.com/DanRosauer/phylospatial
source("R scripts/calc_div_end_22May20.R")

#diff bgs
bg.names<-c("randomBg", "biasBg")

#since we decided to use the results only at the level of Western Ghats, cropping the prediction raster to the western ghats extent early on itself

for(h in 1:length(bg.names)){
#create a list to store the raster values
spmat<-list()

for(i in 1:length(sp)){
#assigning folder name
fold<-paste0("results/", sp[i], "/")
#load the Maxent output of the ith species
sdm.file<-list.files(path=fold, pattern=paste0("*", bg.names[h], "_map.tif"))
sdm<-raster(paste0(fold, sdm.file))

#cropping the prediction raster to the extent of Western Ghats
sdm<-mask(sdm, wg)

#convert the raster layer into a linear vector of raster values
sdm.vec<-rasterToPoints(sdm)

#checking to see if there are any NA values
print(paste0("The no. of NA values is ", length(which(is.na(sdm.vec[,3])))))

#removing NA values - which would be the mask and assigning this to the list
spmat[[i]]<-sdm.vec[,3]
}

#cbinding the list
spmat<-as.data.frame(Reduce(cbind, spmat))

#giving this dataframe column names equivalent to the species from where these values were derived
colnames(spmat)<-sp

#subset each matrix by genus (all species belonging to a single genus) in one matrix - the Phylogenetic Endemism values will be calculated genus wise and summed across genera (following Fenker et al, 2020)
rhy.mat<-spmat[,grep("Rhysida", colnames(spmat))]
eth.mat<-spmat[,grep("Ethmostigmus", colnames(spmat))]
dig.mat<-spmat[,grep("Digitipes", colnames(spmat))]

#loading the tree file
rhy.tree<-phylo4(read.tree("dataframes/Rhysida.tre"))
eth.tree<-phylo4(read.tree("dataframes/Ethmo.tre"))
dig.tree<-phylo4(read.tree("dataframes/Digitipes.tre"))

#inputing the tip label remap code
remap<-read.csv("dataframes/tiplabels.csv", header=FALSE)
colnames(remap)<-c("tip_name", "sp_name")

#removing extra quotes from tip labels
remap$tip_name<-str_replace_all(remap$tip_name, "\"", "")

#convert from factor to character
remap[,1]<-as.character(remap[,1])
remap[,2]<-as.character(remap[,2])

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
rhy.tree2<-rename.tips(rhy.tree)
eth.tree2<-rename.tips(eth.tree)
dig.tree2<-rename.tips(dig.tree)

#remove within species edges - keep the longest branch from an intra-specific clade
rhy.tree2<-subset(rhy.tree2, tips.exclude=c("Rhysida sp 2", "Rhysida ikhalama", "Rhysida immarginata"))

dig.tree2<-subset(dig.tree2, tips.exclude=c("Digitipes jonesii 2", "Digitipes coonoorensis - 1", "Digitipes barnabasi 2", "Digitipes barnabasi 3"))

#rename the tip labels with final_name
rename.tips.fin<-function(subtree){
for(i in 1:nTips(subtree)) {
tip_label<-labels(subtree)[i]
new_label<-remap.sp$final_name[remap.sp$sp_name==tip_label]
labels(subtree)[i]<-new_label
}
return(subtree)
}

#renaming the tip labels again so that they match the names from spatial data
rhy.tree2<-rename.tips.fin(rhy.tree2)
eth.tree2<-rename.tips.fin(eth.tree2)
dig.tree2<-rename.tips.fin(dig.tree2)

#ensure that the tree tips match the spatial names
check<-function(mat, tree){
spatial_names<-colnames(mat)
labels<-as.character(tipLabels(tree))
on_tree<-intersect(spatial_names,labels)
not_on_tree<-setdiff(spatial_names,labels)
if(length(not_on_tree)==0){
mat<-mat[,labels]
print("All the spatial data names are present in the phylogenetic tree")
return(mat)
}else{
print("Check spatial labels and tip labels")
return(mat)
}
}

#checking if labels match and changing the order of the columns in the spatial matrix to match the order of labels in the tree
rhy.mat2<-check(rhy.mat, rhy.tree2)
eth.mat2<-check(eth.mat, eth.tree2)
dig.mat2<-check(dig.mat, dig.tree2)

#make a list of mat and tree files from these species
mat<-list(rhy.mat2, eth.mat2, dig.mat2)
tree<-list(rhy.tree2, eth.tree2, dig.tree2)

#initialize lists where output of the functions will be saved
sr<-matrix(0, nrow=nrow(mat[[1]]), ncol=length(mat))
we<-matrix(0, nrow=nrow(mat[[1]]), ncol=length(mat))
pd<-matrix(0, nrow=nrow(mat[[1]]), ncol=length(mat))
pe<-matrix(0, nrow=nrow(mat[[1]]), ncol=length(mat))

#calculating diversity and endemism measures genus-wise
for(i in 1:length(mat)){
sr[,i]<-calc_SR(tree[[i]], mat[[i]], "probability")
we[,i]<-calc_WE(tree[[i]], mat[[i]], "probability")
pd[,i]<-calc_PD(tree[[i]], mat[[i]], "probability")
pe[,i]<-calc_PE(tree[[i]], mat[[i]], "probability")
}

#summing these measures across genera - converting it into a list
div<-list(rowSums(sr), rowSums(we), rowSums(pd), rowSums(pe))

#create function to scale by row sums 
prop.r<-function(x){
return(x/sum(x))
}

#apply to each of the diversity and endemism measures
div.prop<-list(t(apply(sr, 1, prop.r)), t(apply(we, 1, prop.r)), t(apply(pd, 1, prop.r)), t(apply(pe, 1, prop.r)))

div.names<-c("SR", "WE", "PD", "PE")

#getting the longitude and longitude of the raster cells which have non-NA values - confirmed that the order is the same
coord<-rasterToPoints(mask(raster(paste0(fold, sdm.file)), wg))
coord<-coord[,c(1,2)]

#running a loop for plotting diversity measures
for(i in 1:length(div)){
#cbinding the coordinates to the result
colnames(ras)<-c("Longitude", "Latitude", "div")

rhy.prop<-as.data.frame(cbind(coord, div.prop[[i]][,1]))
eth.prop<-as.data.frame(cbind(coord, div.prop[[i]][,2]))
dig.prop<-as.data.frame(cbind(coord, div.prop[[i]][,3]))

colnames(rhy.prop)<-c("Longitude", "Latitude", "div.prop")
colnames(eth.prop)<-colnames(rhy.prop)
colnames(dig.prop)<-colnames(dig.prop)

#write the data to disk
writeRaster(rasterFromXYZ(ras), file=paste0("results/", div.names[i]), format="GTiff", overwrite=TRUE)

writeRaster(rasterFromXYZ(rhy.prop), file=paste0("results/Rhysida_", div.names[i], "_prop"), format="GTiff", overwrite=TRUE)

writeRaster(rasterFromXYZ(eth.prop), file=paste0("results/Ethmostigmus_", div.names[i], "_prop"), format="GTiff", overwrite=TRUE)

writeRaster(rasterFromXYZ(dig.prop), file=paste0("results/Digitipes_", div.names[i], "_prop"), format="GTiff", overwrite=TRUE)
}
}

### LOOSE CODE ###
#pdf(file=paste0("results/pe_", bg.names[h], "_log_cs2.pdf"),  height=5, width=8)
#map2<-ggplot(data=ras.pe, aes(x=Longitude, y=Latitude)) + geom_raster(aes(fill=log(PE))) + scale_fill_gradient2(high="red", low="white", midpoint=median(log(ras.pe$PE)))
#print(map2)
#dev.off()

#pdf(file=paste0("results/pe_", bg.names[h], "_cs2.pdf"),  height=5, width=8)
#map4<-ggplot(data=ras.pe, aes(x=Longitude, y=Latitude)) + geom_raster(aes(fill=PE)) + scale_fill_gradient2(high="red", low="white", midpoint=median(ras.pe$PE))
#print(map4)
#dev.off()

### optional plot PE result ###
# 
# library(classInt)
# class_count <- 15
# 
# my.class.fr<-classIntervals(log(output$PE),n=class_count,style="equal")   
# my.col.fr<-findColours(my.class.fr,rainbow(class_count,start=0.1)) # ramp colors
# 
# # Map
# x11()
# plot(output$Longitude,output$Latitude,col=my.col.fr,pch=20, xlab="Longitude", ylab="Latitude", main="Phylogenetic Endemism (PE) for Rhysida in peninsular India")

#my.fill1<-scale_fill_gradientn(name=div.names[i], colours=rev(c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')))

#map[[i]]<-ggplot(data=ras, aes(x=Longitude, y=Latitude)) + geom_raster(aes(fill=div)) + my.fill1 + ggtitle(paste0(bg.names[h], " - ", div.names[i])) + theme_bw()

#my.fill2<-scale_fill_gradientn(name=paste0("log(",div.names[i], ")"), colours=rev(c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')))

#log.map[[i]]<-ggplot(data=ras, aes(x=Longitude, y=Latitude)) + geom_raster(aes(fill=log(div))) + my.fill2 + ggtitle(paste0(bg.names[h], " - log(", div.names[i], ")")) + theme_bw()
#}

#pdf(file=paste0("results/", bg.names[h], "_div.pdf"),  height=10, width=15)
#panel.map<-ggarrange(map[[1]], log.map[[2]], map[[3]], log.map[[4]],  ncol=2, nrow=2)
#print(panel.map)
#dev.off()

#pdf(file=paste0("results/", bg.names[h], "_div.pdf"),  height=10, width=15)
#ggarrange(map[[1]], map[[2]], map[[3]], map[[4]], ncol=2, nrow=2)
#dev.off()    
#}

#plotting the residuals from linear regression of PD~SR and PE~WE
# #looking at a plot of PD~SR
# sr.all<-div[[1]]
# pd.all<-div[[3]]
# 
# plot(x=sr.all, y=pd.all)
# lines(x=sr.all, y=pd.predict, col="red")
# 
# #running the linear regression
# pd.lm<-lm(pd.all~(sr.all+(sr.all)^2))
# summary(pd.lm)
# 
# pd.predict<-predict(pd.lm)
# 
# #plot standardized residuals against SR
# plot(resid(pd.lm)~fitted(pd.lm))

