#export LD_LIBRARY_PATH=/usr/local/lib
#R CMD BATCH filename.r
#setwd("~/Downloads/maxent_analysis")

#set working directory
#setwd("/home/bharti/D/PostDoc data/centipede_spatial")

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

#since we decided to use the results only at the level of Western Ghats, cropping the prediction raster to the western ghats extent early on itself

#create a list to store the raster values
spmat<-data.frame()

for(i in 1:length(sp)){
#assigning folder name
fold<-paste0("results/revisions/", sp[i], "/processed/")
#load the Maxent output of the ith species
sdm.file<-list.files(path=fold, pattern=paste0("*", "_bc2_rank_map.tif"))
sdm<-raster(paste0(fold, sdm.file))

#cropping the prediction raster to the extent of Western Ghats
sdm<-mask(sdm, wg)

#convert the raster layer into a linear vector of raster values
sdm.vec<-rasterToPoints(sdm)

#checking to see if there are any NA values
print(paste0("The no. of NA values is ", length(which(is.na(sdm.vec[,3])))))

#removing NA values - which would be the mask and assigning this to the list
if(i==1){
spmat<-sdm.vec[,3]
}else{
spmat<-cbind(spmat, sdm.vec[,3])
}
}

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
div.names<-c("SR", "WE", "PD", "PE")

#getting the longitude and latitude of the raster cells which have non-NA values - confirmed that the order is the same
coord<-rasterToPoints(mask(raster(paste0(fold, sdm.file)), wg))
coord<-coord[,c(1,2)]

#running a loop for plotting diversity measures
for(i in 1:length(div)){
#cbinding the coordinates to the result
ras<-as.data.frame(cbind(coord, div[[i]]))
colnames(ras)<-c("Longitude", "Latitude", "div")

#write the data to disk
writeRaster(rasterFromXYZ(ras), file=paste0("results/revisions/results_figures/", div.names[i]), format="GTiff", overwrite=TRUE)
}

div<-list()
map.wg<-list()
log.map<-list()

#load the data files 
for(i in 1:length(div.names)){
div[[i]]<-raster(paste0("results/revisions/results_figures/", div.names[i], ".tif"))
ras<-div[[i]]
ras<-as.data.frame(rasterToPoints(ras))
colnames(ras)<-c("Longitude", "Latitude", "div")

if(div.names[i] %in% c("SR", "PD")){
#normal scale
my.fill1<-scale_fill_gradientn(name=div.names[i], colours=rev(c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')))

my.scale<-scale_y_continuous(breaks=seq(8, 24, 2))

#plot for western ghats alone
map.wg[[i]]<-ggplot(data=ras, aes(x=Longitude, y=Latitude)) + geom_tile(aes(fill=div)) + my.fill1 + ggtitle(paste0(div.names[i])) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.scale + coord_fixed()

}else{

#log scale
my.fill2<-scale_fill_gradientn(name=paste0("log(",div.names[i], ")"), colours=rev(c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')))

#plot for western ghats alone
map.wg[[i]]<-ggplot(data=ras, aes(x=Longitude, y=Latitude)) + geom_tile(aes(fill=log(div))) + my.fill2 + ggtitle(paste0("log(",div.names[i],")")) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.scale + coord_fixed()
}
}

#arranging the above graphs
pdf(file=paste0("results/revisions/results_figures/diversity_measures_1km_continuous.pdf"),  height=15, width=12, useDingbats=FALSE)
ggarrange(map.wg[[1]], map.wg[[3]], map.wg[[2]], map.wg[[4]], ncol=2, nrow=2)
dev.off()

####