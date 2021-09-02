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
library(RStoolbox)

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

#######################################
####ANALYZE BIODIVERSE OUTPUT FILES####
#load biodiverse input file
bd.ip<-read.csv("dataframes/biodiverse/biodiverse_16Mar21.csv", header=TRUE)

bd<-read.csv("results/biodiverse/threegenera_results_16Mar21.csv", header=TRUE)
bd.r<-read.csv("results/biodiverse/threegenera_struc_randomization_results_16Mar21.csv", header=TRUE)

#keeping the columns of interest from bd
colnames(bd)
identical(bd$ENDW_RICHNESS, bd$RICHNESS_ALL) #identical - keeping only richness
identical(bd$RICHNESS_ALL, bd$RICHNESS_SET1) #identical - keeping only richness
identical(bd$ENDW_CWE, bd$ENDW_WE) #not the same, ENDW_CWE is weighted by richness - but I don't need this, keeping only ENDW_WE
identical(bd$ENDW_WE, bd$ENDW_SINGLE) #these are the same, retaining only ENDW_WE
#_P is scaled by tree length, _per_taxon is proportion contribution by each taxon 0 _RPD_DIFF2 tells us how much more or less PD is there than expected, _RPD_NULL - gives the value of the denominator used for RPD calculations

#based on the above reference, only keeping columns that are of interest to us
bd<-bd[,c(2,19,6,7,8,11,12,13,15,16,18)]

#renaming columns
colnames(bd)
colnames(bd)<-c("id", "sr", "we", "pd", "pd.p", "pe", "pe.p", "rpd", "rpd.denom", "rpe", "rpe.denom")

#keeping the column names of interest from bd.r
colnames(bd.r) #what we are interested in is the percentiles - ranking of original value against randomizations
bd.r<-bd.r[,c(2,25,26,29,30,31,33,34,36)]
colnames(bd.r)
colnames(bd.r)<-c("id", "p.pd", "p.pd.p", "p.pe", "p.pe.p", "p.rpd", "p.rpd.denom", "p.rpe", "p.rpe.denom") 

#left_join bd.ip and bd
bd.all<-left_join(bd.ip, bd, by="id")
bd.all<-left_join(bd.all, bd.r, by="id")

#eyeballing the distribution of various values
#par(mfrow=c(2,2))
hist((bd.all$p.pe)*100, xlab="Percentile", main="Histogram of PE percentiles")
hist((bd.all$p.rpe)*100, xlab="Percentile", main="Histogram of RPE percentiles")
hist((bd.all$p.pd)*100, xlab="Percentile", main="Histogram of PD percentiles")
hist((bd.all$p.rpd)*100, xlab="Percentile", main="Histogram of RPD percentiles")

#create a new column which tells us if the cells are either greater than 90 percentile or less than 10 percentile
#the break mentioned - anything less than that is put in one category
bd.all$pd.sig<-cut(bd.all$p.pd.p, breaks=c(-0.00001, 0.2, 0.8, 1.00001), labels=c("<20", "20-80", ">80"))
bd.all$pe.sig<-cut(bd.all$p.pe.p, breaks=c(-0.00001, 0.2, 0.8, 1.00001), labels=c("<20", "20-80", ">80"))
bd.all$rpd.sig<-cut(bd.all$p.rpd, breaks=c(-0.00001, 0.2, 0.8, 1.00001), labels=c("<20", "20-80", ">80"))
bd.all$rpe.sig<-cut(bd.all$p.rpe, breaks=c(-0.00001, 0.2, 0.8, 1.00001), labels=c("<20", "20-80", ">80"))

#plot sr, pd, we, pe
#normal scale
my.fill1<-scale_fill_gradientn(name=NULL, colours=rev(c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')))
my.col.sig<-scale_color_manual(name="Percentile", labels=c("<20", ">80"), values=c("#bdbdbd", "#000000"), drop=TRUE)
my.scale<-scale_y_continuous(breaks=seq(8, 24, 2))

plot.raw<-list()

#plot diversity measures for western ghats alone
plot.raw[[1]]<-ggplot(data=bd.all, aes(x=x, y=y)) + geom_raster(aes(fill=sr)) + my.fill1 + ggtitle("TD") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.scale + coord_fixed()

plot.raw[[2]]<-ggplot(data=bd.all, aes(x=x, y=y)) + geom_raster(aes(fill=we)) + my.fill1 + ggtitle("WE") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.scale + coord_fixed()

plot.raw[[3]]<-ggplot(data=bd.all, aes(x=x, y=y)) + geom_raster(aes(fill=pd.p)) + my.fill1 + ggtitle("PD") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.scale + coord_fixed() + geom_point(data=bd.all[bd.all$pd.sig!="20-80",], aes(x=x, y=y, color=pd.sig, size=2)) + my.col.sig

plot.raw[[4]]<-ggplot(data=bd.all, aes(x=x, y=y)) + geom_raster(aes(fill=pe.p)) + my.fill1 + ggtitle("PE") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.scale + coord_fixed() + geom_point(data=bd.all[bd.all$pe.sig!="20-80",], aes(x=x, y=y, color=pe.sig, size=2)) + my.col.sig

plot.raw[[5]]<-ggplot(data=bd.all, aes(x=x, y=y)) + geom_raster(aes(fill=rpd)) + my.fill1 + ggtitle("RPD") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.scale + coord_fixed() + geom_point(data=bd.all[bd.all$rpd.sig!="20-80",], aes(x=x, y=y, color=rpd.sig, size=2)) + my.col.sig

plot.raw[[6]]<-ggplot(data=bd.all, aes(x=x, y=y)) + geom_raster(aes(fill=rpe)) + my.fill1 + ggtitle("RPE") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.scale + coord_fixed() + geom_point(data=bd.all[bd.all$rpe.sig!="20-80",], aes(x=x, y=y, color=rpe.sig, size=2)) + my.col.sig

for(i in 1:length(plot.raw)){
jpeg(file=paste0("results/revisions/results_figures/div_end_", i, ".jpeg"), height=8, width=5, res=300, units="in")
print(plot.raw[[i]])
dev.off()

pdf(file=paste0("results/revisions/results_figures/div_end_", i, ".pdf"), height=11, width=8.5, useDingbats=FALSE)
print(plot.raw[[i]])
dev.off()
}

##################
##Beta diversity##
#I was not able to load the cluster analysis results from Biodiverse properly, so redoing the analysis here
#load biodiverse input file
bd.ip<-read.csv("dataframes/biodiverse/biodiverse_16Mar21.csv", header=TRUE)
bd.mat<-as.matrix(bd.ip[,4:22])
cen.tree<-read.tree("dataframes/biodiverse/threegenera_biodiverse.tre")

#checking if species names are identical in the matrix and the tree
identical(sort(colnames(bd.mat)), sort(labels(cen.tree)))

##SIMPSON INDEX##
beta<-betapart::beta.pair(bd.mat, index.family="sorensen")
sim<-beta$beta.sim
sim.tree<-as.dendrogram(hclust(sim, method="average"))
sim.dend<-as.dendrogram(sim.tree)

phylobeta<-betapart::phylo.beta.pair(bd.mat, cen.tree, index.family="sorensen")
phylosim<-phylobeta$phylo.beta.sim
phylosim.tree<-as.dendrogram(hclust(phylosim, method="average"))
phylosim.dend<-as.dendrogram(phylosim.tree)

#matching the same names on the lat longs of dataframe
sim.clust<-cutree(sim.tree, k=4) 
phylosim.clust<-cutree(phylosim.tree, k=4)
clust.df<-data.frame(id=bd.ip$id, x=bd.ip$x, y=bd.ip$y, sim.clust=sim.clust, phylosim.clust=phylosim.clust)

pdf(file=paste0("results/revisions/results_figures/simpson_upgma_trees.pdf"),  height=8.5, width=11, useDingbats=FALSE)
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

plot.beta.sim[[2]]<-ggplot(data=clust.df, aes(x=x, y=y)) + geom_raster(aes(fill=phylosim.clust)) + ggtitle("PhyloSor(Turn)") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.fill + my.scale + coord_fixed()

pdf(file=paste0("results/revisions/results_figures/simpson_map.pdf"), height=8.5, width=11, useDingbats=FALSE)
ggarrange(plot.beta.sim[[1]], plot.beta.sim[[2]], ncol=2, nrow=1)
dev.off()

for(i in 1:length(plot.beta.sim)){
jpeg(file=paste0("results/revisions/results_figures/beta_div_", i, ".jpeg"), height=8, width=5, res=300, units="in")
print(plot.beta.sim[[i]])
dev.off()
}

#visualizing pairwise dissimilarity on a map - from the gdm vignette
#perform principal component analysis using the simpson/phylosimpson dissimilarity matrix, to split the dissimilarity based on multiple dimensions
pca.sim<-prcomp(sim)
pca.sim.x<-as.data.frame(pca.sim$x)

#add location information to the above dataframe
identical(rownames(pca.sim.x), as.character(bd.ip$id))
pca.sim.x<-cbind(bd.ip[,1:3], pca.sim.x[1:3])

#use the PC axes to generate rasters 
pc1<-rasterFromXYZ(pca.sim.x[,c("x", "y", "PC1")], crs=wgs)
pc2<-rasterFromXYZ(pca.sim.x[,c("x", "y", "PC2")], crs=wgs)
pc3<-rasterFromXYZ(pca.sim.x[,c("x", "y", "PC3")], crs=wgs)  
pcaRast<-stack(pc1, pc2, pc3) 

#plot the results in colour space
pcaRast[[1]]<-(pcaRast[[1]]-pcaRast[[1]]@data@min) /
(pcaRast[[1]]@data@max-pcaRast[[1]]@data@min)*255
pcaRast[[2]]<-(pcaRast[[2]]-pcaRast[[2]]@data@min) /
(pcaRast[[2]]@data@max-pcaRast[[2]]@data@min)*255
pcaRast[[3]]<-(pcaRast[[3]]-pcaRast[[3]]@data@min) /
(pcaRast[[3]]@data@max-pcaRast[[3]]@data@min)*255
plotRGB(pcaRast, r=1, g=3, b=2, axes=TRUE)

my.scale<-scale_y_continuous(breaks=seq(8, 24, 2))

pdf(file=paste0("results/revisions/results_figures/simpson_map_rgb.pdf"), height=8, width=5, useDingbats=FALSE)
ggRGB(pcaRast, r=1, g=3, b=2) + ggtitle("Simpson") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.scale + coord_fixed()
dev.off()

rm(pcaRast)

#visualizing pairwise dissimilarity on a map - from the gdm vignette
#perform principal component analysis using the simpson/phylosimpson dissimilarity matrix, to split the dissimilarity based on multiple dimensions
pca.phylosim<-prcomp(phylosim)
pca.phylosim.x<-as.data.frame(pca.phylosim$x)

#add location information to the above dataframe
identical(rownames(pca.phylosim.x), as.character(bd.ip$id))
pca.phylosim.x<-cbind(bd.ip[,1:3], pca.phylosim.x[1:3])

#use the PC axes to generate rasters 
pc1<-rasterFromXYZ(pca.phylosim.x[,c("x", "y", "PC1")], crs=wgs)
pc2<-rasterFromXYZ(pca.phylosim.x[,c("x", "y", "PC2")], crs=wgs)
pc3<-rasterFromXYZ(pca.phylosim.x[,c("x", "y", "PC3")], crs=wgs)  
pcaRast<-stack(pc1, pc2, pc3) 

#plot the results in colour space
pcaRast[[1]]<-(pcaRast[[1]]-pcaRast[[1]]@data@min) /
(pcaRast[[1]]@data@max-pcaRast[[1]]@data@min)*255
pcaRast[[2]]<-(pcaRast[[2]]-pcaRast[[2]]@data@min) /
(pcaRast[[2]]@data@max-pcaRast[[2]]@data@min)*255
pcaRast[[3]]<-(pcaRast[[3]]-pcaRast[[3]]@data@min) /
(pcaRast[[3]]@data@max-pcaRast[[3]]@data@min)*255
plotRGB(pcaRast, r=1, g=3, b=2, axes=TRUE)

my.scale<-scale_y_continuous(breaks=seq(8, 24, 2))

pdf(file=paste0("results/revisions/results_figures/phylosor_turn_map_rgb.pdf"), height=8, width=5, useDingbats=FALSE)
ggRGB(pcaRast, r=1, g=3, b=2) + ggtitle("PhyloSor(Turn)") + geom_polygon(data=wg.map, aes(x=long, y=lat, group=group), size=0.1, color="black", fill=NA) + theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12), axis.title.x=element_blank(), axis.title.y=element_blank()) + my.scale + coord_fixed()
dev.off()

####