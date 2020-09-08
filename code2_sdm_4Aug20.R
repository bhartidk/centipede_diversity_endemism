###This code takes input locations for each species and all available environmental layers to run two Maxent models - one using background locations from the bias model and the second using randomly chosen background locations. The results for each species are automatically created within a subfolder named by species within the /results folder. AUC, TSS and Kappa calculations, five-fold cross-validation (for species with >5 presence locations) and Jackknifing test for small sample sizes are carried out and results saved within the same sub-folder.

#set working directory - this should lead to the location where the folders /results, /dataframes and /R scripts folders are present
#setwd("/home/bharti/D/PostDoc data/spatial_data")

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
#library(rmaxent)
library(MASS)
#library(combinat)
#library(ROCR)

#setting projection
wgs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#mercator projection
epsg.3395<-"+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

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
pred<-stack(wc_layers, elev, soil.stype)

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
#write.csv(ind.occ[[i]], file=paste0("dataframes/", sp2[i], "_occ.csv"), row.names=FALSE)
}

#input all the Bg data
randomBg<-read.csv("dataframes/randomBg.csv")
biasBg<-read.csv("dataframes/biasBg.csv")
#bufBg<-read.csv("dataframes/bufBg.csv")
bg<-list(randomBg, biasBg)
bg.names<-c("randomBg", "biasBg")

#function for cross-validation
#dividing presence data into 5 groups - using 4/5ths for training and 1/5ths for test. Code from Hijmans and Elith 2017

#For each of the folds , dividing all data into test and training data, combining test data with background data and giving presence codes to the rows. 

#Here k refers to the number of k folds, fold refers to the k fold group, pres is the presence data dataframe and bg is the background data dataframe
maxent_kfold<-function(kf, pres, bg){
evl<-as.data.frame(matrix(nrow=5, ncol=8))
colnames(evl)<-c("Threshold", "Sensitivity", "Specificity", "AUC-train", "AUC-test", "Cohen's Kappa", "TSS", "COR")
group<-kfold(pres, kf)
for(k in 1:kf){
test<-pres[group==k,]
train_pres<-pres[group!=k, ]
train<-rbind(train_pres, bg)
code<-c(rep(1, nrow(train_pres)), rep(0, nrow(bg)))
mod.it<-maxent(x=train,  p=code, factors=c('soil_stype', 'geol_lith'))
evl.it<-evaluate(p=test, a=bg, mod.it)

#saving the AUC value
auc.train<-mod.it@results[grep("Training.AUC*", rownames(mod.it@results))]
auc.test<-evl.it@auc

#finding the threshold which produces Maximum Sum of Sensitivity and Specificity
th<-threshold(evl.it, stat='spec_sens')

#extracting kappa value for threshold with Maximum Sum of Sensitivity and Specificity
kappa<-evl.it@kappa[evl.it@t==th]

#calculating TSS value for threshold with Maximum Sum of Sensitivity and Specificity, the Maxent manual says that sensitivity is also called True Positive Rate (TPR) and specificity is also called True Negative Rate (TNR). TSS is TPR+TNR-1
tpr<-evl.it@TPR[evl.it@t==th]
tnr<-evl.it@TNR[evl.it@t==th]
tss<-tpr+tnr-1

#saving the correlation value
cor.val<-evl.it@cor

#saving all these evaluation metrics to disk
evl[k,]<-c(th, tpr, tnr, auc.train, auc.test, kappa, tss, cor.val)
}
evl
}

##MAXENT MODELS WITH CENTIPEDE DATA##
#running a loop perform the model exercise over species and multiple background locations
for(i in 1:length(ind.occ)){
#create a directory for the species in question
dir.create(paste0("results/", unique(gps$Genus[gps$Species==sp[i]]), "_", sp[i]))
fold<-paste0("results/", unique(gps$Genus[gps$Species==sp[i]]), "_", sp[i])

#taking in occurrence values by species
sp.occ<-as.data.frame(ind.occ[i])

#removing multiple sampling locations in a grid cell
source('R scripts/Eliminate Points in Same Cell of a Raster.r')
sp.occ<-elimCellDups(sp.occ, pred[[1]], longLatFields=c('Longitude', 'Latitude'), priority=NULL)
nrow(sp.occ) 

#converting into a SpatialPoints dataframe
coordinates(sp.occ)<-~Longitude+Latitude
crs(sp.occ)<-wgs

#extract values from pred
spSites<-extract(pred, sp.occ)
spSites<-cbind(coordinates(sp.occ), spSites)
spSites<-spSites[complete.cases(spSites),]
write.csv(spSites, file=paste0(fold, "/", sp2[i], "Sites.csv"), row.names=FALSE)

#create a list where model outputs will be saved for each species
mods<-list()
modMaps<-stack()

for(j in 1:length(bg)){
#selecting the kind of background sites
spBg<-bg[[j]]

#removing geol_lith
#spBg<-spBg[,-ncol(spBg)]

#match the presence data with background data
spPres<-spSites
preds<-pred

#creating input data
trainData<-rbind(spPres, spBg)
presentBg<-c(rep(1, times=nrow(spPres)), rep(0, times=nrow(spBg)))

modIn<-cbind(trainData, presentBg)
write.csv(modIn, file=paste0(fold, "/", sp2[i], "_", bg.names[j], "_input.csv"), row.names=FALSE)

#since there are so few points, maybe doesn't make sense to leave aside test points at all
mod<-dismo::maxent(x=modIn[,3:(ncol(modIn)-1)], p=modIn[,ncol(modIn)], factors=c('soil_stype'), args=c('writebackgroundpredictions=true'), path='/home/bharti/D/PostDoc data/spatial_data/dataframes/maxent_temp')
mods[[j]]<-mod

#saving model results to disk
write.csv(mod@results, paste0(fold, "/", sp2[i], "_", bg.names[j], "_results.csv"))

#saving the model to disk
save(mod, file=paste0(fold, "/", sp2[i], "_", bg.names[j], "_model"))

#predicting the model to the extent of peninsular India
modMap<-predict(mod, preds)
modMaps<-stack(modMaps, modMap)

#saving the prediction raster to disk
writeRaster(modMap, filename=paste0(fold, "/", sp2[i], "_", bg.names[j], "_map"), format="GTiff", overwrite=TRUE)

#plot variable percent contribution 
jpeg(file=paste0(fold, "/", sp2[i], "_", bg.names[j], "_varCon.jpeg"), height=10, width = 10, units = 'in', res=300)
#plot(mod)
per.con<-var.importance(mod)[,c(1,2)]
per.con<-per.con[order(per.con$percent.contribution, decreasing=TRUE),]
per.con[,1]<-factor(per.con[,1], levels=per.con[order(per.con$percent.contribution, decreasing=TRUE), 1])
p<-ggplot(per.con, aes(x=variable, y=percent.contribution)) + geom_point() + ylim(0,100) + theme(axis.text.x=element_text(angle=90, hjust=0))
print(p)
dev.off()

#plot variable permutation importance
jpeg(file=paste0(fold, "/", sp2[i], "_", bg.names[j], "_perImp.jpeg"), height=10, width = 10, units = 'in', res=300)
per.imp<-var.importance(mod)[,c(1,3)]
per.imp<-per.imp[order(per.imp$permutation.importance, decreasing=TRUE),]
per.imp[,1]<-factor(per.imp[,1], levels=per.imp[order(per.imp$permutation.importance, decreasing=TRUE), 1])
p<-ggplot(per.imp, aes(x=variable, y=permutation.importance)) + geom_point() + ylim(0,100) + theme(axis.text.x=element_text(angle=90, hjust=0))
print(p)
dev.off()

#saving the response curve
jpeg(file=paste0(fold, "/", sp2[i], "_", bg.names[j], "_response.jpeg"),  height=10, width=12, units='in', res=300)
modRes<-response(mod, range='p')
dev.off()

#plot predicted map
jpeg(file=paste0(fold, "/", sp2[i], "_", bg.names[j], "_map.jpeg"),  height=10, width=12, units='in', res=300)
plot(modMap, main=paste0(sp[i], "_", bg.names[j]))
points(spPres[,c(1,2)])
dev.off()

#Since I am not setting aside any test points, I am adding the custom JK code from my PhD analysis
env.input<-modIn[,3:(ncol(modIn)-1)]
env.names<-colnames(env.input)
pres.input<-modIn[,ncol(modIn)]
modGain_ex<-as.data.frame(matrix(nrow=ncol(env.input), ncol=3))
colnames(modGain_ex)<-c("var_op", "env_var", "mod_gain" )
modGain_in<-modGain_ex

for(k in 1:ncol(env.input)){
exModel<-maxent(x=env.input[,-k], p=pres.input, factors=c('soil_stype'), args=c('writebackgroundpredictions=true'), path='/home/bharti/D/PostDoc data/spatial_data/dataframes/maxent_temp')

env.in<-as.data.frame(env.input[,k])
colnames(env.in)<-colnames(env.input)[k]

inModel<-maxent(x=env.in, p=pres.input, factors=c('soil_stype'), args=c('writebackgroundpredictions=true'), path='/home/bharti/D/PostDoc data/spatial_data/dataframes/maxent_temp')

#important to remember here that I am using REGULARIZED TRAINING GAIN 
modGain_ex[k,]<-c("without", env.names[k], exModel@results[2]) 
modGain_in[k,]<-c("with.only", env.names[k], inModel@results[2]) 
}

fullModel<-maxent(x=env.input, p=pres.input, factors=c('soil_stype'), args=c('writebackgroundpredictions=true'), path='/home/bharti/D/PostDoc data/spatial_data/dataframes/maxent_temp')
modGain_all<-c("full", "all variables", fullModel@results[2])
modGain<-rbind(modGain_in, modGain_ex, modGain_all)
modGain$mod_gain<-as.numeric(modGain$mod_gain)
modGain$var_op<-factor(modGain$var_op)
modGain$env_var<-factor(modGain$env_var)
modGain$var_op<-factor(modGain$var_op, levels(modGain$var_op)[c(2,3,1)])

jpeg(file = paste0(fold, "/", sp2[i], "_", bg.names[j], "_JK.jpeg"),  height=10, width=10, units='in', res=300)
#area.colour<-c(rep("1", times=ncol(env.input)), "2")
area.colour<-c("#31a354", "#e34a33", "#2b8cbe")
p<-ggplot(modGain, aes(env_var, mod_gain)) + geom_bar(aes(fill=var_op), position="dodge", stat="identity", width=0.75) + xlab("Predictor variable") + ylab("Regularised gain") + scale_fill_manual(values=area.colour, guide=guide_legend(reverse=TRUE)) + theme(legend.position="top") + theme(strip.text=element_text(size=14)) + theme(legend.title=element_blank()) + coord_flip()
print(p)
dev.off()

#evaluate model using the function evaluate
modEval<-evaluate(p=modIn[modIn$presentBg==1,3:(ncol(modIn)-1)], a=modIn[modIn$presentBg==0,3:(ncol(modIn)-1)], mod)

#saving the evaluation metrics to disk
save(modEval, file=paste0(fold, "/", sp2[i], "_", bg.names[j], "_eval"))

#saving AUC to disk
jpeg(file=paste0(fold, "/", sp2[i], "_", bg.names[j], "_ROC.jpeg"),  height=10, width=10, units='in', res=300)
plot(modEval, 'ROC')
dev.off()

#calculating TSS
mod.th<-threshold(modEval, stat='spec_sens')
mod.tpr<-modEval@TPR[modEval@t==mod.th]
mod.tnr<-modEval@TNR[modEval@t==mod.th]
mod.tss<-(mod.tpr+mod.tnr)-1
write.csv(mod.tss, file=paste0(fold, "/", sp2[i], "_", bg.names[j], "_tss.csv"), row.names=FALSE)

#saving probability density of presence and absence points to disk
jpeg(file=paste0(fold, "/", sp2[i], "_", bg.names[j], "_density.jpeg"),  height=10, width =20, units='in', res=300)
par(mfrow=c(1,2))
boxplot(modEval, notch=FALSE)
density(modEval)
dev.off()

#evaluating the models - cross-validation
#running the k-fold partition with my data-set and saving the results
if(nrow(spPres)>5){
modKfold<-maxent_kfold(5, pres=modIn[modIn$presentBg==1,3:(ncol(modIn)-1)], bg=modIn[modIn$presentBg==0,3:(ncol(modIn)-1)])
write.csv(modKfold, file=paste0(fold, "/", sp2[i], "_", bg.names[j], "_kfold.csv"), row.names=FALSE)
}

#carrying out the jackknifing test for small sample sizes
#creating a vector resJK which stores the results of the test points being predicted as present (1) or not (0)
resJK<-rep(NA, times=nrow(spPres))

#create a vector that saves the proportion presence in the map
propPres<-rep(NA, times=nrow(spPres))

#running evaluation as per Pearson, 2007 - running SDM with rare species
#run models dropping each presence location at a time
#create all possible arrangements of certain 0s and certain 1s
comb<-function(n,m){
t(combn(n,m,function(x)replace(numeric(n),x,1)))
}

if(nrow(spPres)<30){
for(k in 1:nrow(spPres)){
#creating input data - removing one location at a time (the jackknifing procedure)
spTrain<-spPres[-k,]
spTest<-as.data.frame(t(spPres[k,]))

#converting spTest into a SpatialPoints dataframe
coordinates(spTest)<-~Longitude+Latitude
crs(spTest)<-wgs

trainData<-rbind(spTrain, spBg)
presentBg<-c(rep(1, times=nrow(spTrain)), rep(0, times=nrow(spBg)))
modIn<-cbind(trainData, presentBg)

#running the jackknife model
modJK<-maxent(x=modIn[,3:(ncol(modIn)-1)], p=modIn[,ncol(modIn)], factors=c('soil_stype'), args=c('writebackgroundpredictions=true'), path='/home/bharti/D/PostDoc data/spatial_data/dataframes/maxent_temp')

#predicting the model to the extent of peninsular India
modJKMap<-dismo::predict(modJK, preds)

#finding the threshold which maximizes sum of sensitivity and specificity
thJK<-modJK@results[grep("*training.sensitivity.plus.specificity.Cloglog*", rownames(modJK@results))]

#creating a binary prediction map, where values greater than thJK are given the value of 1 and less than thJK are given a value of 0
modJKMap[modJKMap<thJK]<-0
modJKMap[modJKMap>thJK | modJKMap==thJK]<-1

#extract the value of left out present point from the binary map
resJK[k]<-extract(modJKMap, spTest)

#find the proportion of cells in the extent that have species present in the binary map
propPres[k]<-freq(x=modJKMap, value=1)/ncell(modJKMap)
}

#observed statistic d
d.obs<-sum(resJK*(1-propPres))

#calculate all possible permutations of resJK
nones<-length(which(resJK==1))
perm<-comb(length(resJK), nones)

#removing the actual result from the randomized result list
perm.str<-apply(apply(perm, c(1,2), as.character), 1, paste, collapse="")
perm<-perm[-which(perm.str==paste(as.character(resJK), collapse = "")),]

#add a condition when the model is never able to predict the heldout point correctly
if(length(perm)==0){
d.pval<-paste0("The JK held out predictions showed no variability. The output was ", paste(comb(length(resJK), nones), collapse=""))

write.csv(d.pval, file=paste0(fold, "/", sp2[i], "_", bg.names[j], "_smallSampleTest.csv"), row.names=FALSE)

next
}

#calculate the sum of product of perm with (1-propPres)
D.exp<-perm %*% (1-propPres)

#find which of these elements in D.exp is greater than d.obs
perm2<-perm[which(D.exp==d.obs | D.exp>d.obs),]
if(is.vector(perm2)){
perm2<-t(as.matrix(perm2))
} 
perm2.p<-perm2

for(m in 1:nrow(perm2)){
for(n in 1:ncol(perm2)){
if(perm[m,n]==0){
perm2.p[m,n]<-1-propPres[n]
}else{
perm2.p[m,n]<-propPres[n]
}
}  
}

#multiply the probabilities in each row and sum these values across rows
d.pval<-paste0("Probability that D under null distribution is greater than or equal to observed d: ", d.obs,  " is ", sum(apply(perm2.p, 1, prod)))

write.csv(d.pval, file=paste0(fold, "/", sp2[i], "_", bg.names[j], "_smallSampleTest.csv"), row.names=FALSE)
}
}

##comparing the models with different backgrounds
#calculating the number of parameters of each model
nparams<-unlist(lapply(mods, get.params))

#running the calc.aicc function for all the models - I need to check the order in which I need to mention lat and long
#Returns all NAs if the number of parameters is larger than the number of observations (occurrence localities).
aicc<-calc.aicc(nparams, spPres[,c(1,2)], modMaps)
aicc$models<-bg.names
aicc<-aicc[order(aicc$AICc), c(5,4,1,2,3)]
write.csv(aicc, file=paste0(fold, "/", sp2[i], "_aicc.csv"), row.names=FALSE)
}
###



###LOOSE CODE###
#replace level numbers for categorical variables with the actual name
#soil.stype.lev<-levels(soil.stype)[[1]]
#colnames(soil.stype.lev)<-c("soil_stype", "soil_stype_cat")
#pred.ext<-left_join(pred.ext, soil.stype.lev, by="soil_stype")

# According to Phillips (2006), the percent contribution of each variable is calculated as follows:
#   
#"While the Maxent model is being trained, it keeps track of which environmental variables are contributing to fitting the model. Each step of the Maxent algorithm increases the gain of the model by modifying the coefficient for a single feature; the program assigns the increase in the gain to the environmental variable(s) that the feature depends on. Converting to percentages at the end of the training process, we get the percent contribution."
# 
# "The percent contribution values are only heuristically defined: they depend on the particular path that the Maxent code uses to get to the optimal solution, and a different algorithm could get to the same solution via a different path, resulting in different percent contribution values. In addition, when there are highly correlated environmental variables, the percent contributions should be interpreted with caution."
# 
# Also according to Phillips (2006), the permutation importance of each variable is calculated as follows:
#   
#   "...for each environmental variable in turn, the values of that variable on training presence and background data are randomly permuted. The model is reevaluated on the permuted data, and the resulting drop in training AUC is shown in the table, normalized to percentages."
# 
# "The permutation importance measure depends only on the final Maxent model, not the path used to obtain it. The contribution for each variable is determined by randomly permuting the values of that variable among the training points (both presence and background) and measuring the resulting decrease in training AUC. A large decrease indicates that the model depends heavily on that variable. Values are normalized to give percentages."

##on Maxent output formats - from the Phillips Maxent Tutorial##

#Maxent supports four output formats for model values: raw, cumulative,logistic and cloglog. First, the raw output is just the Maxent exponential model itself.  Second, the cumulative value corresponding to a raw value of r is the percentage of the Maxent distribution with raw value atmost r.Cumulative output is best interpreted in terms of predicted omission rate: if we set a cumulative threshold of c, the resulting binary prediction would have omission rate c% on samples drawn from the Maxent distribution itself, and we can predicta similar omission rate for samples drawn from the species distribution.  Third, if c is the exponential of the entropy of the maxent distribution, then the logistic value corresponding to a raw value of r is c·r/(1+c·r).  This is a logistic function, because the raw value is an exponential function of the environmental variables.  The cloglog value corresponding to a raw value of r is 1-exp(-c·r).  The four output formats are all monotonically related, but they are scaled differently, and have different interpretations.
#The default output is cloglog, which is the easiest to conceptualize: it gives an estimate between 0 and 1 of probability of presence.  Note that probability of presence depends strongly on details of the sampling design, such as the quadratsize and (for vagile organisms) observation time; cloglog output estimates probability of presence assuming that the sampling design is such that typical presence localities have an expected abundance of one individual per quadrat, which results in a probability of presence of about 0.63.