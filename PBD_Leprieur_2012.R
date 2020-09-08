#######################################################################################################################################
#
# beta.pd.decompo
#
# R function to quantify phylogenetic beta diversity (PBD) and its 'true' turnover and phylogenetic diversity (PD) components based on the UniFrac and PhyloSor indices. 
#
# contacts: albouycamille@gmail.com, julien.debortoli@ird.fr and fabien.leprieur@univ-montp2.fr
#
#
# Arguments:
#   - com: presence/absence community matrix (Sites * Species)
#   - tree: rooted phylogenetic tree in phylo format 
#   - type: index used for decomposition: 
#       - "UniFrac" 
#       - "Phylosor" 
#       - "both" (default) for UniFrac + Phylosor
#   - output.dist: logical, to display results in distance matrix format (default: F)
#   - random: number of randomization to compute null model. Species are randomized across the tips from the tree while holding species richness and compositional beta diversity constant.
#
# Value: returns a list containing
#
#   - "betadiv": matrix with observed values of PBD and its turnover and PD components. 
#                               * PhyloSor and/or UniFrac = observed PBD  											                                            #								* PhyloSor_turn and/or UniFrac_turn = observed phylogenetic turnover component of PBD
#                               * PhyloSor_PD and/or UniFrac_PD = observed PD component of PBD
#                If randomization=T, standardized effect size (SES) values are also given for each index.
#                               * SES_PhyloSor and/or SES_UniFrac = Standardized Effect Size of PBD
#                               * SES_PhyloSor_turn and/or SES_UniFrac_turn = Standardized Effect Size of the turnover component of PBD
#                               * SES_PhyloSor_PD and/or SES_UniFrac_PD = Standardized Effect Size of the PD component of PBD
#
#                If output.dist=T, results are given as a list of distance matrices.
#               
#                               
#
#   - "util.pd": matrix with paired values of phylogenetic diversity (i.e. PD values for each community and PD values for the two communities combined)
#
#######################################################################################################################################

library(picante)

############ Products from a phylogenetic tree to compute Phylogenetic Beta Diversity decomposition (utility function) #######################

beta.pd.utils <- function(com, tree) {
	
  combin  <-  combn(nrow(com),2)
  labcomb <-  apply(combin,2,function(x) paste(rownames(com)[x],collapse="-"))
  pd.obs  <-  pd2(com,tree)[,"PD"]
  com.tot <- t(apply(combin,2,function(x) colSums(com[x,])>0))
  pd.obs.tot <- pd2(com.tot,tree)[,"PD"]
  sum.pd.obs <- apply(combin,2,function(x) sum(pd.obs[x]))
  min.pd.obs <- apply(pd.obs.tot-t(combn(pd.obs,2)),1,min)
  dif.pd.obs <- apply(combin,2,function(x) diff(pd.obs[x]))
     	return(list(pd.obs=pd.obs,pd.obs.tot=pd.obs.tot,sum.pd.obs=sum.pd.obs,min.pd.obs=min.pd.obs,dif.pd.obs=dif.pd.obs,labcomb=labcomb,combin=combin))

  } # end of function beta.pd.utils

############ Phylogenetic Beta Diversity decomposition from the PhyloSor index (utility function) #######################

betasor.pd <- function(betapd) {
	
  PhyloSor <- (2*betapd$pd.obs.tot-betapd$sum.pd.obs)/betapd$sum.pd.obs
  PhyloSor_turn <- betapd$min.pd.obs/(betapd$sum.pd.obs-betapd$pd.obs.tot+betapd$min.pd.obs)
  PhyloSor_PD <- (abs(betapd$dif.pd.obs)/betapd$sum.pd.obs)*((betapd$sum.pd.obs-betapd$pd.obs.tot)/(betapd$sum.pd.obs-betapd$pd.obs.tot+betapd$min.pd.obs))
  return(as.matrix(data.frame(PhyloSor,PhyloSor_turn,PhyloSor_PD,row.names=betapd$labcomb)))

} # end of function betasor.pd

############ Phylogenetic Beta Diversity decomposition from the UniFrac index (utility function) #######################

betajac.pd <- function(betapd) {
	
  UniFrac <- (2*betapd$pd.obs.tot-betapd$sum.pd.obs)/betapd$pd.obs.tot
  UniFrac_turn <- 2*(betapd$min.pd.obs)/(betapd$sum.pd.obs-betapd$pd.obs.tot+2*betapd$min.pd.obs)
  UniFrac_PD <- UniFrac-UniFrac_turn
  return(as.matrix(data.frame(UniFrac,UniFrac_turn,UniFrac_PD,row.names=betapd$labcomb)))
  
} # end offunction betajac.pd

############ General function for Phylogenetic Beta Diversity decomposition #######################

beta.pd.decompo <- function(com, tree,type="both",output.dist=F, random=F) {
  TYPE=c("both","PhyloSor","Unifrac")
  type=pmatch(type,TYPE)
  if(class(tree)!="phylo") stop("### invalid tree's format: \"phylo\" format required ###\n\n")
  if(is.na(type)) stop("### invalid index type's argument ###\n\n")
  if(any(!(colnames(com)%in%tree$tip))) cat("\n### warnings: some species in community matrix not included in the tree ###\n\n")
  
  if(type==2) {
    util.pd <- beta.pd.utils(com,tree)
    betadiv <- betasor.pd(util.pd)
    }
  
  if(type==3) {
    util.pd <- beta.pd.utils(com,tree)
    betadiv <- betajac.pd(util.pd)
    }
  
  if(type==1) {
    util.pd <- beta.pd.utils(com,tree)
    betadiv <- cbind(betasor.pd (util.pd),betajac.pd(util.pd))
    }
 
  if(random) {
    el <- ifelse(type==1,6,3)
    sim <- array(dim=c(dim(combn(nrow(com),2))[2],el,random))
    i=1
    x11(h=2)
    plot(seq(1:random),rep(1,random),type="n",axes=F,xlab="Iteration",ylab="")
    axis(1);box()
    
    while(i<=random)  {
      rdtree <- tree
      rdtree$tip.label <- sample(rdtree$tip.label)
      if(type==2) sim[,,i] <- betasor.pd (beta.pd.utils(com,rdtree))
      
      if(type==3) sim[,,i] <- betajac.pd(beta.pd.utils(com,rdtree))
     
      if(type==1) {
        util.pd.sim <- beta.pd.utils(com,rdtree)
        sim[,,i] <- cbind(betasor.pd (util.pd.sim),betajac.pd(util.pd.sim))
      }
      
      points(i,1,pch=15)
      i=i+1
      
    }
  dev.off()
  
  rand.beta.stat <- apply(sim,c(1,2),function(x) c(mean(x,na.rm=T),sd(x,na.rm=T)))
  betadivpd.z <- (betadiv-rand.beta.stat[1,,])/rand.beta.stat[2,,]
  betadiv <- cbind(betadiv,betadivpd.z)
  
  if(type==1){ 
colnames(betadiv) <- c("PhyloSor","PhyloSor_turn","PhyloSor_PD","UniFrac","UniFrac_turn",
"UniFrac_PD","SES_PhyloSor","SES_PhyloSor_turn","SES_PhyloSor_PD","SES_UniFrac","SES_UniFrac_turn","SES_UniFrac_PD")
  }
  else {
    if(type==2) colnames(betadiv) <- c("PhyloSor","PhyloSor_turn","PhyloSor_PD",
      "SES_PhyloSor","SES_PhyloSor_turn","SES_PhyloSor_PD")
    else colnames(betadiv)=c("UniFrac","UniFrac_turn","UniFrac_PD","SES_UniFrac","SES_UniFrac_turn","SES_UniFrac_PD") 
  }
  if(output.dist) {
    results=list()
    for (i in 1:ncol(betadiv)) {
      results[[i]] <- dist.mat(com,betadiv[,i])
      names(results)[[i]] <- colnames(betadiv)[i]
      }
    return(list(betadiv=results,util.pd=cbind(matrix(util.pd$pd.obs[util.pd$combin],ncol=2,byrow=T,dimnames=list(util.pd$labcomb,c("Com1","Com2"))),pd.obs.tot=util.pd$pd.obs.tot)))
    }
  else return(list(betadiv=betadiv,util.pd=cbind(matrix(util.pd$pd.obs[util.pd$combin],ncol=2,byrow=T,dimnames=list(util.pd$labcomb,c("Com1","Com2"))),pd.obs.tot=util.pd$pd.obs.tot)))
  }
  if(output.dist) {
    results=list() 
    for (i in 1:ncol(betadiv)) {
      results[[i]] <- dist.mat(com,betadiv[,i])
      names(results)[[i]] <- colnames(betadiv)[i]
      }
    return(list(betadiv=results,util.pd=cbind(matrix(util.pd$pd.obs[util.pd$combin],ncol=2,byrow=T,dimnames=list(util.pd$labcomb,c("Com1","Com2"))),pd.obs.tot=util.pd$pd.obs.tot)))
    }
  else return(list(betadiv=betadiv,util.pd=cbind(matrix(util.pd$pd.obs[util.pd$combin],ncol=2,byrow=T,dimnames=list(util.pd$labcomb,c("Com1","Com2"))),pd.obs.tot=util.pd$pd.obs.tot)))

} # end of function beta.pd.decompo

############ Paired matrix to distance matrix conversion (utility function) #######################

dist.mat <- function(com,pair) {
	
  ncom <- nrow(com)
  distmat <- matrix(nrow=ncom,ncol=ncom,0,dimnames=list(rownames(com),rownames(com)))
  st <- c(0,cumsum(seq(ncom-1,2)))+1
  end <- cumsum(seq(ncom-1,1))
  for (i in 1:(ncom-1)) distmat[i,(ncom:(seq(1,ncom)[i]))]=c(pair[end[i]:st[i]],0)
  distmat <- as.dist(t(distmat))
  return(distmat)
  
} # end of function dist.mat

############ Phylogenetic Diversity Faith (adapted from pd function in picante library) #######################

pd2 <- function (samp, tree, include.root = TRUE) {
	
  if (is.null(tree$edge.length)) {
    stop("Tree has no branch lengths, cannot compute pd")
    }
  species <- colnames(samp)
  tree  <- node.age(tree)
  PDout <- apply(samp,1, function(x) {
    present <- species[x > 0]
    treeabsent <- tree$tip.label[which(!(tree$tip.label %in%present))]
    if (length(present) == 0) {
      PD <- 0
      }
    else if (length(present) == 1) {
      if (!is.rooted(tree) || !include.root) {
        warning("Rooted tree and include.root=TRUE argument required to calculate PD of single-species sampunities. Single species sampunity assigned PD value of NA.")
        PD <- NA
        }
      else {
        PD <- tree$ages[which(tree$edge[, 2] ==
        which(tree$tip.label == present))]
        }
      }  
    else if (length(treeabsent) == 0) {
      PD <- sum(tree$edge.length)
      }
    else {
      sub.tree <- drop.tip(tree, treeabsent)
      if (include.root) {
        if (!is.rooted(tree)) {
          stop("Rooted tree required to calculate PD with include.root=TRUE argument")
          }
        sub.tree <- node.age(sub.tree)
        sub.tree.depth <- max(sub.tree$ages)
        orig.tree.depth <- max(tree$ages[which(tree$edge[,2] %in% which(tree$tip.label %in% present))])
        PD <- sum(sub.tree$edge.length) + (orig.tree.depth - sub.tree.depth)
        }
      else {
        PD <- sum(sub.tree$edge.length)
        }
      }   
    SR <- length(present)  
    PDout <- c(PD,SR)
    } )         
  PDout <- t(PDout)
  rownames(PDout) <- rownames(samp)
  colnames(PDout) <- c("PD","SR")   
  return(PDout)  
} # end of function pd2