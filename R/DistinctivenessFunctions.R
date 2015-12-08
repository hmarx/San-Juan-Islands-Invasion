###### Distinctiveness Functions #######
###### 6 Nov 2014 ######################
###### Hannah E. Marx ##################

### For analysis of phylogenetic and functional distinctivness 
###### to nearest native (DNNS/NNFD) 
###### and mean native community (MDNS/MFD)
### Also, calculate standardized effect size compared to null model randomizing invasive species occurences
### Used for Marx et al. 2015 Diveristy and Distributions

## Data (analysis.R)
## Phylogeny: phylo object with all species in "source pool"
## Community data: data.frame with species names matching phylogeny, presence/absence of speices in each community: 
######## row = species names; columns = 1 (presnce), or 0 (absence)
## Trait data: data.frame with species names mathcing phylogeny, trait values 
######## row = species names (must match format of names in phylogeny); columns = trait value


##################################### PHYLOGENETIC DISTINCTIVENESS ##################################### 
############################ Mean Nearest Native Phylogenetic Distance (DNNS)
############################ Mean Phylogenetic Distance to Native Community (MDNS)

## phy = phylogenetic tree of species pool
## communinty = occurence data for species (rownames match phy), in each community (colnames): invasive species = i, native spcecies = n, absence = 0 
## col = the community to calcualte metrics for

phyloDistinct  <- function(phy, community, col){
  dataPruned <- community[which(!community[col] == "0"), col,  drop=F]   ## prune community step
  data3 <- as.list(unlist(dataPruned))
  names(data3) <- rownames(dataPruned)
  dataPruned <- data3
  dist.data <- data.frame() # make an empty data frame
  if (length(dataPruned) == 1){  #if there is only one species, return NAs
    tmp <- cbind("NA", "NA", "NA","NA", "NA", "NA", "NA")
    dist.data <- rbind(dist.data, tmp)
    colnames(dist.data) <- c("Species.Status", "Nearest.native", "MinDist.Nearest.native", "MeanDist.NativeCommunity", "Nearest.invasive", "MinDist.Nearest.invasive", "MeanDist.InvasiveCommunity") # add column names
    rownames(dist.data) <- names(dataPruned) 
    return(dist.data)
  } 
  n.sp <- length(dataPruned[dataPruned == "n"]) # get all the species present in the community which are native
  i.sp <- length(dataPruned[dataPruned == "i"]) # get all the species present in the community which are invasive
  print(c(col, "Total native species =", n.sp), quote=F) ## print the number of native species on the island
  print(c(col, "Total invasive species =", i.sp), quote=F) ## print the number of native species on the island
  tmp <- treedata(phy, dataPruned) ## check the names 
  new.phy <- tmp$phy ## prune the tree to include just species in the community of interest 
  com.data <- tmp$data   ## subset the community data, to make sure in the right order
  n.tips <- rownames(com.data)[com.data == "n"] # get all the tips which are native
  print(c(col, "Total native tips =", length(n.tips)), quote=F) ## print the number of natives after pruning to tree
  i.tips <- rownames(com.data)[com.data == "i"] # get all the tips which are invasive
  print(c(col, "Total invasive tips =", length(i.tips)), quote=F) #### print the number of invasive after pruning to tree
  dist <- cophenetic.phylo(new.phy) # calculate cophenetic distance for community 
  for (i in 1:ncol(dist)){
    i.tips.new <- i.tips[i.tips != colnames(dist)[i]] # remove self comparison for invasives
    n.tips.new <- n.tips[n.tips != colnames(dist)[i]] # remove self comparison for natives
    ifelse(colnames(dist)[i] %in% i.tips, sps <- "i", sps <- "n") #lable taxa as Native/introduced
    
    tmp.i <- min(dist[which(rownames(dist) %in% i.tips.new),i]) # distance to nearest invasive species
    tmp.n <- min(dist[which(rownames(dist) %in% n.tips.new),i]) # distance to nearest native species 
    
    near.n <- n.tips.new[which(dist[n.tips.new,i] == tmp.n)] # name of nearest native
    near.i <- i.tips.new[which(dist[i.tips.new,i] == tmp.i)] # name of nearest invasive
    
    tmp.i.mean <- mean(dist[which(rownames(dist) %in% i.tips.new),i]) # mean distance to invasive community
    tmp.n.mean <- mean(dist[which(rownames(dist) %in% n.tips.new),i]) # mean distance to native community
    
    if (length(near.n) > 1){
      nn <- data.frame()
      for (k in 1:length(near.n)){
        t <- dataPruned[names(dataPruned)==near.n[k]]
        nn <- paste(nn, near.n[k], sep=".")
      }
      
    } else {
      if (length(near.n) == 0){
        nn <- "NA"
      } else{
        nn <- names(dataPruned[names(dataPruned)==near.n])
      }
    }
    
    if (length(near.i) > 1){
      ii <- data.frame()
      for (k in 1:length(near.i)){
        t <- dataPruned[names(dataPruned)==near.i[k]]
        ii <- paste(ii, near.i[k], sep=".")
      }
      
    } else {
      if (length(near.i) == 0){
        ii <- "NA"
      } else {
        ii <- names(dataPruned[names(dataPruned)==near.i])
      }
    }
    
    if (length(near.i)==0){
      tmp.i <- 0
    }
    if (length(near.n)==0){
      tmp.n <- 0
    }
    
    tmp <- cbind(sps, nn, tmp.n, tmp.n.mean, ii, tmp.i, tmp.i.mean)
    dist.data <- rbind(dist.data, tmp)
  } 
  colnames(dist.data) <- c("Species.Status", "Nearest.native", "MinDist.Nearest.native", "MeanDist.NativeCommunity", "Nearest.invasive", "MinDist.Nearest.invasive", "MeanDist.InvasiveCommunity") # add column names
  rownames(dist.data) <- colnames(dist) # add rownames corresponding to the tip of interest
  return(dist.data)
  
}


############################ Summarize phyloDistinct output
# output = output from phyloDistinct()
summary.DNNS.MDNS<- function(output){
  summ <- data.frame()
  if (nrow(output)==1){
    summ <- rbind("NA", "NA", "NA", "NA", "NA","NA", "NA", "NA", "NA","NA", "NA", "NA")
    rownames(summ) <- c("meanDNNSnatives", "meanDNNSinvasives", "n.natives", "n.invasives", "t.DNNS.p.value", "t.DNNS.conf.int.Lo","t.DNNS.conf.int.Hi", 
                        "meanMDNSnatives", "meanMDNSinvasives", "t.MDNS.p.value", "t.MDNS.conf.int.Lo","t.MDNS.conf.int.Hi")
    return(summ)
  }
  DNNSn <- output[output["Species.Status"]=="n",] # get all native
  DNNSi <- output[output["Species.Status"]=="i",] # get all native
  
  if (nrow(DNNSn) > 1 && nrow(DNNSi) > 1) {
    meanDNNSn <- mean(as.numeric(as.character(DNNSn[,"MinDist.Nearest.native"])))  
    meanDNNSi <- mean(as.numeric(as.character(DNNSi[,"MinDist.Nearest.native"]))) 
    MDNSnn <- mean(as.numeric(as.character(DNNSn[,"MeanDist.NativeCommunity"])))
    MDNSin <- mean(as.numeric(as.character(DNNSi[,"MeanDist.NativeCommunity"])))
    
    t.DNNS <- t.test(as.numeric(as.character(DNNSn[,"MinDist.Nearest.native"])), as.numeric(as.character(DNNSi[,"MinDist.Nearest.native"])), paired=F)
    t.MDNS <- t.test(as.numeric(as.character(DNNSn[,"MeanDist.NativeCommunity"])), as.numeric(as.character(DNNSi[,"MeanDist.NativeCommunity"])), paired=F)
    
    summ <- rbind(meanDNNSn, meanDNNSi, nrow(DNNSn), nrow(DNNSi), t.DNNS$p.value, t.DNNS$conf.int[1], t.DNNS$conf.int[2], MDNSnn, MDNSin, t.MDNS$p.value, t.MDNS$conf.int[1], t.MDNS$conf.int[2])
    
  } else if (nrow(DNNSn) > 1 && nrow(DNNSi) == 1) {
    meanDNNSn <- mean(as.numeric(as.character(DNNSn[,"MinDist.Nearest.native"])))  
    meanDNNSi <- mean(as.numeric(as.character(DNNSi[,"MinDist.Nearest.native"]))) 
    MDNSnn <- mean(as.numeric(as.character(DNNSn[,"MeanDist.NativeCommunity"])))
    MDNSin <- mean(as.numeric(as.character(DNNSi[,"MeanDist.NativeCommunity"])))
    
    t.DNNS <- "NA"
    t.MDNS <- "NA"
    
    summ <- rbind(meanDNNSn, meanDNNSi, nrow(DNNSn), nrow(DNNSi), "NA", "NA", "NA", MDNSnn, MDNSin, "NA", "NA", "NA")
    
  } else {
    summ <- rbind("NA", "NA", "NA", "NA", "NA","NA", "NA", "NA", "NA","NA", "NA", "NA")
    
  }
  
  rownames(summ) <- c("meanDNNSnatives", "meanDNNSinvasives", "n.natives", "n.invasives", "t.DNNS.p.value", "t.DNNS.conf.int.Lo","t.DNNS.conf.int.Hi", 
                      "meanMDNSnatives", "meanMDNSinvasives", "t.MDNS.p.value", "t.MDNS.conf.int.Lo","t.MDNS.conf.int.Hi")
  
  
  return(summ)
}

sig.obs.phyloDiversty <- function(summaryDF, metadataFULL, plottype=c("summary.Bar", "IslSizeDNNS", "IslSizeDNNS", "4islDNNS", "4islMDNS")){
  summaryDF <- merge(summaryDF, metadata, by=0) # add metadata for islands for plotting and discussing
  nisls <- nrow(summaryDF)
  # indicate islands with significant difference in observed means between status groups for each island
  summaryDF[,"Significance"] <- ifelse(summaryDF[,"t.DNNS.p.value"] <= 0.05, 1,
                                         ifelse(summaryDF[,"t.MDNS.p.value"] <= 0.05, 2, 0)) 
  
  #### summarize significance data by size categories
  summaryDF$Size.cat <- as.character(summaryDF$Size.cat)
  summaryDF$Size.cat <- factor(summaryDF$Size.cat, levels=c("sm", "med", "lg"))
  
  sigDNNS = cbind("Island"=summaryDF[,1], "Size.cat"=as.character(summaryDF[,"Size.cat"]), "Metric"=rep(x = "DNNS", times = nrow(summaryDF)), "Significance"=ifelse(summaryDF[,"t.DNNS.p.value"] <= 0.05, 1,0))
  sigMDNS = cbind("Island"=summaryDF[,1], "Size.cat"=as.character(summaryDF[,"Size.cat"]), "Metric"=rep(x = "MDNS", times = nrow(summaryDF)), "Significance"=ifelse(summaryDF[,"t.MDNS.p.value"] <= 0.05, 1,0))
  obs.phylo <- as.data.frame(rbind(sigDNNS, sigMDNS))
  
  ##### significance in observed difference bewteen mean DNNSi / mean DNNSn 
  nsdifDNNS <- length(summaryDF$t.DNNS.p.value[which(summaryDF$t.DNNS.p.value >= 0.05)]) # 34 islands NS
  sigdifDNNS <- length(summaryDF$t.DNNS.p.value[which(summaryDF$t.DNNS.p.value <= 0.05)]) # 40 isalnds with significant differnece in mean betweenmean MMNPDi / mean MMNPDn 
  summaryDF.DNNS.sig <- summaryDF[which(summaryDF$t.DNNS.p.value <= 0.05),]
  #summaryDF[which(summaryDF$t.DNNS.p.value <= 0.005),] # 3
  per.diff.DNNS <- sigdifDNNS/nisls # 0.5405405 
  nat.greater.DNNS <- sum(ifelse(as.numeric(summaryDF$meanDNNSnatives) > as.numeric(summaryDF$meanDNNSinvasives), 1,0)) #74
  nat.less.DNNS <- sum(ifelse(as.numeric(summaryDF$meanDNNSnatives) < as.numeric(summaryDF$meanDNNSinvasives), 1,0)) #74
  
  sm.sigDNNS <- nrow(merge(metadataFULL[which(metadataFULL$Size.cat == "sm"),], summaryDF.DNNS.sig, by=1)) # 7/74 small islands sig DNNS 0.09459459
  med.sigDNNS <- nrow(merge(metadataFULL[which(metadataFULL$Size.cat == "med"),], summaryDF.DNNS.sig, by=1)) # 20/74 medium isl sig DNNS 0.2702703
  lg.sigDNNS <- nrow(merge(metadataFULL[which(metadataFULL$Size.cat == "lg"),], summaryDF.DNNS.sig, by=1)) # 13/74 large isl sig DNNS 0.1756757
  
  #### significance in observed difference bewteen mean MDNSin / mean MDNSnn
  nsdifMDNS <- length(summaryDF$t.MDNS.p.value[which(summaryDF$t.MDNS.p.value >= 0.05)]) # 53 NS
  sigdifMDNS <- length(summaryDF$t.MDNS.p.value[which(summaryDF$t.MDNS.p.value <= 0.05)]) # 21 significant differneces between mean MDNS
  summaryDF.MDNS.sig <- summaryDF[which(summaryDF$t.MDNS.p.value <= 0.05),]
  #summaryDF[which(summaryDF$t.MDNS.p.value <= 0.005),] # 5
  per.diff.MDNS <- sigdifMDNS/nisls # 0.2837838
  nat.greater.MDNS <-sum(ifelse(as.numeric(summaryDF$meanMDNSnatives) > as.numeric(summaryDF$meanMDNSinvasives), 1,0)) #70
  nat.less.MDNS <-sum(ifelse(as.numeric(summaryDF$meanMDNSnatives) < as.numeric(summaryDF$meanMDNSinvasives), 1,0)) #4
  
  sm.sigMDNS <- nrow(merge(metadataFULL[which(metadataFULL$Size.cat == "sm"),], summaryDF.MDNS.sig, by=1)) # 2/74 small islands sig MDNS 0.02702703
  med.sigMDNS <- nrow(merge(metadataFULL[which(metadataFULL$Size.cat == "med"),], summaryDF.MDNS.sig, by=1)) # 4/74 med islands sig MDNS 0.05405405
  lg.sigMDNS <- nrow(merge(metadataFULL[which(metadataFULL$Size.cat == "lg"),], summaryDF.MDNS.sig, by=1)) # 15/74 large isl sig MDNS 0.2027027
  
  dnns <- c(nisls, sigdifDNNS, nsdifDNNS, per.diff.DNNS, nat.greater.DNNS, nat.less.DNNS, sm.sigDNNS, med.sigDNNS, lg.sigDNNS)
  mdns <- c(nisls, sigdifMDNS, nsdifMDNS, per.diff.MDNS, nat.greater.MDNS, nat.less.MDNS, sm.sigMDNS, med.sigMDNS, lg.sigMDNS)
  
  sigobs <- rbind(dnns, mdns)
  colnames(sigobs) <- c("nIslands", "sigDiff", "ns", "%Sig", "mean_N>I", "mean_N<I", "sigSmall", "sigMed", "sigLg")

  
  if (plottype[1] == "summary.Bar"){
    ### Bar plot of significance 
    p <- ggplot(obs.phylo, aes(x=Significance, fill=factor(Size.cat))) + 
      geom_bar(stat="bin") +
      scale_fill_manual(name="Size Category",
                        breaks=c("sm", "med", "lg"),
                        values=c("sm"="dodgerblue4", "med"="orangered3", "lg"="gold1"),
                        labels=c("sm"="small", "med"="medium", "lg"="large")) +
      theme_bw() +
      scale_x_discrete(name="", labels=c("0"="NS", "1"="Significant")) +
      facet_wrap(~ Metric) +
      ggtitle("Observed difference of \n Phylogenetic Distances\n") 
    return(list(p,sigobs))
  }
  
  
  ########################  Plot observed phylogenetic distinctiveness for each island, ordered by increasing island size
  ### append metadata for island area to column in list of DNNS/MDNS for each island
  list <- phyloObs
  list.meta <- list()
  for (i in 1:length(list)){ 
    tmp <- metadata[names(list[i]), "Area.m2"]
    rep(x=tmp, times=nrow(list[[i]]))
    newlist <- mapply(cbind, list[i], "Area"=tmp, SIMPLIFY=F) 
    list.meta[[i]] <- newlist
  }
  head(list.meta)
  
  ### melt into datafrome for plotting
  phyloObs_melt  <- suppressMessages(melt(data=list.meta))
  #head(phyloObs_melt)
  #dim(phyloObs_melt) # 4835   11
  if (plottype[1] == "none"){
    return(list(sigobs,phyloObs_melt))
    
  }
  
  if (plottype[1] =="IslSizeDNNS"){
    ## Regression of differnence in means to island size
    model <- lm(as.numeric(as.character(MinDist.Nearest.native)) ~ value , data=phyloObs_melt)
    summary(model)
    anova(model)
    r2.DNNS <- paste("R^2 = ", signif(summary(model)$r.squared, 3), sep="") #Explained variation / Total variation
    p.DNNS <- paste("p-value = ", signif(anova(model)[[5]][1], 3), "***", sep="")
    
    head(phyloObs_melt)
    DNNS <- ggplot(phyloObs_melt, aes(x=reorder(factor(L2),value), y=as.numeric(as.character(MinDist.Nearest.native))), position=position_dodge(width=1))#, col=c("magenta1", "green3"))
    DNNS <- DNNS + geom_boxplot(aes(fill = Species.Status), width = 1)
    DNNS <- DNNS + geom_smooth(method="lm", se=T, color="black", aes(group=1))
    DNNS <- DNNS + coord_cartesian(ylim=c(.5, 5000)) 
    DNNS <- DNNS + scale_x_discrete("Island (increasing size)") #, breaks=seq(0, 80, 10)) 
    DNNS <- DNNS + scale_y_log10("Log DNNS") #+ coord_fixed(ratio=4) 
    DNNS <- DNNS + theme_bw() 
    DNNS <- DNNS + scale_fill_manual(values=c("i"= "magenta1", "n"="green3"), labels=c("i"= "MMNPDi", "n"="MMNPDn")) ##breaks=rev(factor(SJnew$status)),
    DNNS <- DNNS + guides(fill=guide_legend(title=""))
    DNNS <- DNNS + theme(legend.position="top")
    DNNS <- DNNS + theme(axis.text.x = element_text(angle = -45, hjust = 0))
    DNNS <- DNNS + ggtitle("Phylogenetic Distance to the Nearest Native Species (DNNS)") + theme(plot.title=element_text(size=rel(1.5)))
    DNNS <- DNNS + annotate("text", label=r2.DNNS, x=3, y=0.9, size=4) #y=max(as.numeric(as.character(SJ_NN_meltNEW$MinDist.Nearest.native))-20)
    DNNS <- DNNS + annotate("text", label=p.DNNS, x=4, y=0.7, size=4) 
    DNNS <- DNNS + theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
    return(list(DNNS,sigobs))
    
  }

  if (plottype[1] =="IslSizeMDNS"){
    ## Regression of differnence in means to island size
    modelMDNS <- lm(as.numeric(as.character(MeanDist.NativeCommunity)) ~ value , data=phyloObs_melt)
    summary(modelMDNS)
    anova(modelMDNS)
    r2.MDNS <- paste("R^2 = ", signif(summary(modelMDNS)$r.squared, 3)) #Explained variation / Total variation
    p.MDNS <- paste("p-value = ", signif(anova(modelMDNS)[[5]][1], 3), "**", sep="")
    
    MDNS <- ggplot(phyloObs_melt, aes(x=reorder(factor(L2),value), y=as.numeric(as.character(MeanDist.NativeCommunity))), position=position_dodge(width=1))#, col=c("magenta1", "green3"))
    MDNS <- MDNS + geom_boxplot(aes(fill = Species.Status), width = 1)
    MDNS <- MDNS + geom_smooth(method="lm", se=T, color="black", aes(group=1))
    MDNS <- MDNS + coord_cartesian(ylim=c(100, 2000))
    MDNS <- MDNS + scale_x_discrete("Island (increasing size)") #, breaks=seq(0, 80, 10)) 
    MDNS <- MDNS + scale_y_log10("Log MDNS")  
    MDNS <- MDNS + theme_bw() 
    MDNS <- MDNS + scale_fill_manual(values=c("i"= "magenta1", "n"="green3"), labels=c("i"= "MDNSNi", "n"="MDNSNn")) ##breaks=rev(factor(SJnew$status)),
    MDNS <- MDNS + guides(fill=guide_legend(title=""))
    MDNS <- MDNS + theme(legend.position="top")
    MDNS <- MDNS + theme(axis.text.x = element_text(angle = -45, hjust = 0))
    MDNS <- MDNS + ggtitle("Mean Phylogenetic Distance to Native Community (MDNS)") + theme(plot.title=element_text(size=rel(1.5)))
    MDNS <- MDNS + annotate("text", label=r2.MDNS, x=3, y=130, size=4) #y=max(as.numeric(as.character(SJ_NN_meltNEW$MinDist.Nearest.native))-20)
    MDNS <- MDNS + annotate("text", label=p.MDNS, x=3.5, y=120, size=4) 
    MDNS <- MDNS + theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
    return(list(MDNS,sigobs))
    
  }
  
  if (plottype[1] =="4islDNNS"){
    four.islands.phyloObs_melt <- phyloObs_melt[phyloObs_melt$L2 %in% four.islands, ]
    DNNS.four <- ggplot(four.islands.phyloObs_melt, aes(x=reorder(factor(L2),value), y=as.numeric(as.character(MinDist.Nearest.native))), position=position_dodge(width=1))#, col=c("magenta1", "green3"))
    DNNS.four <- DNNS.four+ geom_boxplot(aes(fill = Species.Status), width = 1)
    #DNNS.four+ geom_smooth(method="lm", se=T, color="black", aes(group=1))
    DNNS.four <- DNNS.four+ coord_cartesian(ylim=c(.5, 5000))
    DNNS.four <- DNNS.four+ scale_x_discrete("Island (increasing size)") #, breaks=seq(0, 80, 10)) 
    DNNS.four <- DNNS.four+ scale_y_log10("Log DNNS") #+ coord_fixed(ratio=4) 
    DNNS.four <- DNNS.four+ theme_bw() 
    DNNS.four <- DNNS.four+ scale_fill_manual(values=c("i"= "magenta1", "n"="green3"), labels=c("i"= "MMNPDi", "n"="MMNPDn")) ##breaks=rev(factor(SJnew$status)),
    DNNS.four <- DNNS.four+ guides(fill=guide_legend(title=""))
    DNNS.four <- DNNS.four+ theme(legend.position="top")
    DNNS.four <- DNNS.four+ theme(axis.text.x = element_text(angle = -45, hjust = 0))
    DNNS.four <- DNNS.four+ ggtitle("Phylogenetic Distance to the Nearest Native Species (DNNS)") + theme(plot.title=element_text(size=rel(1.5)))
    #DNNS.four+ annotate("text", label=r2.DNNS, x=5, y=0.09, size=4) #y=max(as.numeric(as.character(SJ_NN_meltNEW$MinDist.Nearest.native))-20)
    #DNNS.four+ annotate("text", label=p.DNNS, x=5, y=0.05, size=4) 
    DNNS.four <- DNNS.four+ theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
    return(list(DNNS.four,sigobs))
  
  }
  
  if (plottype[1] =="4islMDNS"){
    four.islands.phyloObs_melt <- phyloObs_melt[phyloObs_melt$L2 %in% four.islands, ]
    MDNS.four <- ggplot(four.islands.phyloObs_melt, aes(x=reorder(factor(L2),value), y=as.numeric(as.character(MeanDist.NativeCommunity))), position=position_dodge(width=1))#, col=c("magenta1", "green3"))
    MDNS.four <- MDNS.four + geom_boxplot(aes(fill = Species.Status), width = 1)
    MDNS.four <- MDNS.four  + coord_cartesian(ylim=c(100, 2000))
    #MDNS.four <- MDNS.four  + geom_smooth(method="lm", se=T, color="black", aes(group=1))
    MDNS.four <- MDNS.four  + scale_x_discrete("Island (increasing size)") #, breaks=seq(0, 80, 10)) 
    MDNS.four <- MDNS.four  + scale_y_log10("Log MDNSN")  
    MDNS.four <- MDNS.four  + theme_bw() 
    MDNS.four <- MDNS.four  + scale_fill_manual(values=c("i"= "magenta1", "n"="green3"), labels=c("i"= "MDNSi", "n"="MDNSn")) ##breaks=rev(factor(SJnew$status)),
    MDNS.four <- MDNS.four  + guides(fill=guide_legend(title=""))
    MDNS.four <- MDNS.four  + theme(legend.position="top")
    MDNS.four <- MDNS.four  + theme(axis.text.x = element_text(angle = -45, hjust = 0))
    MDNS.four <- MDNS.four  + ggtitle("Mean Phylogenetic Distance to Native Community (MDNS)") + theme(plot.title=element_text(size=rel(1.5)))
    #MDNS.four <- MDNS.four  + annotate("text", label=r2.MDNS, x=5, y=100, size=4) #y=max(as.numeric(as.character(SJ_NN_meltNEW$MinDist.Nearest.native))-20)
    #MDNS.four <- MDNS.four  + annotate("text", label=p.MDNS, x=5, y=90, size=4) 
    MDNS.four <- MDNS.four  + theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm")) 
    return(list(MDNS.four,sigobs))
    
  }
  
  
  
}


##################################### FUNCTIONAL DISTINCTIVENESS ##################################### 
######################## Difference in Observed Trait Values ###########################################

#### Prune trait values down for each community and calculate difference in observed mean trait values 
#### for native and invasive species
# phy= phylogeny of species pool 
# community = community data matrix (rownames = species, colnames= communtiy names)
# traits = trait matirx, with first column labeled "Status"
# OneTrait = the name of one trait from trait matrix (rownames = community names, col= trait)
# col = column number indicating the community of interest 

summary.trait.meas  <- function(phy, community, traits, OneTrait, col){
  # Prune community
  dataPruned <- (community)[which(!community[col] == "0"), col]  ## prune tree step
  n.sp <- length(dataPruned[dataPruned == "n"]) # get all the tips which are native
  i.sp <- length(dataPruned[dataPruned == "i"]) # get all the tips which are invasive
  tot.sp <- length(dataPruned)
  if (length(dataPruned) == 1){
    summ <- rbind("NA", "NA", "NA", "NA", "NA","NA", "NA", "NA", "NA","NA", "NA", "NA",
                  "NA", "NA", "NA","NA", "NA", "NA", "NA", "NA")
    rownames(summ) <- c("Total.tips.with.trait","Percent.of.total.tips.with.trait", 
                        "Median.tot","Min.tot", "Max.tot",
                        "Total.native.tips.with.trait",  "Percent.of.native.tips.with.trait", 
                        "Mean.nat", "Median.nat", "Min.nat", "Max.nat",
                        "Total.invasive.tips.with.trait", "Percent.of.invasive.tips.with", 
                        "Mean.inv", "Median.inv", "Min.inv", "Max.inv", 
                        "p.value.mean.trait", "conf.int.Lo", "conf.int.Hi")
    colnames(summ) <- names(OneTrait)
    return(summ)
  } 
  
  tmp <- treedata(phy, dataPruned) ## check the names 
  new.phy <- tmp$phy ## prune the tree
  com.data <- tmp$data   ## subset the community data, to make sure in the right order
  tot.tips <- length(com.data)
  n.tips <- length((com.data)[com.data == "n"]) # get all the tips which are native
  i.tips <- length((com.data)[com.data == "i"]) # get all the tips which are invasive
  
  # merge Communty and trait data for one trait
  comTraittmp <- merge(com.data, OneTrait, all.x=T, by=0)
  comTrait <- comTraittmp[,2:3]
  rownames(comTrait) <- comTraittmp$Row.names
  
  # Prune trait data
  comTraitPruned <- (comTrait)[which(!comTrait[names(OneTrait)] == "NA"),] ## prune data to trait
  if (nrow(comTraitPruned) == 0){
    summ <- rbind("NA", "NA", "NA", "NA", "NA","NA", "NA", "NA", "NA","NA", "NA", "NA",
                  "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
    
    
    rownames(summ) <- c("Total.tips.with.trait","Percent.of.total.tips.with.trait", 
                        "Median.tot","Min.tot", "Max.tot",
                        "Total.native.tips.with.trait",  "Percent.of.native.tips.with.trait", 
                        "Mean.nat", "Median.nat", "Min.nat", "Max.nat",
                        "Total.invasive.tips.with.trait", "Percent.of.invasive.tips.with", 
                        "Mean.inv", "Median.inv", "Min.inv", "Max.inv", 
                        "p.value.mean.trait", "conf.int.Lo", "conf.int.Hi")
    colnames(summ) <- names(OneTrait)   
    return(summ)
  }
  i.comTraitPruned <- (comTraitPruned)[which(comTraitPruned[1] == "i"),] ## prune data to trait
  n.comTraitPruned <- (comTraitPruned)[which(comTraitPruned[1] == "n"),] ## prune data to trait
  n.sp.trait <- length(rownames(n.comTraitPruned))
  i.sp.trait <- length(rownames(i.comTraitPruned))
  tot.sp.trait <- n.sp.trait + i.sp.trait
  tot.sp.trait.percent <- tot.sp.trait/tot.tips
  n.sp.trait.percent <- n.sp.trait/tot.tips
  i.sp.trait.percent <- i.sp.trait/tot.tips
  
  if (n.sp.trait > 1 && i.sp.trait > 1) {
    comTraitPruned.median <- median(comTraitPruned[,2])
    comTraitPruned.min <- min(comTraitPruned[,2])
    comTraitPruned.max <- max(comTraitPruned[,2])
    
    n.sp.mean <- mean(n.comTraitPruned[,2])
    n.sp.median <- median(n.comTraitPruned[,2])
    n.sp.min <- min(n.comTraitPruned[,2])
    n.sp.max <- max(n.comTraitPruned[,2])
    
    i.sp.mean <- mean(i.comTraitPruned[,2])
    i.sp.median <- median(i.comTraitPruned[,2])
    i.sp.min <- min(i.comTraitPruned[,2])
    i.sp.max <- max(i.comTraitPruned[,2])
    
    t.means <- t.test(as.numeric(as.character(i.comTraitPruned[,2])), as.numeric(as.character(n.comTraitPruned[,2])), paired=F)
    
    summ <- rbind(tot.sp.trait, tot.sp.trait.percent, comTraitPruned.median, comTraitPruned.min, comTraitPruned.max,
                  n.sp.trait, n.sp.trait.percent, n.sp.mean, n.sp.median, n.sp.min, n.sp.max,
                  i.sp.trait, i.sp.trait.percent, i.sp.mean, i.sp.median, i.sp.min, i.sp.max,
                  t.means$p.value, t.means$conf.int[1], t.means$conf.int[2])
    
  } else {
    summ <- rbind("NA", "NA", "NA", "NA", "NA","NA", "NA", "NA", "NA","NA", "NA", "NA",
                  "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
    
  }
  
  rownames(summ) <- c("Total.tips.with.trait","Percent.of.total.tips.with.trait", 
                      "Median.tot","Min.tot", "Max.tot",
                      "Total.native.tips.with.trait",  "Percent.of.native.tips.with.trait", 
                      "Mean.nat", "Median.nat", "Min.nat", "Max.nat",
                      "Total.invasive.tips.with.trait", "Percent.of.invasive.tips.with", 
                      "Mean.inv", "Median.inv", "Min.inv", "Max.inv", 
                      "p.value.mean.trait", "conf.int.Lo", "conf.int.Hi")
  
  colnames(summ) <- names(OneTrait)
  return(summ)
  
}

sig.obs.functionDiversty <- function(obsdiffSum, obstraitSum, metadata, traitName, plottype=c("summary.Bar")){
  nns <- length(obstraitSum$p.value.mean.trait[which(obstraitSum$p.value.mean.trait >= 0.05)]) # 65 NS
  nsig <- length(obstraitSum$p.value.mean.trait[which(obstraitSum$p.value.mean.trait <= 0.05)]) # 9 significant differneces between means 
  sigSum <- (obstraitSum[which(obstraitSum$p.value <= 0.05),]) # 9 significant differneces between means 
  percentSig <- nsig/nisls #9/74 #0.1216216
  nat.greater <- suppressWarnings(sum(ifelse(as.numeric(obstraitSum$Median.nat) > as.numeric(obstraitSum$Median.inv), 1,0), na.rm = T)) #66
  nat.less <- suppressWarnings(sum(ifelse(as.numeric(obstraitSum$Median.nat) < as.numeric(obstraitSum$Median.inv), 1,0), na.rm = T)) #66
  
  sig.sm <- nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigSum, by=0)) # 0 small islands sig 
  sig.med <- nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigSum, by=0)) # 6/74 med islands sig  0.08108108
  sig.lg <- nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigSum, by=0)) # 3/74 large isl sig  0.04054054
  
  tr.dif <- cbind(nisls, nsig, nns, percentSig, nat.greater, nat.less, sig.sm, sig.med, sig.lg)
  rownames(tr.dif) <- traitName
  colnames(tr.dif) <- c("nIslands", "sigDiff", "ns", "%Sig", "mean_N>I", "mean_N<I", "sigSmall", "sigMed", "sigLg")  
  
  trait.Sum.ranm <- obsdiffSum[obsdiffSum$t.NNFD.p.value != "NA",]
  tmp <- transform(merge(obstraitSum, obsdiffSum, by=0), row.names=Row.names,Row.names=NULL)
  trait.tmp <- (merge(tmp, metadata, by=0))
  nisls <- nrow(trait.tmp)
  # indicate islands with significant difference in observed means between status groups for each island
  trait.tmp[,"Significance"] <- ifelse(trait.tmp[,"p.value.mean.trait"] <= 0.05, 1,
                                       ifelse(trait.tmp[,"t.NNFD.p.value"] <= 0.05, 2,
                                              ifelse(trait.tmp[,"t.MFD.p.value"] <= 0.05, 3, 0)))
  
  
  #### summarize significance data by size categories
  trait.tmp$Size.cat <- as.character(trait.tmp$Size.cat)
  trait.tmp$Size.cat <- factor(trait.tmp$Size.cat, levels=c("sm", "med", "lg"))
  
  sigMeas = cbind("Island"=trait.tmp[,1], "Size.cat"=as.character(seedmass[,"Size.cat"]), "Metric"=rep(x = "measure", times = nrow(seedmass)), "Significance"=ifelse(seedmass[,"p.value.mean.trait"] <= 0.05, 1,0))
  sigNNFD = cbind("Island"=trait.tmp[,1], "Size.cat"=as.character(trait.tmp[,"Size.cat"]), "Metric"=rep(x = "NNFD", times = nrow(trait.tmp)), "Significance"=ifelse(trait.tmp[,"t.NNFD.p.value"] <= 0.05, 1,0))
  sigMFD = cbind("Island"=trait.tmp[,1], "Size.cat"=as.character(trait.tmp[,"Size.cat"]), "Metric"=rep(x = "MFD", times = nrow(trait.tmp)), "Significance"=ifelse(trait.tmp[,"t.MFD.p.value"] <= 0.05, 1,0))
  obs.func.diff <- as.data.frame(rbind(sigMeas, sigNNFD, sigMFD))
  
  ##### significance in observed difference bewteen mean NNFDi / mean NNFDn 
  nsdifNNFD <- length(trait.Sum.ranm$t.NNFD.p.value[which(trait.Sum.ranm$t.NNFD.p.value >= 0.05)]) # 34 islands NS
  sigdifNNFD <- length(trait.Sum.ranm$t.NNFD.p.value[which(trait.Sum.ranm$t.NNFD.p.value <= 0.05)]) # 40 isalnds with significant differnece in mean betweenmean MMNPDi / mean MMNPDn 
  trait.tmp.NNFD.sig <- trait.Sum.ranm[which(trait.Sum.ranm$t.NNFD.p.value <= 0.05),]
  #trait.tmp[which(trait.tmp$t.NNFD.p.value <= 0.005),] # 3
  per.diff.NNFD <- sigdifNNFD/nisls # 0.5405405 
  nat.greater.NNFD <- sum(ifelse(as.numeric(trait.Sum.ranm$meanNNFDnatives) > as.numeric(trait.Sum.ranm$meanNNFDinvasives), 1,0)) #74
  nat.less.NNFD <- sum(ifelse(as.numeric(trait.Sum.ranm$meanNNFDnatives) < as.numeric(trait.Sum.ranm$meanNNFDinvasives), 1,0)) #74
  
  sm.sigNNFD <- nrow(merge(metadata[which(metadata$Size.cat == "sm"),], trait.tmp.NNFD.sig, by=0)) # 7/74 small islands sig NNFD 0.09459459
  med.sigNNFD <- nrow(merge(metadata[which(metadata$Size.cat == "med"),], trait.tmp.NNFD.sig, by=0)) # 20/74 medium isl sig NNFD 0.2702703
  lg.sigNNFD <- nrow(merge(metadata[which(metadata$Size.cat == "lg"),], trait.tmp.NNFD.sig, by=0)) # 13/74 large isl sig NNFD 0.1756757
  
  #### significance in observed difference bewteen mean MFDin / mean MFDnn
  nsdifMFD <- length(trait.Sum.ranm$t.MFD.p.value[which(trait.Sum.ranm$t.MFD.p.value >= 0.05)]) # 53 NS
  sigdifMFD <- length(trait.Sum.ranm$t.MFD.p.value[which(trait.Sum.ranm$t.MFD.p.value <= 0.05)]) # 21 significant differneces between mean MFD
  trait.tmp.MFD.sig <- trait.Sum.ranm[which(trait.Sum.ranm$t.MFD.p.value <= 0.05),]
  #trait.tmp[which(trait.tmp$t.MFD.p.value <= 0.005),] # 5
  per.diff.MFD <- sigdifMFD/nisls # 0.2837838
  nat.greater.MFD <-sum(ifelse(as.numeric(trait.Sum.ranm$meanMFDnatives) > as.numeric(trait.Sum.ranm$meanMFDinvasives), 1,0)) #70
  nat.less.MFD <-sum(ifelse(as.numeric(trait.Sum.ranm$meanMFDnatives) < as.numeric(trait.Sum.ranm$meanMFDinvasives), 1,0)) #4
  
  sm.sigMFD <- nrow(merge(metadata[which(metadata$Size.cat == "sm"),], trait.tmp.MFD.sig, by=0)) # 2/74 small islands sig MFD 0.02702703
  med.sigMFD <- nrow(merge(metadata[which(metadata$Size.cat == "med"),], trait.tmp.MFD.sig, by=0)) # 4/74 med islands sig MFD 0.05405405
  lg.sigMFD <- nrow(merge(metadata[which(metadata$Size.cat == "lg"),], trait.tmp.MFD.sig, by=0)) # 15/74 large isl sig MFD 0.2027027
  
  NNFD <- c(nisls, sigdifNNFD, nsdifNNFD, per.diff.NNFD, nat.greater.NNFD, nat.less.NNFD, sm.sigNNFD, med.sigNNFD, lg.sigNNFD)
  MFD <- c(nisls, sigdifMFD, nsdifMFD, per.diff.MFD, nat.greater.MFD, nat.less.MFD, sm.sigMFD, med.sigMFD, lg.sigMFD)
  
  sigobs <- rbind(tr.dif, NNFD, MFD)
  colnames(sigobs) <- c("nIslands", "sigDiff", "ns", "%Sig", "mean_N>I", "mean_N<I", "sigSmall", "sigMed", "sigLg")
  
  
  if (plottype[1] == "summary.Bar"){
    ### Bar plot of significance 
    neworder = (c("measure", "NNFD", "MFD"))
    p <- ggplot(arrange(transform(obs.func.diff, Metric=factor(Metric, levels=neworder)), Metric), aes(x=Significance, fill=factor(Size.cat))) + 
      geom_bar(stat="bin") +
      scale_fill_manual(name="Size Category",
                        breaks=c("sm", "med", "lg"),
                        values=c("sm"="dodgerblue4", "med"="orangered3", "lg"="gold1"),
                        labels=c("sm"="small", "med"="medium", "lg"="large")) +
      theme_bw() +
      scale_x_discrete(name="", labels=c("0"="NS", "1"="Significant")) +
      facet_wrap(~ Metric) +
      ggtitle(paste("Observed difference of \n Functional Distances\n", traitName))
    return(list(p,sigobs))
  } else {
    return(sigobs)
  }
  
}


#### Creates of dataframe of just species in phylo and community for that have functional trait data for one trait
# phy= phylogeny of species pool 
# community = community data matrix (rownames = species, colnames= communtiy names)
# traits = trait matirx, with first column labeled "Status"
# OneTrait = the name of one trait from trait matrix (rownames = community names, col= trait)
# col = column number indicating the community of interest 

pruneTrait <- function(phy, community, traits, OneTrait, col){
  ## prune community to include just species present
  dataPruned <- (community)[which(!community[col] == "0"), col]  
  if (length(dataPruned[dataPruned == "n"]) == 1){ # if community has only one species, return
    TMP <- as.data.frame(cbind(paste(names(community[col]), "has only one native species", sep=" "), NA,NA))
    colnames(TMP) <- c("Row.names", "Status", names(OneTrait))
    return(TMP)
  } 
  if (length(dataPruned[dataPruned == "i"]) == 1){
    TMP <- as.data.frame(cbind(paste(names(community[col]), "has only one invasive species", sep=" "), NA,NA))
    colnames(TMP) <- c("Row.names", "Status", names(OneTrait))
    return(TMP)
  } 
  if (length(dataPruned[dataPruned == "i"]) == 0){
    TMP <- as.data.frame(cbind(paste(names(community[col]), "has no invasive species", sep=" "), NA,NA))
    colnames(TMP) <- c("Row.names", "Status", names(OneTrait))
    return(TMP)
  }
  
  tmp <- treedata(phy, dataPruned) ## check the names 
  new.phy <- tmp$phy ## prune the tree
  com.data <- tmp$data   ## subset the community data, to make sure in the right order
  
  # merge Communty and trait data for one trait
  comTrait <- merge(com.data, OneTrait, all.x=T, by=0)
  
  # Prune trait data
  comTraitPruned <- (comTrait)[which(!comTrait[names(OneTrait)] == "NA"),] ## prune data to trait
  # if community has only one species, return
  if (nrow(comTraitPruned) == 0){
    TMP <- as.data.frame(cbind(paste(names(community[col]), "has no values for", names(OneTrait), sep=" "), NA,NA))
    colnames(TMP) <- c("Row.names", "Status", names(OneTrait))
    return(TMP)
  }
  if (length(comTraitPruned[comTraitPruned == "n"]) == 1){
    TMP <- as.data.frame(cbind(paste(names(community[col]), "has only one native with values for", names(OneTrait), sep=" "), NA,NA))
    colnames(TMP) <- c("Row.names", "Status", names(OneTrait))
    return(TMP)
  } 
  if (length(comTraitPruned[comTraitPruned == "i"]) == 1){
    TMP <- as.data.frame(cbind(paste(names(community[col]), "has only one invasive with values for", names(OneTrait), sep=" "), NA,NA))
    colnames(TMP) <- c("Row.names", "Status", names(OneTrait))
    return(TMP)
  } 
  if (length(comTraitPruned[comTraitPruned == "i"]) == 0){
    TMP <- as.data.frame(cbind(paste(names(community[col]), "has no invasives with values for", names(OneTrait), sep=" "), NA,NA))
    colnames(TMP) <- c("Row.names", "Status", names(OneTrait))
    return(TMP)
  }
  if (length(comTraitPruned[comTraitPruned == "n"]) == 0){
    TMP <- as.data.frame(cbind(paste(names(community[col]), "has no natives with values for", names(OneTrait), sep=" "), NA,NA))
    colnames(TMP) <- c("Row.names", "Status", names(OneTrait))
    return(TMP)
  }
  
  #i.comTraitPruned <- (comTraitPruned)[which(comTraitPruned[1] == "i"),] ## prune data to trait
  #n.comTraitPruned <- (comTraitPruned)[which(comTraitPruned[1] == "n"),] ## prune data to trait
  
  colnames(comTraitPruned) <- c("Row.names", "Status", "Trait")
  #print(ggplot(data=as.data.frame(comTraitPruned), aes(x = Trait, fill = as.character(Status))) + geom_density(alpha = 0.5))
  colnames(comTraitPruned) <- c("Row.names", "Status", names(OneTrait))
  rownames(comTraitPruned) <- comTraitPruned$Row.names
  return(comTraitPruned)
}


#### Produces a plot of the distribution of observed functional traits for each island
# list = output of pruneTrait() ; list of trait values for one trait acrross all communities
# metadata = metadata file (island names as rownames)
# meta.data.column.name = the column name of metadata to append to each element in list
# plot.title = main plot title
# y.axis.title = y axis title
melt.trait.to.meta <- function(list, metadata, meta.data.column.name, plot.title, y.axis.title){
  meta.list <- list()
  for (i in 1:length(list)){ 
    tmp <- metadata[as.character(names(list[i])), meta.data.column.name]
    newlist <- mapply(cbind, list[i], meta.data.column.name=tmp, SIMPLIFY=F) 
    meta.list[i] <- newlist
  }
  names(meta.list) <- names(list) 
  length(meta.list)
  meta.list.melt <- melt(meta.list, measure.vars=3, variable.name=4)
  meta.list.melt <- meta.list.melt[which(!meta.list.melt$value == "NA"),]
  head(meta.list.melt)
  
  p <- ggplot(meta.list.melt, aes(x=reorder(factor(L1),meta.data.column.name), y=as.numeric(as.character(value))), position=position_dodge(width=4))#, col=c("magenta1", "green3"))
  p <- p + geom_boxplot(aes(fill = as.character(Status)), width = 1)
  #sla <- sla + geom_smooth(method="lm", se=FALSE, color="black", aes(group=1))
  p <- p + scale_x_discrete("Island (increasing size)") #, breaks=seq(0, 80, 10)) 
  p <- p + scale_y_continuous(y.axis.title) 
  p <- p + theme_bw() 
  p <- p  + scale_fill_manual(values=c("i"= "magenta1", "n"="green3"), labels=c("n"="Native", "i"="Introduced")) ##breaks=rev(factor(SJnew$status)),
  p <- p + guides(fill=guide_legend(title=""))
  p <- p + theme(legend.position="top")
  p <- p + theme(axis.text.x = element_text(angle = -45, hjust = 0))
  p <- p + ggtitle(plot.title) + theme(plot.title=element_text(size=rel(1.5)))
  p
  
  
}

############################  Differnece in trait values between 1) each species and nearest Neighbor (NNFD)
############################                                     2) each species and the native community (MFD)

#### For one trait, calculates Nearest Neighbor Trait Differences (NNFD), 
## and functional diffference for each species to mean trait value of native (MFD.n) and invasive commmunity (MFD.i)
# output = results from phyloDistinct for one community
# traits = results from pruntTrait
# traitname = the name of one trait from trait matrix

functionDistinct <- function(output, traits, traitname){
  diff.data <- data.frame()
  if (nrow(output)==1){
    dat <- cbind("NA", "NA", "NA", "NA")
    diff.data <- rbind(diff.data, dat)
    colnames(diff.data) <-c("Species.Status", paste("NNFD", traitname, sep="_"), "MFD.n",  "MFD.i")
    return(diff.data)
  }
  if (nrow(as.data.frame(output[output["Species.Status"]=="n",]))==1){
    dat <- cbind("NA", "NA", "NA", "NA")
    diff.data <- rbind(diff.data, dat)
    colnames(diff.data) <- c("Species.Status", paste("NNFD", traitname, sep="_"), "MFD.n",  "MFD.i")
    return(diff.data)
  }
  if (nrow(as.data.frame(output[output["Species.Status"]=="i",]))==0){
    dat <- cbind("NA", "NA", "NA", "NA")
    diff.data <- rbind(diff.data, dat)
    colnames(diff.data) <- c("Species.Status", paste("NNFD", traitname, sep="_"), "MFD.n",  "MFD.i")
    return(diff.data)
  }
  output <- as.data.frame(output) ### DNNS results; have info about nearest native species
  
  # merge Communty and trait data for one trait: Row.names, Species.Status,  Nearest.native, leafletSize
  comDNNSTrait <- merge(output[c("Species.Status","Nearest.native")], na.omit(traits[traitname]), by=0) #merge DNNS data and trait data; remove NA
  ### IF THERE IS NO TRAIT DATA FOR NEAREST NATIVE (DNNS) <- NA....
  if (nrow(comDNNSTrait) == 0) {
    dat <- cbind("NA", "NA", "NA", "NA")
    diff.data <- rbind(diff.data, dat)
    colnames(diff.data) <- c("Species.Status", paste("NNFD", traitname, sep="_"), "MFD.n",  "MFD.i")
    return(diff.data)
  }
  
  dim(comDNNSTrait)
  i.comDNNSTrait <- (comDNNSTrait)[which(comDNNSTrait[2] == "i"),] ## just invasive species with trait on island
  n.comDNNSTrait <- (comDNNSTrait)[which(comDNNSTrait[2] == "n"),] ## native species with trait on island
  
  # create difference matrix of trait values == value of difference in trait values between each species 
  Table1 <- (comDNNSTrait[4])
  trait.dist <- t(outer(Table1[,1], Table1[,1], `-`))
  rownames(trait.dist) <- comDNNSTrait$Row.names
  colnames(trait.dist) <- comDNNSTrait$Row.names
  #dim(trait.dist)
  #head(trait.dist)
  
  for (i in 1:nrow(trait.dist)){
    i.traits.new <- i.comDNNSTrait[i.comDNNSTrait$Row.names != comDNNSTrait[i,"Row.names"], "Row.names"] # remove self comparison for invasives
    n.traits.new <- n.comDNNSTrait[n.comDNNSTrait$Row.names != comDNNSTrait[i,"Row.names"], "Row.names"]# remove self comparison for natives
    ifelse(colnames(trait.dist)[i] %in% i.comDNNSTrait$Row.names, ss <- "i", ss <- "n") #lable taxa as Native/introduced
    
    sp.name <- comDNNSTrait[i,"Row.names"]  #name of i
    
    mfd.i <- mean(trait.dist[which(rownames(trait.dist) %in% i.traits.new),i]) # mean (differnece in traits between i and each invasive species in community)
    mfd.n <- mean(trait.dist[which(rownames(trait.dist) %in% n.traits.new),i]) # mean (differnece in traits between i and each native species in community)
    
    near.n <- strsplit(as.character(comDNNSTrait[i,"Nearest.native"]), split=".", fixed=TRUE)[[1]] #get name of nearest native
    #near.i <- strsplit(as.character(comDNNSTrait[i,"Nearest.native"]), split=".", fixed=TRUE)[[1]] #get name of nearest native
    
    ## Calculate differnece in traits between species and nearest native; if near.n > 1, take median of trait  
    if (length(near.n) > 1){
      #print(c(near.n, i))
      tmp <- data.frame()
      for (k in 1:length(near.n)){
        NNFD.tmp <- comDNNSTrait[which(comDNNSTrait$Row.names==near.n[k]), traitname] # trait value nearest native
        tmp <- c(tmp, NNFD.tmp)
      }
      tmp <- as.numeric(tmp, na.rm=T)
      NNFD <- comDNNSTrait[i, traitname] - median(tmp, na.rm=T) # 0.06399687
      if (length(NNFD) == 0) {
        NNFD <- "NA"
        dat <- cbind(ss, NNFD, mfd.n, mfd.i)
        rownames(dat) <- sp.name
      } else {
        dat <- cbind(ss, NNFD, mfd.n, mfd.i)
        rownames(dat) <- sp.name
      }  
    } else{
      NNFD <- comDNNSTrait[i, traitname] -  comDNNSTrait[which(comDNNSTrait$Row.names==near.n), traitname] # trait value for i - nearest.native
      #trait.dist[i, (colnames(trait.dist) ==near.n)] # the same thing
      
      if (length(NNFD) == 0) {
        NNFD <- "NA"
        dat <- cbind(ss, NNFD, mfd.n, mfd.i)
        rownames(dat) <- sp.name
      } else {
        dat <- cbind(ss, NNFD, mfd.n, mfd.i)
        rownames(dat) <- sp.name
        
      }
      
    } 
    
    diff.data <- rbind(diff.data, dat)
  }
  rownames(diff.data) <- comDNNSTrait$Row.names
  colnames(diff.data) <- c("Species.Status", paste("NNFD", traitname, sep="_"), "MFD.n",  "MFD.i")
  diff.data
}

#### Plot functional distinctiveness of one trait for one island
# NNFDoutput = output from functionDistinct() for one trait
# islandname = name of island
# traitname = name of trait
# metric = which metric (a column number from NNFDoutput)
plot.functionDistinct.Obs <- function(NNFDoutput, islandname, traitname, metric){
  p <- qplot(data=NNFDoutput, x=factor(as.character(Species.Status)), y=as.numeric(as.character(NNFDoutput[,metric])))
  p <- p + geom_boxplot(aes(fill=factor(as.character(Species.Status))), width = 1)
  p <- p + scale_x_discrete(" ", breaks=seq(0, 80, 10)) 
  p <- p + ylab("Log NNFD (species on island)") 
  p <- p + theme_bw() 
  p <- p + scale_fill_manual(values=c("i"= "magenta1", "n"="green3"), labels=c("i"="Introduced", "n" ="Native"))
  p <- p + guides(fill=guide_legend(title=""))
  p <- p + theme(legend.position="top")
  p <- p + ggtitle(paste(islandname, traitname)) + theme(plot.title=element_text(size=rel(1.5)))
  return(p)
}

#### Plot distinctiveness of one trait, across multiple islands, by increasing island size 
## appends metadata to list all islands NNFDoutput, melts into dataframe, and output box plot
#### Specific for Mean Functional Distance Differences (NNFD)
# list = results from functionDistinct() for one community, each list element is one community
# metadata = metadata file (island names as rownames)
# meta.data.column.name = the column name of metadata to append to each element in list
# plot.title = main plot title 
# # y.axis.title = y axis title
melt.NNFD.to.meta <- function(list, metadata, meta.data.column.name, plot.title, y.axis.title){
  meta.list <- list()
  for (i in 1:length(list)){ 
    tmp <- metadata[as.character(names(list[i])), meta.data.column.name]
    newlist <- mapply(cbind, list[i], meta.data.column.name=tmp, SIMPLIFY=F) 
    meta.list[i] <- newlist
  }
  names(meta.list) <- names(list) 
  length(meta.list)
  meta.list.melt <- melt(meta.list, measure.vars=2)
  meta.list.melt <- meta.list.melt[which(!meta.list.melt$value == "NA"),]
  head(meta.list.melt)
  
  p <- ggplot(meta.list.melt, aes(x=reorder(factor(L1),meta.data.column.name), y=as.numeric(as.character(value))), position=position_dodge(width=4))#, col=c("magenta1", "green3"))
  p <- p + geom_boxplot(aes(fill = as.character(Species.Status)), width = 1)
  #sla <- sla + geom_smooth(method="lm", se=FALSE, color="black", aes(group=1))
  p <- p + scale_x_discrete("Island (increasing size)") #, breaks=seq(0, 80, 10)) 
  p <- p + scale_y_continuous(y.axis.title) 
  p <- p + theme_bw() 
  p <- p  + scale_fill_manual(values=c("i"= "magenta1", "n"="green3"), labels=c("n"="Native", "i"="Introduced")) ##breaks=rev(factor(SJnew$status)),
  p <- p + guides(fill=guide_legend(title=""))
  p <- p + theme(legend.position="top")
  p <- p + theme(axis.text.x = element_text(angle = -45, hjust = 0))
  p <- p + ggtitle(plot.title) + theme(plot.title=element_text(size=rel(1.5)))
  p <- p + theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
  p
  
}

#### Plot distinctiveness of one trait, across multiple islands, by increasing island size 
## appends metadata to list all islands functionDistinct() output, melts into dataframe, and output box plot
#### Specific for Mean Functional Distance Differences (MFD)
# list = results from functionDistinct() for one community, each list element is one community
# metadata = metadata file (island names as rownames)
# meta.data.column.name = the column name of metadata to append to each element in list
# plot.title = main plot title 
# # y.axis.title = y axis title
melt.MFD.to.meta <- function(list, metadata, meta.data.column.name, plot.title, y.axis.title){
  meta.list <- list()
  for (i in 1:length(list)){ 
    tmp <- metadata[as.character(names(list[i])), meta.data.column.name]
    newlist <- mapply(cbind, list[i], meta.data.column.name=tmp, SIMPLIFY=F) 
    meta.list[i] <- newlist
  }
  names(meta.list) <- names(list) 
  length(meta.list)
  meta.list.melt <- melt(meta.list, measure.vars=2)
  meta.list.melt <- meta.list.melt[which(!meta.list.melt$value == "NA"),]
  head(meta.list.melt)
  
  p <- ggplot(meta.list.melt, aes(x=reorder(factor(L1),meta.data.column.name), y=as.numeric(as.character(MFD.n))), position=position_dodge(width=4))#, col=c("magenta1", "green3"))
  p <- p + geom_boxplot(aes(fill = as.character(Species.Status)), width = 1)
  #sla <- sla + geom_smooth(method="lm", se=FALSE, color="black", aes(group=1))
  p <- p + scale_x_discrete("Island (increasing size)") #, breaks=seq(0, 80, 10)) 
  p <- p + scale_y_continuous(y.axis.title) 
  p <- p + theme_bw() 
  p <- p  + scale_fill_manual(values=c("i"= "magenta1", "n"="green3"), labels=c("n"="Native", "i"="Introduced")) ##breaks=rev(factor(SJnew$status)),
  p <- p + guides(fill=guide_legend(title=""))
  p <- p + theme(legend.position="top")
  p <- p + theme(axis.text.x = element_text(angle = -45, hjust = 0))
  p <- p + ggtitle(plot.title) + theme(plot.title=element_text(size=rel(1.5)))
  p <- p + theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
  p
  
  
}

#### Summarizes NNFD, MFD across each community for each trait
#outputTraitDistance = output from functionDistinct() for each trait (list of each community)
summary.function.diff <- function(outputTraitDistance){
  outputTraitDistance.naomit <- (outputTraitDistance)[which(!outputTraitDistance[2] == "NA"), ]
  summ <- data.frame()
  if (nrow(outputTraitDistance.naomit)==0){
    summ <- rbind("NA", "NA", "NA", "NA", "NA","NA", "NA", "NA", "NA","NA", "NA", "NA")
    rownames(summ) <- c("meanNNFDnatives", "meanNNFDinvasives", "n.natives", "n.invasives", "t.NNFD.p.value", "t.NNFD.conf.int.Lo","t.NNFD.conf.int.Hi", 
                        "meanMFDnatives", "meanMFDinvasives", "t.MFD.p.value", "t.MFD.conf.int.Lo","t.MFD.conf.int.Hi")  #colnames(summ) <- names(traits[col])
    return(summ)
  }
  natives <- outputTraitDistance.naomit[outputTraitDistance.naomit["Species.Status"]=="n",] # get all native
  invasives <- outputTraitDistance.naomit[outputTraitDistance.naomit["Species.Status"]=="i",] # get all native
  
  if (nrow(natives) > 1 && nrow(invasives) > 1){
    meanNNFDn <- mean(as.numeric(as.character(natives[,2])))
    meanMFDn <- mean(as.numeric(as.character(natives[,3])))
    
    meanNNFDi <- mean(as.numeric(as.character(invasives[,2])))
    meanMFDi <- mean(as.numeric(as.character(invasives[,3])))
    
    t.NNFD <- t.test(as.numeric(as.character(natives[,2])), as.numeric(as.character(invasives[,2])), paired=F)
    t.MFD <- t.test(as.numeric(as.character(natives[,3])), as.numeric(as.character(invasives[,3])), paired=F)
    
    summ <- rbind(meanNNFDn, meanNNFDi, nrow(natives), nrow(invasives), t.NNFD$p.value, t.NNFD$conf.int[1], t.NNFD$conf.int[2], meanMFDn, meanMFDi, t.MFD$p.value, t.MFD$conf.int[1], t.MFD$conf.int[2])
    
  } else if (nrow(natives) > 1 && nrow(invasives) == 1){
    meanNNFDn <- mean(as.numeric(as.character(natives[,2])))
    meanMFDn <- mean(as.numeric(as.character(natives[,3])))
    meanNNFDi <- mean(as.numeric(as.character(invasives[,2])))
    meanMFDi <- mean(as.numeric(as.character(invasives[,3])))
    
    t.NNFD <- "NA"
    t.MFD <- "NA"
    
    summ <- rbind(meanNNFDn, meanNNFDi, nrow(natives), nrow(invasives), "NA", "NA", "NA", meanMFDn, meanMFDi, "NA", "NA", "NA")
    
  } else {
    summ <- rbind("NA", "NA", "NA", "NA", "NA","NA", "NA", "NA", "NA","NA", "NA", "NA")
    
  }
  
  rownames(summ) <- c("meanNNFDnatives", "meanNNFDinvasives", "n.natives", "n.invasives", "t.NNFD.p.value", "t.NNFD.conf.int.Lo","t.NNFD.conf.int.Hi", 
                      "meanMFDnatives", "meanMFDinvasives", "t.MFD.p.value", "t.MFD.conf.int.Lo","t.MFD.conf.int.Hi")  #colnames(summ) <- names(traits[col])
  return(summ)
}


########## Randomize occurence matrix of invasive speices to create a null distributon of invaded communities
## Random = invasive presence across community (&  phylo. distance of invasive communtiy, DNNS invasives)
## Constant = native community matrix (& phylodistances of native community, DNNS/NNFD natives, MDNS/ MFD native-native)

#For the randomizations:
#- for each column in community matrix, radomly shuffle presence of invasive, keeping the observed proportion the same (native species remain unchanged)
#- repeat 1000x to create a list of community matrices
#- Then, for each column in communtiy matrix/
#--- calculate mean of observed metric
#--- iterate across list of simulated communities to generate list of simulated metrics, and mean of each 
#--- plot the distribution of simulted means. This is the null distribution.
#--- mean obs.metric - mean(simulated means metric) / SD (simulated means metric) == ses.metric
#--- compare observed mean metric ( = test statistic) to null, see if significantly different by rank .

#### Generate a list of random communities  
# phy = community phylogeny 
# com = community matrix, rownames = species, colnames = communtities (observed)
# traits = dataframe with rownames = species, colname == "Status" has "i" coded for invasive spcies, "n" for native
# N = the number of null communities to generate

randomizeCommunity <- function(phy, com, traits, N){
  community.matrix.list <- list()
  inv.com <- rownames(traits[traits$Status=="i",])
  nat.com <- rownames(traits[traits$Status=="n",])
  ###Simulation 
  rep <- 0
  while (rep < N){
    rep <- rep + 1
    sim.comm <- data.frame(row.names=rownames(com), stringsAsFactors=FALSE)
    for (i in 1:ncol(com)){
      exclude.natives <- com[which(names(com[,i]) %in% nat.com),i] # taxa that are in larger native community
      
      shuffle.tips <- com[which(!names(com[,i]) %in% nat.com),i] # get all the taxa which are not native
      u <- unlist(shuffle.tips)
      l2 <- relist(u[sample(length(u))],skeleton=shuffle.tips)
      tmp <- c(l2, exclude.natives)
      df <- data.frame(matrix(unlist(tmp), byrow=T), row.names=names(tmp))
      #new.com <- merge(com[i], t(as.data.frame(l2, stringsAsFactors=FALSE)), by=0, all.x=T)
      new.com <- merge(com[i], df, by=0, all.x=T, sort=F)
      
      #l3 <- as.matrix(new.com)
      #l3[is.na(l3)] <- "n"
      #new.com.final <- l3[,3, drop=F]
      new.com.final <- as.data.frame(new.com[, 3, drop=F], row.names=new.com[,1])
      #new.com.final <- as.data.frame(cbind(new.com[, 3], rownames= new.com[,1], drop=F))
      sim.comm <- cbind(sim.comm, new.com.final)
      
    }
    sim.comm
    colnames(sim.comm) = colnames(com) 
    community.matrix.list[[rep]] <- sim.comm
  }
  return(community.matrix.list)
}


##### Simulate DNNS, MDNS  #####################
# phy = community phylogeny
# com = observed community matrix
# island = name of community (island); = colnames in com 
# traits = trait dataset; rownames = species names, colnames[,1] = "Status"
# N = number of communities to simulate 
sim.meanDNNS.MDNS <- function(phy, com, island, traits, N){
  tmp <- treedata(phy, com) 
  new.phy <- tmp$phy ## prune the tree to include just species in larger species pool to draw from
  com.data <- as.data.frame(tmp$data)   ## subset the community data, to make sure in the right order
  #dim(com.data)
  #head(com.data)
  ## Create a distribution of ranomized communities, summarize 
  rand.inv.com <- randomizeCommunity(com=com.data[island], traits=traits, N=N)
  
  sim.mean <- data.frame()
  
  for (i in 1:length(rand.inv.com)){
    comm.island <- data.frame(lapply(rand.inv.com[[i]], as.character), stringsAsFactors=F)
    #length(which(comm.island[1] == "n"))
    rownames(comm.island) <- rownames(rand.inv.com[[i]])
    All.Sim.com.Dist <- phyloDistinct(phy=phy, community=comm.island, col=(names(comm.island)))#apply DNNS funciton across all communities ...find.NN
    All.Sim.com.DistsummaryDNNS <- summary.DNNS.MDNS(All.Sim.com.Dist) #apply summary funciton across all communities...summary.DNNS
    
    tmp <- c(All.Sim.com.DistsummaryDNNS["n.natives",], All.Sim.com.DistsummaryDNNS["n.invasives",], All.Sim.com.DistsummaryDNNS["meanDNNSinvasives",], All.Sim.com.DistsummaryDNNS["meanDNNSnatives",], All.Sim.com.DistsummaryDNNS["meanMDNSinvasives", ], All.Sim.com.DistsummaryDNNS["meanMDNSnatives", ])
    sim.mean <- rbind(sim.mean, tmp)
  }
  colnames(sim.mean) <- c("n.native.tips", "n.invasive.tips", "meanDNNSinvasives", "meanDNNSnatives", "meanMDNSinvasives", "meanMDNSnatives")
  return(sim.mean)
  
}



############ Standardized effect size of Phylogenetic Distances ######################
### null distribution ==  from randomizeCommunity.R

#### A function to summarize community data on an island 
# phy = community phylogeny
# com = observed community matrix
# col = island name
commSummary  <- function(phy, community, col){
  dataPruned <- (community)[which(!community[col] == "0"), col]  ## prune tree step
  n.sp <- length(dataPruned[dataPruned == "n"]) # get all the tips which are native
  i.sp <- length(dataPruned[dataPruned == "i"]) # get all the tips which are invasive
  tot.sp <- length(dataPruned)
  if (length(dataPruned) == 1){
    summ <- rbind(n.sp, i.sp, tot.sp, 0, 0, 0, 0, 0,0,0)
    rownames(summ) <- c("Total native species", "Total invasive species", "Total species", 
                        "Total native tips", "Total invasive tips", "Total tips",
                        "Percent native species", "Percent invasive species", "Percent native tips", "Percent invasive tips")
    colnames(summ) <- names(community[col])
    return(summ)
  } # if community has only one species, return
  #print(c(names(community[col]), "Total native species =", n.sp), quote=F) ## print the number of native species on the island
  #print(c(names(community[col]), "Total invasive species =", i.sp), quote=F) ## print the number of native species on the island
  tmp <- treedata(phy, dataPruned) ## check the names 
  new.phy <- tmp$phy ## prune the tree
  com.data <- tmp$data   ## subset the community data, to make sure in the right order
  tot.tips <- length(com.data)
  n.tips <- rownames(com.data)[com.data == "n"] # get all the tips which are native
  #print(c(names(community[col]), "Total native tips =", length(n.tips)), quote=F) ## print the number of natives after pruning to tree
  i.tips <- rownames(com.data)[com.data == "i"] # get all the tips which are invasive
  #print(c(names(community[col]), "Total invasive tips =", length(i.tips)), quote=F) #### print the number of invasive after pruning to tree
  #Percents:
  n.percent <- as.numeric(n.sp)/as.numeric(tot.sp)
  i.percent <- as.numeric(i.sp)/as.numeric(tot.sp)
  
  n.tips.percent <- as.numeric(length(n.tips))/as.numeric(tot.tips)
  i.tips.percent <- as.numeric(length(i.tips))/as.numeric(tot.tips)
  
  summ <- rbind(n.sp, i.sp, tot.sp, length(n.tips), length(i.tips), tot.tips, n.percent, i.percent, n.tips.percent, i.tips.percent)
  rownames(summ) <- c("Total native species", "Total invasive species", "Total species", 
                      "Total native tips", "Total invasive tips", "Total tips",
                      "Percent native species", "Percent invasive species", "Percent native tips", "Percent invasive tips")
  colnames(summ) <- names(community[col])
  return(summ)
  
}

#### Standardized effect size of Phylogenetic Distances
## p.rank.DNNS.inv <- min(DNNS.inv.rankLo, DNNS.inv.rankHi) / (N + 1) 
## proportion of simulated means as or more extreme than observed
# phy = community phylogeny
# com = observed community matrix
# col = island name
# simOneIsland = output from sim.meanDNNS.MDNS() ; a list of the randomized means, each element = one community
# N = number of simulations done in simOneIsland()
ses.PhyloDist <- function(phy, com, island, simOneIsland, N){
  ## Summary of community
  comsummary <- commSummary(phy=phy, community=com, col=island)
  ## Ouput of phyloDistinct
  output <- phyloDistinct(phy=phy, community=com, col=island)
  
  if (output[1,"MinDist.Nearest.native"] == "NA"){
    p <- as.data.frame(cbind(island, "has only one species", "NA","NA", "NA","NA","NA", "NA","NA","NA", "NA", "NA"))
    rownames(p) <- island
    colnames(p) <- c("island", "status", "ntax", "DNNS.obs", "DNNS.rand.mean", "DNNS.rand.sd", "DNNS.obs.rankLow", "DNNS.obs.p.Low", "DNNS.obs.rankHi", "DNNS.obs.p.Hi", "DNNS.obs.z", "runs")
    return(p)
  }
  if (length(which(output[,"Species.Status"]=="n")) == 0){
    p <- as.data.frame(cbind(island, "has no native species", "NA","NA", "NA","NA","NA", "NA","NA","NA", "NA", "NA"))
    rownames(p) <- island
    colnames(p) <- c("island", "status", "ntax", "DNNS.obs", "DNNS.rand.mean", "DNNS.rand.sd", "DNNS.obs.rankLow", "DNNS.obs.p.Low", "DNNS.obs.rankHi", "DNNS.obs.p.Hi", "DNNS.obs.z", "runs")
    return(p)
  }
  if (length(which(output[,"Species.Status"]=="i")) == 0){
    p <- as.data.frame(cbind(island, "has no invasive species", "NA","NA", "NA","NA","NA", "NA","NA","NA", "NA", "NA"))
    rownames(p) <- island
    colnames(p) <- c("island", "status", "ntax", "DNNS.obs", "DNNS.rand.mean", "DNNS.rand.sd", "DNNS.obs.rankLow", "DNNS.obs.p.Low", "DNNS.obs.rankHi", "DNNS.obs.p.Hi", "DNNS.obs.z", "runs")
    return(p)
  }
  if (length(which(output[,"Species.Status"]=="n")) == 1){
    p <- as.data.frame(cbind(island, "has only one native species", "NA","NA", "NA","NA","NA", "NA","NA","NA", "NA", "NA"))
    rownames(p) <- island
    colnames(p) <- c("island", "status", "ntax", "DNNS.obs", "DNNS.rand.mean", "DNNS.rand.sd", "DNNS.obs.rankLow", "DNNS.obs.p.Low", "DNNS.obs.rankHi", "DNNS.obs.p.Hi", "DNNS.obs.z", "runs")
    return(p)
  }
  # observed mean DNNS /MDNS for natives, invasives
  outputSummary <- summary.DNNS.MDNS(output)
  
  # simulated mean DNNS / MDNS by randomizing invasive species in each community
  #simOneIsland <- sim.meanDNNS.MDNS(phy, com, island, traits, N)
  simOneIsland[island]

  ##### ses.DNNS.invasives  
  #mean of randomization 
  meanDNNSinvasives.inv.dist <- as.numeric(na.omit(simOneIsland[island][[1]][[3]])) #meanDNNSinvasives
  mean.sim.DNNS.inv <- mean(meanDNNSinvasives.inv.dist, na.rm=T)
  median.sim.DNNS.inv <- median(meanDNNSinvasives.inv.dist, na.rm=T)
  sd.sim.DNNS.inv <- sd(meanDNNSinvasives.inv.dist, na.rm=T)
  se.sim.DNNS.inv <- mean.sim.DNNS.inv / sqrt(N)
  #observed mean
  obs.DNNS.i <- as.numeric(as.character(outputSummary["meanDNNSinvasives",]))
  #Z-value assuming that the null hypothesis is true 
  z.DNNS.inv <- (obs.DNNS.i - mean.sim.DNNS.inv) / sd.sim.DNNS.inv 
  # quantile / ranks = frequency that oberved pattern is greater than the randomized values
  DNNS.inv.rankHi <- sum(meanDNNSinvasives.inv.dist >= as.numeric(obs.DNNS.i)) #simulated means as or more extremely lower than observed
  DNNS.inv.rankLo <- sum(meanDNNSinvasives.inv.dist <= as.numeric(obs.DNNS.i))  #simulated means as or more extremely lower than observed
  # p-values of observed mean vs. distribution of randomized means 
  p.rank.DNNS.inv <- min(DNNS.inv.rankLo, DNNS.inv.rankHi) / (N + 1) #proportion of simulated means as or more extreme than observed
  p.rank.DNNS.inv.oneminusdiff <- 1 - (abs(DNNS.inv.rankLo - DNNS.inv.rankHi) / (N + 1)) #proportion of simulated means as or more extreme than observed
  p.rank.DNNS.inv.Lo <- 1 - (DNNS.inv.rankLo / (N + 1)) #proportion of simulated means as or more extreme than observed
  p.rank.DNNS.inv.Hi <- 1 - (DNNS.inv.rankHi / (N + 1)) #proportion of simulated means as or more extreme than observed
  # p-values of z-value
  one.tailed.p.z.DNNS.inv <- pnorm(-abs(z.DNNS.inv))  #p-value of your sample is the lowest alpha level you could have used for your test and still rejected the null hypothesis given your sample. 
  
  ##### ses.MDNS.inv.nat
  #mean of randomization 
  meanMDNSinvasives.inv.dist <- as.numeric(na.omit(simOneIsland[island][[1]][[5]])) #meanMDNSinvasives
  mean.sim.MDNS.inv <- mean(meanMDNSinvasives.inv.dist, na.rm=T)
  median.sim.MDNS.inv <- median(meanMDNSinvasives.inv.dist, na.rm=T)
  sd.sim.MDNS.inv <- sd(meanMDNSinvasives.inv.dist, na.rm=T)
  se.sim.MDNS.inv <- mean.sim.MDNS.inv / sqrt(N)
  #observed mean
  obs.MDNS.i <- as.numeric(as.character(outputSummary["meanMDNSinvasives",]))
  #Z-value assuming that the null hypothesis is true 
  z.MDNS.inv <- (obs.MDNS.i - mean.sim.MDNS.inv) / sd.sim.MDNS.inv 
  # quantile / ranks = frequency that oberved pattern is greater than the randomized values
  MDNS.inv.rankHi <- sum(meanMDNSinvasives.inv.dist >= as.numeric(obs.MDNS.i)) #simulated means as or more extremely lower than observed
  MDNS.inv.rankLo <- sum(meanMDNSinvasives.inv.dist <= as.numeric(obs.MDNS.i))  #simulated means as or more extremely lower than observed
  # p-values of observed mean vs. distribution of randomized means 
  p.rank.MDNS.inv <- min(MDNS.inv.rankLo, MDNS.inv.rankHi) / (N + 1) #proportion of simulated means as or more extreme than observed
  p.rank.MDNS.inv.oneminusdiff <- 1 - (abs(MDNS.inv.rankLo-MDNS.inv.rankHi) / (N + 1)) #proportion of simulated means as or more extreme than observed
  p.rank.MDNS.inv.Lo <- 1 - (MDNS.inv.rankLo / (N + 1)) #proportion of simulated means as or more extreme than observed
  p.rank.MDNS.inv.Hi <- 1 - (MDNS.inv.rankHi / (N + 1)) #proportion of simulated means as or more extreme than observed
  # p-values of z-value
  one.tailed.p.z.MDNS.inv <- pnorm(-abs(z.MDNS.inv))  #p-value of your sample is the lowest alpha level you could have used for your test and still rejected the null hypothesis given your sample. 
  
  DNNSp.i <- rbind(island, "DNNS_inv", outputSummary["n.invasives",], obs.DNNS.i, mean.sim.DNNS.inv, median.sim.DNNS.inv, sd.sim.DNNS.inv, se.sim.DNNS.inv,
                   z.DNNS.inv, DNNS.inv.rankLo, DNNS.inv.rankHi, 
                   p.rank.DNNS.inv, p.rank.DNNS.inv.oneminusdiff, p.rank.DNNS.inv.Lo, p.rank.DNNS.inv.Hi,  one.tailed.p.z.DNNS.inv,  N)
  
  MDNSp.i <- rbind(island,  "MDNS_inv_nat", outputSummary["n.invasives",], obs.MDNS.i, mean.sim.MDNS.inv, median.sim.MDNS.inv, sd.sim.MDNS.inv, se.sim.MDNS.inv, 
                  z.MDNS.inv, MDNS.inv.rankLo, MDNS.inv.rankHi,
                  p.rank.MDNS.inv, p.rank.MDNS.inv.oneminusdiff, p.rank.MDNS.inv.Lo, p.rank.MDNS.inv.Hi, one.tailed.p.z.MDNS.inv,  N)
  
  output.DNNS <- cbind(DNNSp.i, MDNSp.i)
  
  rownames(output.DNNS) <- c("island", "Metric", "n.inv.tax", "mean.obs.metric", "mean.randomized.null", "median.randomized.null", "sd.randomized.null", "se.randomized.null", 
                             "obs.z", "rankLow", "rankHi", 
                             "p.value.ranks","p.value.ranks.oneminusdiff", "p.value.ranks.lo","p.value.ranks.hi", "one.tailed.p.value.z.", "runs")

  colnames(output.DNNS) <- c(paste(island), paste(island))
  
  return(as.data.frame(output.DNNS))
  
}

sum.sesPhyloDist <- function(plottype=c("NullObsIntervalDNNS","NullObsIntervalMDNS", "ses", "summary.Bar"), simPhyloOut, metadata){
  
  ## Read in simulation file
  sdf <- read.csv(file=simPhyloOut, as.is=T, row.names=1)
  head(sdf[, 1:6])
  list.sdf <- list()
  for (i in 1:ncol(sdf))
    try({
      newlist <- (cbind(sdf[, c(1:6)]))
      head(newlist)
      colnames(newlist) <- c("n.native.tips", "n.invasive.tips", "meanDNNSinvasives", "meanDNNSnatives", "meanMDNSinvasives", "meanMDNSnatives")
      list.sdf[[i]] <- newlist
      sdf <- sdf[, -c(1:6)]
      
      
    }, silent=TRUE)
  names(list.sdf) <- SJ_islands
  
  ## remove islands that don't make sense for the null hyp.
  #remove.islands.sim <- c("All_SanJuanIslands", "Unnamed_west_of_Castle_Island") ## no null species pool, only one invasive species
  #list.sdf.new  <- list.sdf[-which(names(list.sdf) %in% remove.islands.sim)]
  #length(list.sdf.new) # 72
  
  ## Append metadata, observed means to each element in list of communities 
  #head(metadata)
  listIslands <- list.sdf#list.sdf.new
  list.meta.null.distrib.SJ <- list()
  for (i in 1:length(listIslands)){ 
    tmp <- metadata[as.character(names(listIslands[i])), "Area.m2"]
    tmp.DNNS.i <- rep(summ.DNNS.SJ[as.character(names(summ.DNNS.SJ[i]))][[1]][[2]], times=nrow(listIslands[[i]])) # observed DNNS.i
    tmp.MDNS.in <- rep(summ.DNNS.SJ[as.character(names(summ.DNNS.SJ[i]))][[1]][[9]], times=nrow(listIslands[[i]])) # observed MDNS.i to mean native community
    newlist <- mapply(cbind, "islands"=names(listIslands[i]), listIslands[i], "Area.m2"=tmp, "Obs.meanDNNSinvasives"=tmp.DNNS.i, "Obs.meanMDNSinvasives"= tmp.MDNS.in, SIMPLIFY=F)   
    list.meta.null.distrib.SJ[i] <- newlist[1]
  }
  names(list.meta.null.distrib.SJ) <- names(list.sdf)
  #head(list.meta.null.distrib.SJ[[1]])
  
  ## Melt list of simulated mean DNNS/MDNS and observed means into dataframe, and remove NA
  sim.null.distrib.melt <- melt(list.meta.null.distrib.SJ, measure.vars="islands")
  sim.null.distrib.melt <- sim.null.distrib.melt[which(!sim.null.distrib.melt$value == "NA"),] 
  #head(sim.null.distrib.melt)
  
  listIslands <- ses.SanJuan.DNNS.MDNS
  ses.DNNS.MDNS <- data.frame()
  for (i in 1:length(listIslands)){ 
    tmp1 <- metadata[as.character(names(listIslands[i])), "Area.m2"]
    tmp2 <- metadata[as.character(names(listIslands[i])), "Size.cat"]
    newlist <-cbind(t(listIslands[[i]]), "Area.m2"=tmp1, "Size.cat"=tmp2) 
    ses.DNNS.MDNS <- rbind(ses.DNNS.MDNS, newlist)
  }
  #head(ses.DNNS.MDNS)
  #write.csv(ses.DNNS.MDNS, file="ses.DNNS.MDNS.MDNSNew2.csv")
  ses.DNNS.MDNS$p.value.ranks <- as.numeric(as.character(ses.DNNS.MDNS$p.value.ranks)) 
  ses.DNNS.MDNS$obs.z <- as.numeric(as.character(ses.DNNS.MDNS$obs.z))
  ses.DNNS.MDNS <- na.omit(ses.DNNS.MDNS)
  ses.DNNS.MDNS <- ses.DNNS.MDNS[!is.infinite(as.numeric(as.character(ses.DNNS.MDNS$obs.z))),]
  
  ## add a column indicating significant islands for plotting
  sig = (ses.DNNS.MDNS[,"p.value.ranks"] <= 0.05)
  neg = (as.numeric(as.character(ses.DNNS.MDNS[,"obs.z"])) <= 0)
  ses.DNNS.MDNS <- cbind(ses.DNNS.MDNS, sig)
  
  ## code for signficance of observed means, and for direction of pattern (clusering / overdispersed)
  ses.DNNS.MDNS[,"Significance"] <- ifelse(ses.DNNS.MDNS[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.DNNS.MDNS[,"obs.z"])) <= 0, 1,
                                           ifelse(ses.DNNS.MDNS[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.DNNS.MDNS[,"obs.z"])) >= 0, 2, 
                                                  ifelse(ses.DNNS.MDNS[,"p.value.ranks"] >= 0.05 & as.numeric(as.character(ses.DNNS.MDNS[,"obs.z"])) <= 0, 3, 4)))
  
  
  #### summarize significance data by size categories
  nisls <- length(unique(ses.DNNS.MDNS$island))
  ### ses MMNPDi 
  DNNS_inv <- ses.DNNS.MDNS[which(ses.DNNS.MDNS$Metric == "DNNS_inv"), ] 
  sigdifDNNS <- nrow(DNNS_inv[which(as.numeric(as.character(DNNS_inv$p.value.ranks)) <= 0.05), ]) # 22
  nsdifDNNS <- nrow(DNNS_inv[which(as.numeric(as.character(DNNS_inv$p.value.ranks)) >= 0.05), ]) #50 NS
  per.diff.DNNS <- sigdifDNNS/nisls  # 0.3055556
  
  clustered.DNNS <- length(which(as.numeric(as.character(DNNS_inv$p.value.ranks)) <= 0.05 &
                                   as.numeric(as.character(DNNS_inv[,"obs.z"])) < 0 )) ## 22 DNNS clustered, significant
  dispersed.DNNS <- length(which(as.numeric(as.character(DNNS_inv$p.value.ranks)) <= 0.05 &
                                   as.numeric(as.character(DNNS_inv[,"obs.z"])) > 0 )) ## 0 DNNS overdispersed, significant
  
  sigDNNS_inv <- (DNNS_inv[which(as.numeric(as.character(DNNS_inv$p.value.ranks)) <= 0.05), ])
  sm.sigDNNS <- length(which(sigDNNS_inv$Size.cat == "sm")) # 10/72 small islands sig 0.1388889
  med.sigDNNS <- length(which(sigDNNS_inv$Size.cat == "med"))  # 6/72 med islands sig  0.08333333
  lg.sigDNNS <- length(which(sigDNNS_inv$Size.cat == "lg"))  # 6/72 large isl sig  0.08333333
  
  ## ses MDNSin 
  MDNS_inv_nat <- ses.DNNS.MDNS[which(ses.DNNS.MDNS$Metric == "MDNS_inv_nat"), ] 
  sigdifMDNS <- nrow(MDNS_inv_nat[which(as.numeric(as.character(MDNS_inv_nat$p.value.ranks)) <= 0.05), ]) #27 significant
  nsdifMDNS <- nrow(MDNS_inv_nat[which(as.numeric(as.character(MDNS_inv_nat$p.value.ranks)) >= 0.05), ]) #45 NS
  per.diff.MDNS <- sigdifMDNS/nisls  # 0.375
  
  clustered.MDNS <- length(which(as.numeric(as.character(MDNS_inv_nat$p.value.ranks)) <= 0.05 &
                                   as.numeric(as.character(MDNS_inv_nat[,"obs.z"])) < 0 )) ## 5 MDNS clustered, significant
  dispersed.MDNS <- length(which(as.numeric(as.character(MDNS_inv_nat$p.value.ranks)) <= 0.05 &
                                   as.numeric(as.character(MDNS_inv_nat[,"obs.z"])) > 0 )) ## 22 MDNS overdispersed, significant
  
  sigMDNS_inv <- (MDNS_inv_nat[which(as.numeric(as.character(MDNS_inv_nat$p.value.ranks)) <= 0.05), ]) #45 NS
  sm.sigMDNS <- length(which(sigMDNS_inv$Size.cat == "sm")) #7/72  0.09722222
  med.sigMDNS <-length(which(sigMDNS_inv$Size.cat == "med")) #16/72   0.2222222
  lg.sigMDNS <-length(which(sigMDNS_inv$Size.cat == "lg")) #4/72   0.05555556
  
  ses.dnns <- c(nisls, sigdifDNNS, nsdifDNNS, per.diff.DNNS, clustered.DNNS, dispersed.DNNS, sm.sigDNNS, med.sigDNNS, lg.sigDNNS)
  ses.mdns <- c(nisls, sigdifMDNS, nsdifMDNS, per.diff.MDNS, clustered.MDNS, dispersed.MDNS, sm.sigMDNS, med.sigMDNS, lg.sigMDNS)
  
  sigobs <- rbind(ses.dnns, ses.mdns)
  colnames(sigobs) <- c("nIslands", "sigDiff", "ns", "%Sig", "sig.clustered", "sig.even", "sigSmall", "sigMed", "sigLg")
  
  ###### Plot distribution of null means, with observed mean, for each island 
  if(plottype[1] == "NullObsIntervalDNNS"){
    p <- ggplot(sim.null.distrib.melt, aes(x=reorder((L1), Area.m2), y=as.numeric(as.character(meanDNNSinvasives))), 
                position=position_dodge(width=1)) 
    p <- p + geom_boxplot(aes(fill=factor(as.character(variable))), width = 1)
    p <- p + scale_x_discrete("Island (increasing size)") 
    p <- p + scale_y_log10("log Null distribution mean(DNNS)")
    p <- p + theme_bw() 
    p <- p + scale_fill_manual(values="grey", labels="")
    p <- p + guides(fill=guide_legend(title="Null distibuiton : random invasive occurrence"))
    p <- p + theme(legend.position="top")
    p <- p + geom_point(aes(y = Obs.meanDNNSinvasives), shape=1, color="magenta1") 
    p <- p + theme(axis.text.x = element_text(angle = -45, hjust = 0))
    p <- p + ggtitle("Null Distribution and observed mean(DNNS)\n Increasing Island Size") + theme(plot.title=element_text(size=rel(1.5)))
    p <- p + theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
    return(list(p, sigobs))
    
  }
  
  if(plottype[1] == "NullObsIntervalMDNS"){
    p <- ggplot(sim.null.distrib.melt, aes(x=reorder((L1), Area.m2), y=as.numeric(as.character(meanMDNSinv_nat))), 
                position=position_dodge(width=1)) 
    p <- p + geom_boxplot(aes(fill=factor(as.character(variable))), width = 1)
    p <- p + scale_x_discrete("Island (increasing size)") #, breaks=seq(0, 80, 10)) 
    p <- p + scale_y_log10("log Null distribution mean(MDNS inv-nat)")
    p <- p + theme_bw() 
    p <- p + scale_fill_manual(values="grey", labels="")
    p <- p + guides(fill=guide_legend(title=" Null distibuiton : random invasive occurrence"))
    p <- p + theme(legend.position="top")
    p <- p + geom_point(aes(y = Obs.meanMDNSinv_nat), shape=1, color="magenta1") 
    p <- p + theme(axis.text.x = element_text(angle = -45, hjust = 0))
    p <- p + ggtitle("Null Distribution and observed mean(MDNS inv-nat)\n Increasing Island Size") + theme(plot.title=element_text(size=rel(1.5)))
    p <- p + theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
    return(list(p,sigobs))
    
  }
  
  
  if(plottype[1] == "ses"){
      p <- ggplot(ses.DNNS.MDNS, aes(x=reorder(factor(island),as.numeric(as.character(Area.m2))), 
                                   y=as.numeric(as.character((obs.z))), color=factor(sig), shape=Metric, fill=factor(sig))) +
      geom_point(size=10) +
      coord_cartesian(ylim=c(-5, 5)) + 
      scale_y_continuous("standardized effect size", breaks=seq(-8, 8, 2)) +
      scale_x_discrete("Island (increasing size)") +
      scale_fill_manual(values=alpha(c("FALSE"="white", "TRUE"= "black"), .3), guide="none") +
      scale_color_manual(name =" ",values=c("FALSE"= "grey", "TRUE"="black"), guide="none") +
      scale_shape_manual(name =" ",values=c("DNNS_inv"= 21, "MDNS_inv_nat"= 24), labels=c("DNNS_inv"="DNNS i  ", "MDNS_inv_nat"="MDNS inv_nat")) +
      theme_bw() +
      theme(legend.position="top") +
      geom_abline(intercept = 0, slope = 0, colour = "grey", size = .5) +
      geom_vline(xintercept = c(18.5, 54.5)) +
      theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
      ggtitle("Significane of Phylogenetic Distances for\n invasive species to nearest native (DNNS),\n and native community (MDNS)\n") +
      theme(plot.title=element_text(size=rel(1.5))) +
      theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
    return(list(p, sigobs))
    
  }
  
  
  
  ################################## ses PHYLO summary final
  sesDNNS <- ses.DNNS.MDNS[ses.DNNS.MDNS[,"Metric"] =="DNNS_inv",]
  sesMDNS <- ses.DNNS.MDNS[ses.DNNS.MDNS[,"Metric"] !="DNNS_inv",]
  
  sigDNNSpos = cbind("Island"=rownames(sesDNNS), "Size.cat"=as.character(sesDNNS[,"Size.cat"]), 
                     "Metric"=rep(x = "DNNS inv > 0", times = nrow(sesDNNS)), 
                     "Significance"=ifelse(sesDNNS[,"p.value.ranks"] <= 0.05 & 
                                             as.numeric(as.character(sesDNNS[,"obs.z"])) >= 0, 1,0))
  sigDNNSneg = cbind("Island"=rownames(sesDNNS), "Size.cat"=as.character(sesDNNS[,"Size.cat"]), 
                     "Metric"=rep(x = "DNNS inv < 0", times = nrow(sesDNNS)), 
                     "Significance"=ifelse(sesDNNS[,"p.value.ranks"] <= 0.05 & 
                                             as.numeric(as.character(sesDNNS[,"obs.z"])) <= 0, 1,0))
  sigMDNSpos = cbind("Island"=rownames(sesMDNS), "Size.cat"=as.character(sesMDNS[,"Size.cat"]), 
                     "Metric"=rep(x = "MDNS inv > 0", times = nrow(sesMDNS)), 
                     "Significance"=ifelse(sesMDNS[,"p.value.ranks"] <= 0.05 & 
                                             as.numeric(as.character(sesMDNS[,"obs.z"])) >= 0, 1,0))
  sigMDNSneg = cbind("Island"=rownames(sesMDNS), "Size.cat"=as.character(sesMDNS[,"Size.cat"]), 
                     "Metric"=rep(x = "MDNS inv < 0", times = nrow(sesMDNS)), 
                     "Significance"=ifelse(sesMDNS[,"p.value.ranks"] <= 0.05 & 
                                             as.numeric(as.character(sesMDNS[,"obs.z"])) <= 0, 1,0))
  sesPhylo <- as.data.frame(rbind(sigDNNSpos, sigDNNSneg, sigMDNSpos, sigMDNSneg))
  
  neworder <- c("DNNS inv > 0", "MDNS inv > 0", "DNNS inv < 0", "MDNS inv < 0")
  sesPhylo2 <- arrange(transform(sesPhylo, Metric=factor(Metric,levels=neworder)),Metric)
  
  
  if(plottype[1] == "summary.Bar"){
    p <- ggplot(sesPhylo2, aes(x=Significance, fill=factor(Size.cat))) + 
      geom_bar(stat="bin") +
      scale_fill_manual(name="Size Category",
                        breaks=c("sm", "med", "lg"),
                        values=c("sm"="dodgerblue4", "med"="orangered3", "lg"="gold1")) +
      theme_bw() +
      scale_x_discrete(name="", labels=c("0"="NS", "1"="Significant")) +
      facet_wrap(~Metric, ncol = 2) +
      ggtitle("SES Phylogenetic Distances \n") 
    return(list(p, sigobs))    
  
  }

}

#################################### Simulate NNFD, MFD  #################################### 
# phy = community phylogeny
# com = observed community matrix
# island = name of community (island); = colnames in com 
# traits = trait dataset; rownames = species names, colnames[,1] = "Status"
# traitname = name of trait to calculate distance of (= colname in traits)
# N = number of communities to simulate 
sim.meanNNFD.MFD <- function(phy, com, island, traits, traitname, N){
  tmp <- treedata(phy, com) 
  new.phy <- tmp$phy ## prune the tree to include just species in larger species pool to draw from
  com.data <- as.data.frame(tmp$data)   ## subset the community data, to make sure in the right order
  #dim(com.data)
  #head(com.data)
  ## Create a distribution of ranomized communities, summarize 
  rand.inv.com <- randomizeCommunity(com=com.data[island], traits=traits, N=N)
  
  sim.mean <- data.frame()
  
  for (i in 1:length(rand.inv.com)){
    comm.island <- data.frame(lapply(rand.inv.com[[i]], as.character), stringsAsFactors=F)
    #length(which(comm.island[1] == "n"))
    rownames(comm.island) <- rownames(rand.inv.com[[i]])
    All.Sim.com.Dist <- phyloDistinct(phy=phy, community=comm.island, col=(names(comm.island)))#apply DNNS funciton across all communities ...find.NN
    All.Sim.com.Functional.Dist <- functionDistinct(output=All.Sim.com.Dist, traits=traits, traitname=traitname)
    if (nrow(All.Sim.com.Functional.Dist) == 0){
      tmp <- c("NA", "NA", "NA", "NA",  "NA", "NA")
    } else {
      All.Sim.com.Funct.summaryNNFD <- summary.function.diff(All.Sim.com.Functional.Dist) #apply summary funciton across all communities...summary.DNNS
      
      tmp <- c(All.Sim.com.Funct.summaryNNFD["n.natives",], All.Sim.com.Funct.summaryNNFD["n.invasives",], 
               All.Sim.com.Funct.summaryNNFD["meanNNFDinvasives",], All.Sim.com.Funct.summaryNNFD["meanNNFDnatives",], 
               All.Sim.com.Funct.summaryNNFD["meanMFDinvasives", ], All.Sim.com.Funct.summaryNNFD["meanMFDnat_nat", ])
    }
    
    sim.mean <- rbind(sim.mean, tmp)
  }
  colnames(sim.mean) <- c("n.native.tips", "n.invasive.tips", "meanNNFDinvasives", "meanNNFDnatives", "meanMFDinvasives", "meanMFDnat_nat")
  return(sim.mean)
  
}

############ Standardized effect size of Functional Distances ######################
### null distribution ==  randomizeCommunity()
## One islands; one trait

# phy = community phylogeny
# com=SJcommNewSim ## names of communities to simulate
# island = island name (e.g.g "Reeflet_Island")
# simOneIslandOneTrait = output of sim.meanNNFD.MFD()
# outputDNNS = phyloObs[["Reeflet_Island"]] 
# traits=SJtraitLog ## trait dataset
# traitname = trait name (e.g. "maxHeight")
# N = number of community randomizations 
ses.FunctionDist <- function(phy, com, island, simOneIslandOneTrait, outputDNNS, traits, traitname, N){
  
  ## Ouput of phyloDistinct
  output <- functionDistinct(outputDNNS, traits, traitname)
  p <- as.data.frame(rbind(island,"NNFD_inv", "NA","NA", "NA","NA","NA", "NA","NA","NA", "NA", "NA","NA", "NA"))
  q <- as.data.frame(rbind(island,"MFD_inv_nat", "NA","NA", "NA","NA","NA", "NA","NA","NA", "NA", "NA","NA", "NA"))
  
  rownames(p) <- c("island", "Metric", "n.inv.tax", "mean.obs.metric", "mean.randomized.null", "median.randomized.null", "sd.randomized.null", "se.randomized.null", 
                   "obs.z", "rankLow", "rankHi", "p.value.ranks", "one.tailed.p.value.z.", "runs")
  
  s <- cbind(p, q)
  colnames(s) <- c(paste(island), paste(island))
  
  if (length(na.omit(as.numeric(as.character(output[,2])))) == 1){
    return(s)
  } else if (length(which(output[,"Species.Status"]=="n")) == 0){
    return(s)
  } else if (length(which(output[,"Species.Status"]=="i")) == 0){
    return(s)
  } else if (length(which(output[,"Species.Status"]=="n")) == 1){
    return(s)
  } else {
    # observed mean NNFD /MFD for natives, invasives
    outputSummary <- summary.function.diff(output) 
    
    # simulated mean NNFD / MFD by randomizing invasive species in each community
    simOneIslandOneTrait[island]
    
    ##### ses.NNFD.invasives
    #mean of randomization 
    meanNNFDinvasives.inv.dist <- as.numeric(na.omit(simOneIslandOneTrait[island][[1]][[3]])) #meanNNFDinvasives
    mean.sim.NNFD.inv <- mean(meanNNFDinvasives.inv.dist, na.rm=T)
    median.sim.NNFD.inv <- median(meanNNFDinvasives.inv.dist, na.rm=T)
    sd.sim.NNFD.inv <- sd(meanNNFDinvasives.inv.dist, na.rm=T)
    se.sim.NNFD.inv <- mean.sim.NNFD.inv / sqrt(N)
    #observed mean
    obs.NNFD.i <- as.numeric(as.character(outputSummary["meanNNFDinvasives",]))
    #Z-value assuming that the null hypothesis is true 
    z.NNFD.inv <- (obs.NNFD.i - mean.sim.NNFD.inv) / sd.sim.NNFD.inv 
    # quantile / ranks = frequency that oberved pattern is greater than the randomized values
    NNFD.inv.rankHi <- sum(meanNNFDinvasives.inv.dist >= as.numeric(obs.NNFD.i)) #simulated means as or more extremely lower than observed
    NNFD.inv.rankLo <- sum(meanNNFDinvasives.inv.dist <= as.numeric(obs.NNFD.i))  #simulated means as or more extremely lower than observed
    # p-values of observed mean vs. distribution of randomized means 
    p.rank.NNFD.inv <- min(NNFD.inv.rankLo, NNFD.inv.rankHi) / (N + 1) #proportion of simulated means as or more extreme than observed
    
    # p-values of z-value
    one.tailed.p.z.NNFD.inv <- pnorm(-abs(z.NNFD.inv))  #p-value of your sample is the lowest alpha level you could have used for your test and still rejected the null hypothesis given your sample. 
    
    ##### ses.MFD.inv.nat
    #mean of randomization 
    meanMFDinvasives.inv.dist <- as.numeric(na.omit(simOneIslandOneTrait[island][[1]][[5]])) #meanMFDinvasives
    mean.sim.MFD.inv <- mean(meanMFDinvasives.inv.dist, na.rm=T)
    median.sim.MFD.inv <- median(meanMFDinvasives.inv.dist, na.rm=T)
    sd.sim.MFD.inv <- sd(meanMFDinvasives.inv.dist, na.rm=T)
    se.sim.MFD.inv <- mean.sim.MFD.inv / sqrt(N)
    #observed mean
    obs.MFD.i <- as.numeric(as.character(outputSummary["meanMFDinvasives",]))
    #Z-value assuming that the null hypothesis is true 
    z.MFD.inv <- (obs.MFD.i - mean.sim.MFD.inv) / sd.sim.MFD.inv 
    # quantile / ranks = frequency that oberved pattern is greater than the randomized values
    MFD.inv.rankHi <- sum(meanMFDinvasives.inv.dist >= as.numeric(obs.MFD.i)) 
    MFD.inv.rankLo <- sum(meanMFDinvasives.inv.dist <= as.numeric(obs.MFD.i)) 
    # p-values of observed mean vs. distribution of randomized means 
    p.rank.MFD.inv <- min(MFD.inv.rankLo, MFD.inv.rankHi) / (N + 1) #proportion of simulated means as or more extremely lower than observed
    #p.rank.MFD.inv.oneminusdiff <- 1 - (abs(MFD.inv.rankLo-MFD.inv.rankLo) / (N + 1)) #proportion of simulated means as or more extreme than observed
    
    # p-values of z-value
    one.tailed.p.z.MFD.inv <- pnorm(-abs(z.MFD.inv))  
    NNFDp.i <- rbind(island, "NNFD_inv", outputSummary["n.invasives",], obs.NNFD.i, mean.sim.NNFD.inv, median.sim.NNFD.inv, sd.sim.NNFD.inv, se.sim.NNFD.inv,
                     z.NNFD.inv, NNFD.inv.rankLo, NNFD.inv.rankHi, p.rank.NNFD.inv, one.tailed.p.z.NNFD.inv,  N)
    
    MFDp.i <- rbind(island,  "MFD_inv_nat", outputSummary["n.invasives",], obs.MFD.i, mean.sim.MFD.inv, median.sim.MFD.inv, sd.sim.MFD.inv, se.sim.MFD.inv, 
                    z.MFD.inv, MFD.inv.rankLo, MFD.inv.rankHi,p.rank.MFD.inv, one.tailed.p.z.MFD.inv,  N)
    
    output.NNFD <- cbind(NNFDp.i, MFDp.i)
    rownames(output.NNFD) <- c("island", "Metric", "n.inv.tax", "mean.obs.metric", "mean.randomized.null", "median.randomized.null", "sd.randomized.null", "se.randomized.null", 
                               "obs.z", "rankLow", "rankHi", "p.value.ranks", "one.tailed.p.value.z.", "runs")
    colnames(output.NNFD) <- c(paste(island), paste(island))
    
    return(as.data.frame(output.NNFD))
    
  }
  
  
}


# output = .csv of sim.meanNNFD.MFD()
# names of islands 
read.nullOutput <- function(sim.output, islands.sim){
  ## Load output of sim.meanNNFD.MFD()
  sdf <- read.csv(file=sim.output, as.is=T, row.names=1)
  head(sdf[, 1:6])
  list.sdf <- list()
  try(
    for (i in 1:ncol(sdf)){
      newlist <- (cbind(sdf[,c(1:6)]))
      colnames(newlist) <- c("n.native.tips", "n.invasive.tips", "meanNNFDinvasives", "meanNNFDnatives", "meanMFDinvasives", "meanMFDnat_nat")
      list.sdf[[i]] <- newlist
      sdf <- sdf[, -c(1:6)]
      
    },   silent=TRUE)
  names(list.sdf) <- islands.sim
  return(list.sdf)
  
}


# plottype = "NullInvOcc" (boxplot showing distribution of null means, with observed mean) 
#             OR "ses.allIslands" (points of z-scores for SES of each observed mean to null for each island, significant islands hightlighted)
#             OR "summary.Bar" (barplot of frequency of significant islands for NNFD / MFD clustered/overdispersed)
# sim.output = directory of .csv file for the output of the null distribution 
# islands.sim = list of the island names to use
# phyloObs = output of phyloDistinct() applies across all communities 
# traits = trait file; rownames = species, column[1] = "Status", other columns = trait values
# traitname = name of the trait to analize
# metadata = rownames = island names; one column == "Area.m2
# outputPath = file extension to save summary output
sum.sesFunctionDist <- function(plottype=c("NullObsIntervalNNFD", "NullObsIntervalMFD", "ses", "summary.Bar"), sim.output, islands.sim, phyloObs, 
                                traits, traitname, metadata, saveout = FALSE, outputPath = NULL){
  ## Observed Values
  obs.NNFD <- lapply(islands.sim, function(x) functionDistinct(output=phyloObs[[x]], traits, traitname)) 
  names(obs.NNFD) <- islands.sim 
  ## Summary Observed Values
  summ.NNFD  <- lapply(islands.sim, function(x) summary.function.diff(obs.NNFD[[x]])) 
  names(summ.NNFD) <- islands.sim 
  
  ## Read in output of means from simulated communites  (list elements = communities )
  simIslands <- read.nullOutput(sim.output, islands.sim)
  ## append a column of metadata to each element in list
  list.meta.null.distrib <- list()
  for (i in 1:length(simIslands)){ 
    tmp <- metadata[as.character(names(simIslands[i])), "Area.m2"]
    tmp.NNFD.i <- rep(summ.NNFD[as.character(names(summ.NNFD[i]))][[1]][[2]], times=nrow(simIslands[[i]]))
    tmp.MFD.in <- rep(summ.NNFD[as.character(names(summ.NNFD[i]))][[1]][[9]], times=nrow(simIslands[[i]]))
    newlist <- mapply(cbind, "islands"=names(simIslands[i]), simIslands[i], "Area.m2"=tmp, "Obs.meanNNFDinvasives"=tmp.NNFD.i, "Obs.meanMFDinvasives"= tmp.MFD.in, SIMPLIFY=F) 
    list.meta.null.distrib[i] <- newlist[1]
  }
  names(list.meta.null.distrib) <- islands.sim
  #summary(list.meta.null.distrib[1])
  #head(list.meta.null.distrib[[71]])
  
  ## melt for plotting
  sim.null.distrib.melt <- melt.list(list.meta.null.distrib, measure.vars="islands")
  sim.null.distrib.melt <- sim.null.distrib.melt[which(!sim.null.distrib.melt$value == "NA"),]
  #head(sim.null.distrib.melt)
  
  #### Summarize simualted means, standardized effect size 
  ses.SanJuan.NNFD.MFD <- lapply(islands.sim, function(x) ses.FunctionDist(phy=SJfinalTree, com=SJcommNewSim, island=x, 
                                                                           simOneIslandOneTrait=simIslands, outputDNNS=phyloObs[[x]], traits=SJtraitLog, traitname=traitname, N=1000))
  names(ses.SanJuan.NNFD.MFD) <- islands.sim
  #length(ses.SanJuan.NNFD.MFD)
  if (saveout == TRUE) {
    write.csv(ses.SanJuan.NNFD.MFD, file=outputPath)
  }
  
  listIslands <- ses.SanJuan.NNFD.MFD
  ses.SJ.NNFD <- data.frame()
  for (i in 1:length(listIslands)){ 
    if (length(listIslands[[i]]) != 2){
      newlist <- NULL
    } else {
      tmp1 <- metadata[as.character(names(listIslands[i])), "Area.m2"]
      tmp2 <- metadata[as.character(names(listIslands[i])), "Size.cat"]
      newlist <-cbind(t(listIslands[[i]]), "Area.m2"=tmp1, "Size.cat"=tmp2) 
    }
    ses.SJ.NNFD <- rbind(ses.SJ.NNFD, newlist)
    
  }
  #head(ses.SJ.NNFD)
  #dim(ses.SJ.NNFD) #144
  ses.SJ.NNFD$obs.z <- as.numeric(as.character(ses.SJ.NNFD$obs.z)) 
  ses.SJ.NNFD$obs.z[!is.finite(ses.SJ.NNFD$obs.z)] <- NA
  ses.SJ.NNFD <- na.omit(ses.SJ.NNFD)
  ses.SJ.NNFD$p.value.ranks <- as.numeric(as.character(ses.SJ.NNFD$p.value.ranks)) 
  ses.SJ.NNFD.sig <- subset(ses.SJ.NNFD, p.value.ranks <= 0.05)
  sig = (ses.SJ.NNFD[,"p.value.ranks"] <= 0.05)
  ses.SJ.NNFD <- cbind(ses.SJ.NNFD, sig)
  
  ses.SJ.NNFD[,"Significance"] <- ifelse(ses.SJ.NNFD[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.SJ.NNFD[,"obs.z"])) <= 0, 1,
                                         ifelse(ses.SJ.NNFD[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.SJ.NNFD[,"obs.z"])) >= 0, 2, 
                                                ifelse(ses.SJ.NNFD[,"p.value.ranks"] >= 0.05 & as.numeric(as.character(ses.SJ.NNFD[,"obs.z"])) <= 0, 3, 4)))
  
  
  #### summarize significance data by size categories
  ses.SJ.NNFD$Size.cat <- as.character(ses.SJ.NNFD$Size.cat)
  ses.SJ.NNFD$Size.cat <- factor(ses.SJ.NNFD$Size.cat, levels=c("sm", "med", "lg"))
  
  ##### ses MMNPDi 
  NNFD_inv <- na.omit(ses.SJ.NNFD[which(ses.SJ.NNFD$Metric == "NNFD_inv"), ] )
  sigNNFD_inv <- (NNFD_inv[which(as.numeric(as.character(NNFD_inv$p.value.ranks)) <= 0.05), ]) #7
  ## Significant by size category
  sm.NNFD <- length(which(sigNNFD_inv$Size.cat == "sm")) #3/71  0.04225352
  med.NNFD <- length(which(sigNNFD_inv$Size.cat == "med")) #3/71   0.04225352
  lg.NNFD <- length(which(sigNNFD_inv$Size.cat == "lg")) #1/71   0.01408451
  
  
  #### ses MDNSin 
  MFD_inv_nat <- ses.SJ.NNFD[which(ses.SJ.NNFD$Metric == "MFD_inv_nat"), ] 
  sigMFD_inv <- (MFD_inv_nat[which(as.numeric(as.character(MFD_inv_nat$p.value.ranks)) <= 0.05), ]) #13
  ## Significant by size category
  sm.MFD <- length(which(sigMFD_inv$Size.cat == "sm")) #2/71  0.02816901
  med.MFD <- length(which(sigMFD_inv$Size.cat == "med")) #7/71   0.09859155
  lg.MFD <- length(which(sigMFD_inv$Size.cat == "lg")) #4/71   0.05633803

  nisls <- length(unique(ses.SJ.NNFD$island))
  
  dispersed.NNFD <- length(which(as.numeric(as.character(sigNNFD_inv$obs.z)) > 0))
  clustered.NNFD <- length(which(as.numeric(as.character(sigNNFD_inv$obs.z)) < 0))
  
  dispersed.MFD <- length(which(as.numeric(as.character(sigMFD_inv$obs.z)) > 0))
  clustered.MFD <- length(which(as.numeric(as.character(sigMFD_inv$obs.z)) < 0))
  
  ses.nnfd <- c(nisls, nrow(sigNNFD_inv), abs(nisls - nrow(sigNNFD_inv)),  nrow(sigNNFD_inv) / nrow(NNFD_inv), clustered.NNFD, dispersed.NNFD, sm.NNFD, med.NNFD, lg.NNFD)
  ses.mfd <- c(nisls, nrow(sigMFD_inv), abs(nisls - nrow(sigMFD_inv)), nrow(sigMFD_inv) / nrow(MFD_inv_nat), clustered.MFD, dispersed.MFD, sm.MFD, med.MFD, lg.MFD)
  
  sigobs <- rbind(ses.nnfd, ses.mfd)
  colnames(sigobs) <- c("nIslands", "sigDiff", "ns", "%Sig", "sig.clustered", "sig.even", "sigSmall", "sigMed", "sigLg")
  
  rownames(sigobs) <- c(paste("sesNNFD", traitname), paste("sesMFD", traitname))
  
  #return(sigobs)
  

  if(plottype[1] == "NullObsIntervalNNFD"){
    p1 <- ggplot(sim.null.distrib.melt, aes(x=reorder((L1), Area.m2), y=as.numeric(as.character(meanNNFDinvasives))), 
                 position=position_dodge(width=1)) +
      geom_boxplot(aes(fill=factor(as.character(variable))), width = 1) +
      scale_x_discrete("Island (increasing size)") +
      scale_y_discrete("Null distribution mean(NNFD)") +
      theme_bw() +
      scale_fill_manual(values="grey", labels="") + 
      guides(fill=guide_legend(title=" Null distibuiton : random invasive occurrence")) +
      theme(legend.position="top") +
      geom_point(aes(y = Obs.meanNNFDinvasives), shape=1, color="magenta1") +
      theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
      ggtitle(paste("Null distribution and observed NNFD of\n", traitname, sep=" ")) + 
      theme(plot.title=element_text(size=rel(1.5))) +
      theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm")) 
    return(list(sigobs, p1))
  }
  
  if(plottype[1] == "NullObsIntervalNNFD"){
    p2 <- ggplot(sim.null.distrib.melt, aes(x=reorder((L1), Area.m2), y=as.numeric(as.character(meanMFDinvasives))), 
                 position=position_dodge(width=1)) +
      geom_boxplot(aes(fill=factor(as.character(variable))), width = 1) +
      scale_x_discrete("Island (increasing size)") +
      scale_y_log10("Null distribution mean(MFD inv-nat)") +
      theme_bw() +
      scale_fill_manual(values="grey", labels="") +
      guides(fill=guide_legend(title=" Null distibuiton : random invasive occurrence")) +
      theme(legend.position="top") +
      geom_point(aes(y = Obs.meanMPFDinv_nat), shape=1, color="magenta1") +
      theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
      ggtitle(paste("Null Distribution and observed MFD of\n", traitname, sep=" ")) + 
      theme(plot.title=element_text(size=rel(1.5))) +
      theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
    return(list(sigobs, p2))
  }

  if(plottype[1] == "ses"){
    x1 <- length(which(ses.SJ.NNFD$Size.cat == "sm"))/2 + .5
    x2 <- x1 + length(which(ses.SJ.NNFD$Size.cat == "med"))/2
    p3 <- ggplot(ses.SJ.NNFD, aes(x=reorder(factor(island),as.numeric(as.character(Area.m2))), 
                                   y=as.numeric(as.character((obs.z))), color=factor(sig), shape=Metric, fill=factor(sig))) +
      geom_point(size=10) +  
      coord_cartesian(ylim=c(-5, 5)) + 
      scale_y_continuous("standardized effect size", breaks=seq(-8, 8, 2)) +
      scale_x_discrete("Island (increasing size)") +
      scale_fill_manual(values=alpha(c("FALSE"="white", "TRUE"= "black"), .3), guide="none") +
      scale_color_manual(name =" ",values=c("FALSE"= "grey", "TRUE"="black"), guide="none") +
      scale_shape_manual(name =" ",values=c("NNFD_inv"= 21, "MFD_inv_nat"= 24), labels=c("NNFD_inv"="NNFD i  ", "MFD_inv_nat"="MFD_inv_nat")) +
      theme_bw() +
      theme(legend.position="top") +
      geom_abline(intercept = 0, slope = 0, colour = "grey", size = .5) +
      geom_vline(xintercept = c(x1, x2)) +
      theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
      ggtitle(paste("Significane of", traitname, "difference for\n invasive species to nearest native (NNFD i),\n and native community (MFD_inv_nat)\n", sep=" ")) +
      theme(plot.title=element_text(size=rel(1.5))) +
      theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
    return(list(sigobs, p3))    
  }
  
  if(plottype[1] == "summary.Bar"){
    ################################## ses FUNCTIONAL summary final #################################################################### 
    sesDNNS <- ses.SJ.NNFD[ses.SJ.NNFD[,"Metric"] =="NNFD_inv",]
    sesMDNS <- ses.SJ.NNFD[ses.SJ.NNFD[,"Metric"] =="MFD_inv_nat",]
    
    sesDNNSPos = cbind("Island"=rownames(sesDNNS), "Size.cat"=as.character(sesDNNS[,"Size.cat"]), 
                        "Metric"=rep(x = "NNFD inv > 0", times = nrow(sesDNNS)), 
                        "Significance"=ifelse(sesDNNS[,"p.value.ranks"] <= 0.05 & 
                                                as.numeric(as.character(sesDNNS[,"obs.z"])) >= 0, 1,0))
    sesDNNSNeg = cbind("Island"=rownames(sesDNNS), "Size.cat"=as.character(sesDNNS[,"Size.cat"]), 
                        "Metric"=rep(x = "NNFD inv < 0", times = nrow(sesDNNS)), 
                        "Significance"=ifelse(sesDNNS[,"p.value.ranks"] <= 0.05 & 
                                                as.numeric(as.character(sesDNNS[,"obs.z"])) <= 0, 1,0))
    sesMDNSPos = cbind("Island"=rownames(sesMDNS), "Size.cat"=as.character(sesMDNS[,"Size.cat"]), 
                      "Metric"=rep(x = "MFD inv > 0", times = nrow(sesMDNS)), 
                      "Significance"=ifelse(sesMDNS[,"p.value.ranks"] <= 0.05 & 
                                              as.numeric(as.character(sesMDNS[,"obs.z"])) >= 0, 1,0))
    sesMDNSNeg = cbind("Island"=rownames(sesMDNS), "Size.cat"=as.character(sesMDNS[,"Size.cat"]), 
                      "Metric"=rep(x = "MFD inv < 0", times = nrow(sesMDNS)), 
                      "Significance"=ifelse(sesMDNS[,"p.value.ranks"] <= 0.05 & 
                                              as.numeric(as.character(sesMDNS[,"obs.z"])) <= 0, 1,0))
    sesFunctional <- as.data.frame(rbind(sesDNNSPos, sesDNNSNeg, sesMDNSPos, sesMDNSNeg))
    
    length(which(sesMDNS[,"p.value.ranks"] <= 0.05 & 
                   as.numeric(as.character(sesMDNS[,"obs.z"])) <= 0))
    
    neworder <- c("NNFD inv > 0", "MFD inv > 0", "NNFD inv < 0", "MFD inv < 0")
    sesFunctional2 <- arrange(transform(sesFunctional, Metric=factor(Metric,levels=neworder)),Metric)
    
    p4 <- ggplot(sesFunctional2, aes(x=Significance, fill=factor(Size.cat))) + 
      geom_bar(stat="bin") +
      scale_fill_manual(name="Size Category",
                        breaks=c("sm", "med", "lg"),
                        values=c("sm"="dodgerblue4", "med"="orangered3", "lg"="gold1")) +
      theme_bw() +
      scale_x_discrete(name="", labels=c("0"="NS", "1"="Significant")) +
      facet_wrap(~Metric, ncol = 2) +
      ggtitle(paste(traitname, "SES Functional Difference\n",sep=" "))
    return(list(sigobs, p4))
  }
  

}



#######################################################################################################
## builds dataframe for new plots of observed and expected phylogenetic and functional distinctiveness 
getSesFunctionalDataframe <- function(sim.output, islands.sim, phyloObs, traits, traitname, metadata){
  obs.NNFD <- lapply(islands.sim, function(x) functionDistinct(output=phyloObs[[x]], traits, traitname)) 
  names(obs.NNFD) <- islands.sim 
  ## Summary Observed Values
  summ.NNFD  <- lapply(islands.sim, function(x) functionObsSum(obs.NNFD[[x]])) 
  names(summ.NNFD) <- islands.sim 
  
  ## Read in output of means from simulated communites  (list elements = communities )
  simIslands <- read.nullOutput(sim.output, islands.sim)
  
  ## append a column of metadata to each element in list
  list.meta.null.distrib <- list()
  for (i in 1:length(simIslands)){ 
    tmp <- metadata[as.character(names(simIslands[i])), "Area.m2"]
    tmp.NNFD.i <- rep(summ.NNFD[as.character(names(summ.NNFD[i]))][[1]][[2]], times=nrow(simIslands[[i]]))
    tmp.MFD.in <- rep(summ.NNFD[as.character(names(summ.NNFD[i]))][[1]][[9]], times=nrow(simIslands[[i]]))
    newlist <- mapply(cbind, "islands"=names(simIslands[i]), simIslands[i], "Area.m2"=tmp, "Obs.meanNNFDinvasives"=tmp.NNFD.i, "Obs.meanMPFDinv_nat"= tmp.MFD.in, SIMPLIFY=F) 
    list.meta.null.distrib[i] <- newlist[1]
  }
  names(list.meta.null.distrib) <- islands.sim
  #summary(list.meta.null.distrib[1])
  #head(list.meta.null.distrib[[71]])
  
  ## melt for plotting
  sim.null.distrib.melt <- melt.list(list.meta.null.distrib, measure.vars="islands")
  #head(sim.null.distrib.melt)
  #sim.null.distrib.melt <- sim.null.distrib.melt[which(!sim.null.distrib.melt$value == "NA"),]
  #### Summarize simualted means, standardized effect size 
  ses.SanJuan.NNFD.MFD <- lapply(islands.sim, function(x) ses.FunctionDist(phy=SJfinalTree, com=SJcommNewSim, island=x, 
                                                                           simOneIslandOneTrait=simIslands, outputDNNS=phyloObs[[x]], traits=SJtraitLog, traitname=traitname, N=1000))
  names(ses.SanJuan.NNFD.MFD) <- islands.sim
  
  listIslands <- ses.SanJuan.NNFD.MFD
  ses.SJ.NNFD <- data.frame()
  for (i in 1:length(listIslands)){ 
    if (length(listIslands[[i]]) != 2){
      newlist <- NULL
    } else {
      tmp1 <- metadata[as.character(names(listIslands[i])), "Area.m2"]
      tmp2 <- metadata[as.character(names(listIslands[i])), "Size.cat"]
      newlist <-cbind(t(listIslands[[i]]), "Area.m2"=tmp1, "Size.cat"=tmp2) 
    }
    ses.SJ.NNFD <- rbind(ses.SJ.NNFD, newlist)
    
  }
  
  
  ses.SJ.NNFD$obs.z <- as.numeric(as.character(ses.SJ.NNFD$obs.z)) 
  ses.SJ.NNFD$obs.z[!is.finite(ses.SJ.NNFD$obs.z)] <- NA
  ses.SJ.NNFD$p.value.ranks <- as.numeric(as.character(ses.SJ.NNFD$p.value.ranks)) 
  
  ses.SJ.NNFD[,"sig"] <- ifelse(ses.SJ.NNFD[,"p.value.ranks"] <= 0.05 & is.finite(ses.SJ.NNFD[,"obs.z"]) & !is.na(ses.SJ.NNFD[,"obs.z"]), TRUE, FALSE)
  ses.SJ.NNFD$obs.z[!is.finite(ses.SJ.NNFD$obs.z)] <- NA
  
  x <- ses.SJ.NNFD[ses.SJ.NNFD$Metric == "NNFD_inv",]
  x <- cbind(x, "Metric2"=rep(x = paste("NNFD_inv", traitname), times = nrow(x)))
  y <- ses.SJ.NNFD[ses.SJ.NNFD$Metric == "MFD_inv_nat",]
  y <- cbind(y, "Metric2"=rep(x = paste("MFD_inv_nat", traitname), times = nrow(y)))
  ses.SJ.NNFD.f <- rbind(x, y)
  return(ses.SJ.NNFD.f)
}

