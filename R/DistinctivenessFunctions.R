###### Distinctiveness Functions #######
###### 6 Nov 2014 ########
###### Hannah E. Marx #######

### For analysis of phylogenetic and functional distinctivness to nearest native (MNNPD/MNNFD) 
###### and mean native community (MPD/MFD)
### Also, standardized effect size compared to null model randomizing invasive species occurences
### Used for Marx et al. 2015

## Data (analysis.R)
## Phylogeny: phylo object with all species in "source pool"
## Community data: data.frame with species names matching phylogeny, presence/absence of speices in each community: 
# row = species names; columns = 1 (presnce), or 0 (absence)
## Trait data: data.frame with species names mathcing phylogeny, trait values 
# row = species names (must match format of names in phylogeny); columns = trait value

##################################### PHYLOGENETIC DISTINCTIVENESS ##################################### 
############################ Mean Nearest Native Phylogenetic Distance (MNNPD)
############################ Mean Phylogenetic Distance to Native Community (MPD)
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


#### Summarize phyloDistinct output
# output = output from phyloDistinct()
summary.MNNPD.MPD<- function(output){
  summ <- data.frame()
  if (nrow(output)==1){
    summ <- rbind("NA", "NA", "NA", "NA", "NA","NA", "NA", "NA", "NA","NA", "NA", "NA")
    rownames(summ) <- c("meanMNNPDnatives", "meanMNNPDinvasives", "n.natives", "n.invasives", "t.MNNPD.p.value", "t.MNNPD.conf.int.Lo","t.MNNPD.conf.int.Hi", 
                        "meanMPDnat_nat", "meanMPDinv_nat", "t.MPD.p.value", "t.MPD.conf.int.Lo","t.MPD.conf.int.Hi")
    return(summ)
  }
  MNNPDn <- output[output["Species.Status"]=="n",] # get all native
  MNNPDi <- output[output["Species.Status"]=="i",] # get all native
  
  if (nrow(MNNPDn) > 1 && nrow(MNNPDi) > 1) {
    meanMNNPDn <- mean(as.numeric(as.character(MNNPDn[,"MinDist.Nearest.native"])))  
    meanMNNPDi <- mean(as.numeric(as.character(MNNPDi[,"MinDist.Nearest.native"]))) 
    MPDnn <- mean(as.numeric(as.character(MNNPDn[,"MeanDist.NativeCommunity"])))
    MPDin <- mean(as.numeric(as.character(MNNPDi[,"MeanDist.NativeCommunity"])))
    
    t.MNNPD <- t.test(as.numeric(as.character(MNNPDn[,"MinDist.Nearest.native"])), as.numeric(as.character(MNNPDi[,"MinDist.Nearest.native"])), paired=F)
    t.MPD <- t.test(as.numeric(as.character(MNNPDn[,"MeanDist.NativeCommunity"])), as.numeric(as.character(MNNPDi[,"MeanDist.NativeCommunity"])), paired=F)
    
    summ <- rbind(meanMNNPDn, meanMNNPDi, nrow(MNNPDn), nrow(MNNPDi), t.MNNPD$p.value, t.MNNPD$conf.int[1], t.MNNPD$conf.int[2], MPDnn, MPDin, t.MPD$p.value, t.MPD$conf.int[1], t.MPD$conf.int[2])
    
  } else if (nrow(MNNPDn) > 1 && nrow(MNNPDi) == 1) {
    meanMNNPDn <- mean(as.numeric(as.character(MNNPDn[,"MinDist.Nearest.native"])))  
    meanMNNPDi <- mean(as.numeric(as.character(MNNPDi[,"MinDist.Nearest.native"]))) 
    MPDnn <- mean(as.numeric(as.character(MNNPDn[,"MeanDist.NativeCommunity"])))
    MPDin <- mean(as.numeric(as.character(MNNPDi[,"MeanDist.NativeCommunity"])))
    
    t.MNNPD <- "NA"
    t.MPD <- "NA"
    
    summ <- rbind(meanMNNPDn, meanMNNPDi, nrow(MNNPDn), nrow(MNNPDi), "NA", "NA", "NA", MPDnn, MPDin, "NA", "NA", "NA")
    
  } else {
    summ <- rbind("NA", "NA", "NA", "NA", "NA","NA", "NA", "NA", "NA","NA", "NA", "NA")
    
  }
  
  rownames(summ) <- c("meanMNNPDnatives", "meanMNNPDinvasives", "n.natives", "n.invasives", "t.MNNPD.p.value", "t.MNNPD.conf.int.Lo","t.MNNPD.conf.int.Hi", 
                      "meanMPDnat_nat", "meanMPDinv_nat", "t.MPD.p.value", "t.MPD.conf.int.Lo","t.MPD.conf.int.Hi")
  return(summ)
}


######################## Difference in Observed Trait Values ################################

#### Prune trait values down for each community and calculate difference in observed mean trait values for native and invasive species
# phy= phylogeny of species pool 
# community = community data matrix (rownames = species, colnames= communtiy names)
# traits = trait matirx, with first column labeled "Status"
# OneTrait = the name of one trait from trait matrix (rownames = community names, col= trait)
# col = column number indicating the community of interest 
commTraitSummary  <- function(phy, community, traits, OneTrait, col){
  # Prune community
  dataPruned <- (community)[which(!community[col] == "0"), col]  ## prune tree step
  n.sp <- length(dataPruned[dataPruned == "n"]) # get all the tips which are native
  i.sp <- length(dataPruned[dataPruned == "i"]) # get all the tips which are invasive
  tot.sp <- length(dataPruned)
  if (length(dataPruned) == 1){
    summ <- rbind("NA", "NA", "NA", "NA", "NA","NA", "NA", "NA", "NA","NA", "NA", "NA",
                  "NA", "NA", "NA","NA", "NA", "NA")
    
    rownames(summ) <- c(paste("Total.tips.with", names(OneTrait), sep="."),paste("Percent.of.total.tips.with", names(OneTrait), sep="."), 
                        paste("Median.tot", names(OneTrait), sep="."),paste("Min.tot", names(OneTrait), sep="."), paste("Max.tot", names(OneTrait),sep="."),
                        paste("Total.native.tipswith", names(OneTrait), sep="."),  paste("Percent.of.native.tips.with", names(OneTrait), sep="."), 
                        paste("Median.nat", names(OneTrait), sep="."), paste("Min.nat", names(OneTrait), sep="."), paste("Max.nat", names(OneTrait), sep="."),
                        paste("Total.invasive.tips.with",names(OneTrait), sep="."), paste("Percent.of.invasive.tips.with", names(OneTrait), sep="."), 
                        paste("Median.inv",names(OneTrait), sep="."), paste("Min.inv", names(OneTrait), sep="."), paste("Max.inv", names(OneTrait), sep="."), 
                        paste("p.value", names(OneTrait), sep="."), paste("conf.int.Lo", names(OneTrait), sep="."), paste("conf.int.Hi", names(OneTrait), sep="."))
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
                  "NA", "NA", "NA", "NA", "NA", "NA")
    
    rownames(summ) <- c(paste("Total.tips.with", names(OneTrait), sep="."),paste("Percent.of.total.tips.with", names(OneTrait), sep="."), 
                        paste("Median.tot", names(OneTrait), sep="."),paste("Min.tot", names(OneTrait), sep="."), paste("Max.tot", names(OneTrait),sep="."),
                        paste("Total.native.tipswith", names(OneTrait), sep="."),  paste("Percent.of.native.tips.with", names(OneTrait), sep="."), 
                        paste("Median.nat", names(OneTrait), sep="."), paste("Min.nat", names(OneTrait), sep="."), paste("Max.nat", names(OneTrait), sep="."),
                        paste("Total.invasive.tips.with",names(OneTrait), sep="."), paste("Percent.of.invasive.tips.with", names(OneTrait), sep="."), 
                        paste("Median.inv",names(OneTrait), sep="."), paste("Min.inv", names(OneTrait), sep="."), paste("Max.inv", names(OneTrait), sep="."), 
                        paste("p.value", names(OneTrait), sep="."), paste("conf.int.Lo", names(OneTrait), sep="."), paste("conf.int.Hi", names(OneTrait), sep="."))
    
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
    
    n.sp.median <- median(n.comTraitPruned[,2])
    n.sp.min <- min(n.comTraitPruned[,2])
    n.sp.max <- max(n.comTraitPruned[,2])
    
    i.sp.median <- median(i.comTraitPruned[,2])
    i.sp.min <- min(i.comTraitPruned[,2])
    i.sp.max <- max(i.comTraitPruned[,2])
    
    t.means <- t.test(as.numeric(as.character(i.comTraitPruned[,2])), as.numeric(as.character(n.comTraitPruned[,2])), paired=F)
    
    summ <- rbind(tot.sp.trait, tot.sp.trait.percent, comTraitPruned.median, comTraitPruned.min, comTraitPruned.max,
                  n.sp.trait, n.sp.trait.percent, n.sp.median, n.sp.min, n.sp.max,
                  i.sp.trait, i.sp.trait.percent, i.sp.median, i.sp.min, i.sp.max,
                  t.means$p.value, t.means$conf.int[1], t.means$conf.int[2])
    
  } else {
    summ <- rbind("NA", "NA", "NA", "NA", "NA","NA", "NA", "NA", "NA","NA", "NA", "NA",
                  "NA", "NA", "NA", "NA", "NA", "NA")
    
  }
  
  rownames(summ) <- c(paste("Total.tips.with", names(OneTrait), sep="."),paste("Percent.of.total.tips.with", names(OneTrait), sep="."), 
                      paste("Median.tot", names(OneTrait), sep="."),paste("Min.tot", names(OneTrait), sep="."), paste("Max.tot", names(OneTrait),sep="."),
                      paste("Total.native.tipswith", names(OneTrait), sep="."),  paste("Percent.of.native.tips.with", names(OneTrait), sep="."), 
                      paste("Median.nat", names(OneTrait), sep="."), paste("Min.nat", names(OneTrait), sep="."), paste("Max.nat", names(OneTrait), sep="."),
                      paste("Total.invasive.tips.with",names(OneTrait), sep="."), paste("Percent.of.invasive.tips.with", names(OneTrait), sep="."), 
                      paste("Median.inv",names(OneTrait), sep="."), paste("Min.inv", names(OneTrait), sep="."), paste("Max.inv", names(OneTrait), sep="."), 
                      paste("p.value", names(OneTrait), sep="."), paste("conf.int.Lo", names(OneTrait), sep="."), paste("conf.int.Hi", names(OneTrait), sep="."))
  
  colnames(summ) <- names(OneTrait)
  return(summ)
  
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

############################  Differnece in trait values between 1) each species and nearest Neighbor (MNNFD)
############################                                     2) each species and the native community (MFD)

#### For one trait, calculates Nearest Neighbor Trait Differences (MNNFD), 
## and functional diffference for each species to mean trait value of native (MFD.n) and invasive commmunity (MFD.i)
# output = results from phyloDistinct for one community
# traits = results from pruntTrait
# OneTrait = the name of one trait from trait matrix
# col = column of the community of interest 
functionDistinct <- function(output, traits, traitname){
  diff.data <- data.frame()
  if (nrow(output)==1){
    dat <- cbind("NA", "NA", "NA", "NA")
    diff.data <- rbind(diff.data, dat)
    colnames(diff.data) <-c("Species.Status", paste("MNNFD", traitname, sep="_"), "MFD.n",  "MFD.i")
    return(diff.data)#(print("has only one species; go fuck yourself"))
  }
  if (nrow(as.data.frame(output[output["Species.Status"]=="n",]))==1){
    dat <- cbind("NA", "NA", "NA", "NA")
    diff.data <- rbind(diff.data, dat)
    colnames(diff.data) <- c("Species.Status", paste("MNNFD", traitname, sep="_"), "MFD.n",  "MFD.i")
    return(diff.data)#(print("has only one native species; go fuck yourself"))
  }
  if (nrow(as.data.frame(output[output["Species.Status"]=="i",]))==0){
    dat <- cbind("NA", "NA", "NA", "NA")
    diff.data <- rbind(diff.data, dat)
    colnames(diff.data) <- c("Species.Status", paste("MNNFD", traitname, sep="_"), "MFD.n",  "MFD.i")
    return(diff.data)#(print("has only one native species; go fuck yourself"))
  }
  output <- as.data.frame(output) ### MNNPD results; have info about nearest native species
  
  # merge Communty and trait data for one trait: Row.names, Species.Status,  Nearest.native, leafletSize
  comMNNPDTrait <- merge(output[c("Species.Status","Nearest.native")], na.omit(traits[traitname]), by=0) #merge MNNPD data and trait data; remove NA
  if (nrow(comMNNPDTrait) == 0) {
    dat <- cbind("NA", "NA", "NA", "NA")
    diff.data <- rbind(diff.data, dat)
    colnames(diff.data) <- c("Species.Status", paste("MNNFD", traitname, sep="_"), "MFD.n",  "MFD.i")
    return(diff.data)#(print("has only one native species; go fuck yourself"))
  }
  
  dim(comMNNPDTrait)
  i.comMNNPDTrait <- (comMNNPDTrait)[which(comMNNPDTrait[2] == "i"),] ## just invasive species with trait on island
  n.comMNNPDTrait <- (comMNNPDTrait)[which(comMNNPDTrait[2] == "n"),] ## native species with trait on island
  
  # create distance matrix of trait values == absolute value of difference in trait values between each species 
  trait.dist <- as.matrix(dist(comMNNPDTrait[4], upper=T))
  rownames(trait.dist) <- comMNNPDTrait$Row.names
  colnames(trait.dist) <- comMNNPDTrait$Row.names
  dim(trait.dist)
  head(trait.dist)
  
  for (i in 1:nrow(trait.dist)){
    i.traits.new <- i.comMNNPDTrait[i.comMNNPDTrait$Row.names != comMNNPDTrait[i,"Row.names"], "Row.names"] # remove self comparison for invasives
    n.traits.new <- n.comMNNPDTrait[n.comMNNPDTrait$Row.names != comMNNPDTrait[i,"Row.names"], "Row.names"]# remove self comparison for natives
    ifelse(colnames(trait.dist)[i] %in% i.comMNNPDTrait$Row.names, ss <- "i", ss <- "n") #lable taxa as Native/introduced
    
    sp.name <- comMNNPDTrait[i,"Row.names"]  #name of i
    
    mfd.i <- mean(trait.dist[which(rownames(trait.dist) %in% i.traits.new),i]) # mean (differnece in traits between i and each invasive species in community)
    mfd.n <- mean(trait.dist[which(rownames(trait.dist) %in% n.traits.new),i]) # mean (differnece in traits between i and each native species in community)
    
    near.n <- strsplit(as.character(comMNNPDTrait[i,"Nearest.native"]), split=".", fixed=TRUE)[[1]] #get name of nearest native
    #near.i <- strsplit(as.character(comMNNPDTrait[i,"Nearest.native"]), split=".", fixed=TRUE)[[1]] #get name of nearest native
    
    ## Calculate differnece in traits between species and nearest native; if near.n > 1, take median of trait  
    if (length(near.n) > 1){
      #print(c(near.n, i))
      tmp <- data.frame()
      for (k in 1:length(near.n)){
        MNNFD.tmp <- comMNNPDTrait[which(comMNNPDTrait$Row.names==near.n[k]), traitname] # trait value nearest native
        tmp <- c(tmp, MNNFD.tmp)
      }
      tmp <- as.numeric(tmp, na.rm=T)
      MNNFD <- comMNNPDTrait[i, traitname] - median(tmp, na.rm=T) # 0.06399687
      if (length(MNNFD) == 0) {
        MNNFD <- "NA"
        dat <- cbind(ss, MNNFD, mfd.n, mfd.i)
        rownames(dat) <- sp.name
      } else {
        dat <- cbind(ss, MNNFD, mfd.n, mfd.i)
        rownames(dat) <- sp.name
      }  
    } else{
      MNNFD <- comMNNPDTrait[i, traitname] -  comMNNPDTrait[which(comMNNPDTrait$Row.names==near.n), traitname] # trait value for i - nearest.native
      #trait.dist[i, (colnames(trait.dist) ==near.n)] # the same thing
      
      if (length(MNNFD) == 0) {
        MNNFD <- "NA"
        dat <- cbind(ss, MNNFD, mfd.n, mfd.i)
        rownames(dat) <- sp.name
      } else {
        dat <- cbind(ss, MNNFD, mfd.n, mfd.i)
        rownames(dat) <- sp.name
        
      }
      
    } 
    
    diff.data <- rbind(diff.data, dat)
  }
  rownames(diff.data) <- comMNNPDTrait$Row.names
  colnames(diff.data) <- c("Species.Status", paste("MNNFD", traitname, sep="_"), "MFD.n",  "MFD.i")
  diff.data
}

#### Plot functional distinctiveness of one trait for one island
# MNNFDoutput = output from functionDistinct() for one trait
# islandname = name of island
# traitname = name of trait
# metric = which metric (a column number from MNNFDoutput)
plot.functionDistinct.Obs <- function(MNNFDoutput, islandname, traitname, metric){
  p <- qplot(data=MNNFDoutput, x=factor(as.character(Species.Status)), y=as.numeric(as.character(MNNFDoutput[,metric])))
  p <- p + geom_boxplot(aes(fill=factor(as.character(Species.Status))), width = 1)
  p <- p + scale_x_discrete(" ", breaks=seq(0, 80, 10)) 
  p <- p + ylab("Log MNNFD (species on island)") 
  p <- p + theme_bw() 
  p <- p + scale_fill_manual(values=c("i"= "magenta1", "n"="green3"), labels=c("i"="Introduced", "n" ="Native"))
  p <- p + guides(fill=guide_legend(title=""))
  p <- p + theme(legend.position="top")
  p <- p + ggtitle(paste(islandname, traitname)) + theme(plot.title=element_text(size=rel(1.5)))
  return(p)
}

#### Plot distinctiveness of one trait, across multiple islands, by increasing island size 
## appends metadata to list all islands MNNFDoutput, melts into dataframe, and output box plot
#### Specific for Mean Functional Distance Differences (MNNFD)
# list = results from functionDistinct() for one community, each list element is one community
# metadata = metadata file (island names as rownames)
# meta.data.column.name = the column name of metadata to append to each element in list
# plot.title = main plot title 
# # y.axis.title = y axis title
melt.MNNFD.to.meta <- function(list, metadata, meta.data.column.name, plot.title, y.axis.title){
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

#### Summarizes MNNFD, MFD across each community for each trait
#outputTraitDistance = output from functionDistinct() for each trait (list of each community)
functionObsSum <- function(outputTraitDistance){
  outputTraitDistance.naomit <- (outputTraitDistance)[which(!outputTraitDistance[2] == "NA"), ]
  summ <- data.frame()
  if (nrow(outputTraitDistance.naomit)==0){
    summ <- rbind("NA", "NA", "NA", "NA", "NA","NA", "NA", "NA", "NA","NA", "NA", "NA")
    rownames(summ) <- c("meanMNNFDnatives", "meanMNNFDinvasives", "n.natives", "n.invasives", "t.MNNFD.p.value", "t.MNNFD.conf.int.Lo","t.MNNFD.conf.int.Hi", 
                        "meanMFDnat_nat", "meanMFDinv_nat", "t.MFD.p.value", "t.MFD.conf.int.Lo","t.MFD.conf.int.Hi")  #colnames(summ) <- names(traits[col])
    return(summ)
  }
  natives <- outputTraitDistance.naomit[outputTraitDistance.naomit["Species.Status"]=="n",] # get all native
  invasives <- outputTraitDistance.naomit[outputTraitDistance.naomit["Species.Status"]=="i",] # get all native
  
  if (nrow(natives) > 1 && nrow(invasives) > 1){
    meanMNNFDn <- mean(as.numeric(as.character(natives[,2])))
    meanMFDn <- mean(as.numeric(as.character(natives[,3])))
    
    meanMNNFDi <- mean(as.numeric(as.character(invasives[,2])))
    meanMFDi <- mean(as.numeric(as.character(invasives[,3])))
    
    t.MNNFD <- t.test(as.numeric(as.character(natives[,2])), as.numeric(as.character(invasives[,2])), paired=F)
    t.MFD <- t.test(as.numeric(as.character(natives[,3])), as.numeric(as.character(invasives[,3])), paired=F)
    
    summ <- rbind(meanMNNFDn, meanMNNFDi, nrow(natives), nrow(invasives), t.MNNFD$p.value, t.MNNFD$conf.int[1], t.MNNFD$conf.int[2], meanMFDn, meanMFDi, t.MFD$p.value, t.MFD$conf.int[1], t.MFD$conf.int[2])
    
  } else if (nrow(natives) > 1 && nrow(invasives) == 1){
    meanMNNFDn <- mean(as.numeric(as.character(natives[,2])))
    meanMFDn <- mean(as.numeric(as.character(natives[,3])))
    meanMNNFDi <- mean(as.numeric(as.character(invasives[,2])))
    meanMFDi <- mean(as.numeric(as.character(invasives[,3])))
    
    t.MNNFD <- "NA"
    t.MFD <- "NA"
    
    summ <- rbind(meanMNNFDn, meanMNNFDi, nrow(natives), nrow(invasives), "NA", "NA", "NA", meanMFDn, meanMFDi, "NA", "NA", "NA")
    
  } else {
    summ <- rbind("NA", "NA", "NA", "NA", "NA","NA", "NA", "NA", "NA","NA", "NA", "NA")
    
  }
  
  rownames(summ) <- c("meanMNNFDnatives", "meanMNNFDinvasives", "n.natives", "n.invasives", "t.MNNFD.p.value", "t.MNNFD.conf.int.Lo","t.MNNFD.conf.int.Hi", 
                      "meanMFDnat_nat", "meanMFDinv_nat", "t.MFD.p.value", "t.MFD.conf.int.Lo","t.MFD.conf.int.Hi")  #colnames(summ) <- names(traits[col])
  return(summ)
}


########## Randomize occurence matrix of invasive speices to create a null distributon of invaded communities
## Random = invasive presence across community (&  phylo. distance of invasive communtiy, MNNPD invasives)
## Constant = native community matrix (& phylodistances of native community, MNNPD/MNNFD natives, MPD/ MFD native-native)

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


##### Simulate MNNPD, MPD  #####################
# phy = community phylogeny
# com = observed community matrix
# island = name of community (island); = colnames in com 
# traits = trait dataset; rownames = species names, colnames[,1] = "Status"
# N = number of communities to simulate 
sim.meanMNNPD.MPD <- function(phy, com, island, traits, N){
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
    All.Sim.com.Dist <- phyloDistinct(phy=phy, community=comm.island, col=(names(comm.island)))#apply MNNPD funciton across all communities ...find.NN
    All.Sim.com.DistsummaryMNNPD <- summary.MNNPD.MPD(All.Sim.com.Dist) #apply summary funciton across all communities...summary.MNNPD
    
    tmp <- c(All.Sim.com.DistsummaryMNNPD["n.natives",], All.Sim.com.DistsummaryMNNPD["n.invasives",], All.Sim.com.DistsummaryMNNPD["meanMNNPDinvasives",], All.Sim.com.DistsummaryMNNPD["meanMNNPDnatives",], All.Sim.com.DistsummaryMNNPD["meanMPDinv_nat", ], All.Sim.com.DistsummaryMNNPD["meanMPDnat_nat", ])
    sim.mean <- rbind(sim.mean, tmp)
  }
  colnames(sim.mean) <- c("n.native.tips", "n.invasive.tips", "meanMNNPDinvasives", "meanMNNPDnatives", "meanMPDinv_nat", "meanMPDnat_nat")
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
    #return(print(paste(names(community[col]), "has only one native species; go fuck yourself", sep=" ")))
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
# simOneIsland = output from sim.meanMNNPD.MPD() ; a list of the randomized means, each element = one community
# N = number of simulations done in simOneIsland()
ses.PhyloDist <- function(phy, com, island, simOneIsland, N){
  ## Summary of community
  comsummary <- commSummary(phy=phy, community=com, col=island)
  ## Ouput of phyloDistinct
  output <- phyloDistinct(phy=phy, community=com, col=island)
  
  if (output[1,"MinDist.Nearest.native"] == "NA"){
    p <- as.data.frame(cbind(island, "has only one species", "NA","NA", "NA","NA","NA", "NA","NA","NA", "NA", "NA"))
    rownames(p) <- island
    colnames(p) <- c("island", "status", "ntax", "MNNPD.obs", "MNNPD.rand.mean", "MNNPD.rand.sd", "MNNPD.obs.rankLow", "MNNPD.obs.p.Low", "MNNPD.obs.rankHi", "MNNPD.obs.p.Hi", "MNNPD.obs.z", "runs")
    return(p)
  }
  if (length(which(output[,"Species.Status"]=="n")) == 0){
    p <- as.data.frame(cbind(island, "has no native species", "NA","NA", "NA","NA","NA", "NA","NA","NA", "NA", "NA"))
    rownames(p) <- island
    colnames(p) <- c("island", "status", "ntax", "MNNPD.obs", "MNNPD.rand.mean", "MNNPD.rand.sd", "MNNPD.obs.rankLow", "MNNPD.obs.p.Low", "MNNPD.obs.rankHi", "MNNPD.obs.p.Hi", "MNNPD.obs.z", "runs")
    return(p)
  }
  if (length(which(output[,"Species.Status"]=="i")) == 0){
    p <- as.data.frame(cbind(island, "has no invasive species", "NA","NA", "NA","NA","NA", "NA","NA","NA", "NA", "NA"))
    rownames(p) <- island
    colnames(p) <- c("island", "status", "ntax", "MNNPD.obs", "MNNPD.rand.mean", "MNNPD.rand.sd", "MNNPD.obs.rankLow", "MNNPD.obs.p.Low", "MNNPD.obs.rankHi", "MNNPD.obs.p.Hi", "MNNPD.obs.z", "runs")
    return(p)
  }
  if (length(which(output[,"Species.Status"]=="n")) == 1){
    p <- as.data.frame(cbind(island, "has only one native species", "NA","NA", "NA","NA","NA", "NA","NA","NA", "NA", "NA"))
    rownames(p) <- island
    colnames(p) <- c("island", "status", "ntax", "MNNPD.obs", "MNNPD.rand.mean", "MNNPD.rand.sd", "MNNPD.obs.rankLow", "MNNPD.obs.p.Low", "MNNPD.obs.rankHi", "MNNPD.obs.p.Hi", "MNNPD.obs.z", "runs")
    return(p)
  }
  # observed mean MNNPD /MPD for natives, invasives
  outputSummary <- summary.MNNPD.MPD(output)
  
  # simulated mean MNNPD / MPD by randomizing invasive species in each community
  #simOneIsland <- sim.meanMNNPD.MPD(phy, com, island, traits, N)
  simOneIsland[island]

  ##### ses.MNNPD.invasives  
  #mean of randomization 
  meanMNNPDinvasives.inv.dist <- as.numeric(na.omit(simOneIsland[island][[1]][[3]])) #meanMNNPDinvasives
  mean.sim.MNNPD.inv <- mean(meanMNNPDinvasives.inv.dist, na.rm=T)
  median.sim.MNNPD.inv <- median(meanMNNPDinvasives.inv.dist, na.rm=T)
  sd.sim.MNNPD.inv <- sd(meanMNNPDinvasives.inv.dist, na.rm=T)
  se.sim.MNNPD.inv <- mean.sim.MNNPD.inv / sqrt(N)
  #observed mean
  obs.MNNPD.i <- as.numeric(as.character(outputSummary["meanMNNPDinvasives",]))
  #Z-value assuming that the null hypothesis is true 
  z.MNNPD.inv <- (obs.MNNPD.i - mean.sim.MNNPD.inv) / sd.sim.MNNPD.inv 
  # quantile / ranks = frequency that oberved pattern is greater than the randomized values
  MNNPD.inv.rankHi <- sum(meanMNNPDinvasives.inv.dist >= as.numeric(obs.MNNPD.i)) #simulated means as or more extremely lower than observed
  MNNPD.inv.rankLo <- sum(meanMNNPDinvasives.inv.dist <= as.numeric(obs.MNNPD.i))  #simulated means as or more extremely lower than observed
  # p-values of observed mean vs. distribution of randomized means 
  p.rank.MNNPD.inv <- min(MNNPD.inv.rankLo, MNNPD.inv.rankHi) / (N + 1) #proportion of simulated means as or more extreme than observed
  p.rank.MNNPD.inv.oneminusdiff <- 1 - (abs(MNNPD.inv.rankLo - MNNPD.inv.rankHi) / (N + 1)) #proportion of simulated means as or more extreme than observed
  p.rank.MNNPD.inv.Lo <- 1 - (MNNPD.inv.rankLo / (N + 1)) #proportion of simulated means as or more extreme than observed
  p.rank.MNNPD.inv.Hi <- 1 - (MNNPD.inv.rankHi / (N + 1)) #proportion of simulated means as or more extreme than observed
  # p-values of z-value
  one.tailed.p.z.MNNPD.inv <- pnorm(-abs(z.MNNPD.inv))  #p-value of your sample is the lowest alpha level you could have used for your test and still rejected the null hypothesis given your sample. 
  
  ##### ses.MPD.inv.nat
  #mean of randomization 
  meanMPDinvasives.inv.dist <- as.numeric(na.omit(simOneIsland[island][[1]][[5]])) #meanMPDinvasives
  mean.sim.MPD.inv <- mean(meanMPDinvasives.inv.dist, na.rm=T)
  median.sim.MPD.inv <- median(meanMPDinvasives.inv.dist, na.rm=T)
  sd.sim.MPD.inv <- sd(meanMPDinvasives.inv.dist, na.rm=T)
  se.sim.MPD.inv <- mean.sim.MPD.inv / sqrt(N)
  #observed mean
  obs.MPD.i <- as.numeric(as.character(outputSummary["meanMPDinv_nat",]))
  #Z-value assuming that the null hypothesis is true 
  z.MPD.inv <- (obs.MPD.i - mean.sim.MPD.inv) / sd.sim.MPD.inv 
  # quantile / ranks = frequency that oberved pattern is greater than the randomized values
  MPD.inv.rankHi <- sum(meanMPDinvasives.inv.dist >= as.numeric(obs.MPD.i)) #simulated means as or more extremely lower than observed
  MPD.inv.rankLo <- sum(meanMPDinvasives.inv.dist <= as.numeric(obs.MPD.i))  #simulated means as or more extremely lower than observed
  # p-values of observed mean vs. distribution of randomized means 
  p.rank.MPD.inv <- min(MPD.inv.rankLo, MPD.inv.rankHi) / (N + 1) #proportion of simulated means as or more extreme than observed
  p.rank.MPD.inv.oneminusdiff <- 1 - (abs(MPD.inv.rankLo-MPD.inv.rankHi) / (N + 1)) #proportion of simulated means as or more extreme than observed
  p.rank.MPD.inv.Lo <- 1 - (MPD.inv.rankLo / (N + 1)) #proportion of simulated means as or more extreme than observed
  p.rank.MPD.inv.Hi <- 1 - (MPD.inv.rankHi / (N + 1)) #proportion of simulated means as or more extreme than observed
  # p-values of z-value
  one.tailed.p.z.MPD.inv <- pnorm(-abs(z.MPD.inv))  #p-value of your sample is the lowest alpha level you could have used for your test and still rejected the null hypothesis given your sample. 
  
  MNNPDp.i <- rbind(island, "MNNPD_inv", outputSummary["n.invasives",], obs.MNNPD.i, mean.sim.MNNPD.inv, median.sim.MNNPD.inv, sd.sim.MNNPD.inv, se.sim.MNNPD.inv,
                   z.MNNPD.inv, MNNPD.inv.rankLo, MNNPD.inv.rankHi, 
                   p.rank.MNNPD.inv, p.rank.MNNPD.inv.oneminusdiff, p.rank.MNNPD.inv.Lo, p.rank.MNNPD.inv.Hi,  one.tailed.p.z.MNNPD.inv,  N)
  
  MPDp.i <- rbind(island,  "MPD_inv_nat", outputSummary["n.invasives",], obs.MPD.i, mean.sim.MPD.inv, median.sim.MPD.inv, sd.sim.MPD.inv, se.sim.MPD.inv, 
                  z.MPD.inv, MPD.inv.rankLo, MPD.inv.rankHi,
                  p.rank.MPD.inv, p.rank.MPD.inv.oneminusdiff, p.rank.MPD.inv.Lo, p.rank.MPD.inv.Hi, one.tailed.p.z.MPD.inv,  N)
  
  output.MNNPD <- cbind(MNNPDp.i, MPDp.i)
  
  rownames(output.MNNPD) <- c("island", "Metric", "n.inv.tax", "mean.obs.metric", "mean.randomized.null", "median.randomized.null", "sd.randomized.null", "se.randomized.null", 
                             "obs.z", "rankLow", "rankHi", 
                             "p.value.ranks","p.value.ranks.oneminusdiff", "p.value.ranks.lo","p.value.ranks.hi", "one.tailed.p.value.z.", "runs")

  colnames(output.MNNPD) <- c(paste(island), paste(island))
  
  return(as.data.frame(output.MNNPD))
  
}

#################################### Simulate MNNFD, MFD  #################################### 
# phy = community phylogeny
# com = observed community matrix
# island = name of community (island); = colnames in com 
# traits = trait dataset; rownames = species names, colnames[,1] = "Status"
# traitname = name of trait to calculate distance of (= colname in traits)
# N = number of communities to simulate 
sim.meanMNNFD.MFD <- function(phy, com, island, traits, traitname, N){
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
    All.Sim.com.Dist <- phyloDistinct(phy=phy, community=comm.island, col=(names(comm.island)))#apply MNNPD funciton across all communities ...find.NN
    All.Sim.com.Functional.Dist <- functionDistinct(output=All.Sim.com.Dist, traits=traits, traitname=traitname)
    if (nrow(All.Sim.com.Functional.Dist) == 0){
      tmp <- c("NA", "NA", "NA", "NA",  "NA", "NA")
    } else {
      All.Sim.com.Funct.summaryMNNFD <- functionObsSum(All.Sim.com.Functional.Dist) #apply summary funciton across all communities...summary.MNNPD
      
      tmp <- c(All.Sim.com.Funct.summaryMNNFD["n.natives",], All.Sim.com.Funct.summaryMNNFD["n.invasives",], 
               All.Sim.com.Funct.summaryMNNFD["meanMNNFDinvasives",], All.Sim.com.Funct.summaryMNNFD["meanMNNFDnatives",], 
               All.Sim.com.Funct.summaryMNNFD["meanMFDinv_nat", ], All.Sim.com.Funct.summaryMNNFD["meanMFDnat_nat", ])
    }
    
    sim.mean <- rbind(sim.mean, tmp)
  }
  colnames(sim.mean) <- c("n.native.tips", "n.invasive.tips", "meanMNNFDinvasives", "meanMNNFDnatives", "meanMFDinv_nat", "meanMFDnat_nat")
  return(sim.mean)
  
}

############ Standardized effect size of Functional Distances ######################
### null distribution ==  randomizeCommunity()
## One islands; one trait

# phy = community phylogeny
# com=SJcommNewSim ## names of communities to simulate
# island = island name (e.g.g "Reeflet_Island")
# simOneIslandOneTrait = output of sim.meanMNNFD.MFD()
# outputMNNPD = phyloObs[["Reeflet_Island"]] 
# traits=SJtraitLog ## trait dataset
# traitname = trait name (e.g. "maxHeight")
# N = number of community randomizations 
ses.FunctionDist <- function(phy, com, island, simOneIslandOneTrait, outputMNNPD, traits, traitname, N){
  
  ## Ouput of phyloDistinct
  output <- functionDistinct(outputMNNPD, traits, traitname)
  
  if (length(na.omit(as.numeric(as.character(output[,2])))) == 1){
    p <- as.data.frame(cbind(island, "has only one species", "NA","NA", "NA","NA","NA", "NA","NA","NA", "NA", "NA","NA", "NA"))
    colnames(p) <- c("island", "Metric", "n.inv.tax", "mean.obs.metric", "mean.randomized.null", "median.randomized.null", "sd.randomized.null", "se.randomized.null", 
                               "obs.z", "rankLow", "rankHi", "p.value.ranks", "one.tailed.p.value.z.", "runs")
    rownames(p) <- island
    
    p <- t(p)
    return(p)
    #return(print("This island has only one species; go fuck yourself"))
  } else if (length(which(output[,"Species.Status"]=="n")) == 0){
    p <- as.data.frame(rbind(island, "has no native species", "NA","NA", "NA","NA","NA", "NA","NA","NA", "NA", "NA", "NA", "NA"))
    rownames(p) <- c("island", "Metric", "n.inv.tax", "mean.obs.metric", "mean.randomized.null", "median.randomized.null", "sd.randomized.null", "se.randomized.null", 
                     "obs.z", "rankLow", "rankHi", "p.value.ranks", "one.tailed.p.value.z.", "runs")
    colnames(p) <- island
    p <- t(p)
    return(p)
    #return(print("This island has no native species; go fuck yourself"))
  } else if (length(which(output[,"Species.Status"]=="i")) == 0){
    p <- as.data.frame(rbind(island, "has no invasive species", "NA","NA", "NA","NA","NA", "NA","NA","NA", "NA", "NA", "NA", "NA"))
    rownames(p) <- c("island", "Metric", "n.inv.tax", "mean.obs.metric", "mean.randomized.null", "median.randomized.null", "sd.randomized.null", "se.randomized.null", 
                     "obs.z", "rankLow", "rankHi", "p.value.ranks", "one.tailed.p.value.z.", "runs")
    colnames(p) <- island
    p <- t(p)
    return(p)
  } else if (length(which(output[,"Species.Status"]=="n")) == 1){
    p <- as.data.frame(cbind(island, "has only one native species", "NA","NA", "NA","NA","NA", "NA","NA","NA", "NA", "NA","NA", "NA"))
    colnames(p) <- c("island", "Metric", "n.inv.tax", "mean.obs.metric", "mean.randomized.null", "median.randomized.null", "sd.randomized.null", "se.randomized.null", 
                     "obs.z", "rankLow", "rankHi", "p.value.ranks", "one.tailed.p.value.z.", "runs")
    rownames(p) <- island
    
    p <- t(p)
    return(p)
  } else {
    # observed mean MNNFD /MFD for natives, invasives
    outputSummary <- functionObsSum(output)
    
    # simulated mean MNNFD / MFD by randomizing invasive species in each community
    simOneIslandOneTrait[island]
    
    ##### ses.MNNFD.invasives
    #mean of randomization 
    meanMNNFDinvasives.inv.dist <- as.numeric(na.omit(simOneIslandOneTrait[island][[1]][[3]])) #meanMNNFDinvasives
    mean.sim.MNNFD.inv <- mean(meanMNNFDinvasives.inv.dist, na.rm=T)
    median.sim.MNNFD.inv <- median(meanMNNFDinvasives.inv.dist, na.rm=T)
    sd.sim.MNNFD.inv <- sd(meanMNNFDinvasives.inv.dist, na.rm=T)
    se.sim.MNNFD.inv <- mean.sim.MNNFD.inv / sqrt(N)
    #observed mean
    obs.MNNFD.i <- as.numeric(as.character(outputSummary["meanMNNFDinvasives",]))
    #Z-value assuming that the null hypothesis is true 
    z.MNNFD.inv <- (obs.MNNFD.i - mean.sim.MNNFD.inv) / sd.sim.MNNFD.inv 
    # quantile / ranks = frequency that oberved pattern is greater than the randomized values
    MNNFD.inv.rankHi <- sum(meanMNNFDinvasives.inv.dist >= as.numeric(obs.MNNFD.i)) #simulated means as or more extremely lower than observed
    MNNFD.inv.rankLo <- sum(meanMNNFDinvasives.inv.dist <= as.numeric(obs.MNNFD.i))  #simulated means as or more extremely lower than observed
    # p-values of observed mean vs. distribution of randomized means 
    p.rank.MNNFD.inv <- min(MNNFD.inv.rankLo, MNNFD.inv.rankHi) / (N + 1) #proportion of simulated means as or more extreme than observed
    
    # p-values of z-value
    one.tailed.p.z.MNNFD.inv <- pnorm(-abs(z.MNNFD.inv))  #p-value of your sample is the lowest alpha level you could have used for your test and still rejected the null hypothesis given your sample. 
    
    ##### ses.MFD.inv.nat
    #mean of randomization 
    meanMFDinvasives.inv.dist <- as.numeric(na.omit(simOneIslandOneTrait[island][[1]][[5]])) #meanMFDinvasives
    mean.sim.MFD.inv <- mean(meanMFDinvasives.inv.dist, na.rm=T)
    median.sim.MFD.inv <- median(meanMFDinvasives.inv.dist, na.rm=T)
    sd.sim.MFD.inv <- sd(meanMFDinvasives.inv.dist, na.rm=T)
    se.sim.MFD.inv <- mean.sim.MFD.inv / sqrt(N)
    #observed mean
    obs.MFD.i <- as.numeric(as.character(outputSummary["meanMFDinv_nat",]))
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
    MNNFDp.i <- rbind(island, "MNNFD_inv", outputSummary["n.invasives",], obs.MNNFD.i, mean.sim.MNNFD.inv, median.sim.MNNFD.inv, sd.sim.MNNFD.inv, se.sim.MNNFD.inv,
                     z.MNNFD.inv, MNNFD.inv.rankLo, MNNFD.inv.rankHi, p.rank.MNNFD.inv, one.tailed.p.z.MNNFD.inv,  N)
    
    MFDp.i <- rbind(island,  "MFD_inv_nat", outputSummary["n.invasives",], obs.MFD.i, mean.sim.MFD.inv, median.sim.MFD.inv, sd.sim.MFD.inv, se.sim.MFD.inv, 
                    z.MFD.inv, MFD.inv.rankLo, MFD.inv.rankHi,p.rank.MFD.inv, one.tailed.p.z.MFD.inv,  N)
    
    output.MNNFD <- cbind(MNNFDp.i, MFDp.i)
    rownames(output.MNNFD) <- c("island", "Metric", "n.inv.tax", "mean.obs.metric", "mean.randomized.null", "median.randomized.null", "sd.randomized.null", "se.randomized.null", 
                               "obs.z", "rankLow", "rankHi", "p.value.ranks", "one.tailed.p.value.z.", "runs")
    colnames(output.MNNFD) <- c(paste(island), paste(island))
    
    return(as.data.frame(output.MNNFD))
    
  }
  
  
}










