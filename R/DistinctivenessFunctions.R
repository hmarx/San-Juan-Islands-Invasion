###### Distinctiveness Functions #######
###### 6 Nov 2014 ########
###### Hannah E. Marx #######

### For analysis of phylogenetic and functional distinctivness to nearest native (DNNS/NNFD) 
###### and mean native community (MDNS/MFD)
### Also, standardized effect size compared to null model randomizing invasive species occurences
### Used for Marx et al. 2015

## Data (analysis.R)
## Phylogeny: phylo object with all species in "source pool"
## Community data: data.frame with species names matching phylogeny, presence/absence of speices in each community: 
# row = species names; columns = 1 (presnce), or 0 (absence)
## Trait data: data.frame with species names mathcing phylogeny, trait values 
# row = species names (must match format of names in phylogeny); columns = trait value

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


#### Summarize phyloDistinct output
# output = output from phyloDistinct()
summary.DNNS.MDNS<- function(output){
  summ <- data.frame()
  if (nrow(output)==1){
    summ <- rbind("NA", "NA", "NA", "NA", "NA","NA", "NA", "NA", "NA","NA", "NA", "NA")
    rownames(summ) <- c("meanDNNSnatives", "meanDNNSinvasives", "n.natives", "n.invasives", "t.DNNS.p.value", "t.DNNS.conf.int.Lo","t.DNNS.conf.int.Hi", 
                        "meanMDNSnat_nat", "meanMDNSinv_nat", "t.MDNS.p.value", "t.MDNS.conf.int.Lo","t.MDNS.conf.int.Hi")
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
                      "meanMDNSnat_nat", "meanMDNSinv_nat", "t.MDNS.p.value", "t.MDNS.conf.int.Lo","t.MDNS.conf.int.Hi")
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

############################  Differnece in trait values between 1) each species and nearest Neighbor (NNFD)
############################                                     2) each species and the native community (MFD)

#### For one trait, calculates Nearest Neighbor Trait Differences (NNFD), 
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
    colnames(diff.data) <-c("Species.Status", paste("NNFD", traitname, sep="_"), "MFD.n",  "MFD.i")
    return(diff.data)#(print("has only one species; go fuck yourself"))
  }
  if (nrow(as.data.frame(output[output["Species.Status"]=="n",]))==1){
    dat <- cbind("NA", "NA", "NA", "NA")
    diff.data <- rbind(diff.data, dat)
    colnames(diff.data) <- c("Species.Status", paste("NNFD", traitname, sep="_"), "MFD.n",  "MFD.i")
    return(diff.data)#(print("has only one native species; go fuck yourself"))
  }
  if (nrow(as.data.frame(output[output["Species.Status"]=="i",]))==0){
    dat <- cbind("NA", "NA", "NA", "NA")
    diff.data <- rbind(diff.data, dat)
    colnames(diff.data) <- c("Species.Status", paste("NNFD", traitname, sep="_"), "MFD.n",  "MFD.i")
    return(diff.data)#(print("has only one native species; go fuck yourself"))
  }
  output <- as.data.frame(output) ### DNNS results; have info about nearest native species
  
  # merge Communty and trait data for one trait: Row.names, Species.Status,  Nearest.native, leafletSize
  comDNNSTrait <- merge(output[c("Species.Status","Nearest.native")], na.omit(traits[traitname]), by=0) #merge DNNS data and trait data; remove NA
  if (nrow(comDNNSTrait) == 0) {
    dat <- cbind("NA", "NA", "NA", "NA")
    diff.data <- rbind(diff.data, dat)
    colnames(diff.data) <- c("Species.Status", paste("NNFD", traitname, sep="_"), "MFD.n",  "MFD.i")
    return(diff.data)#(print("has only one native species; go fuck yourself"))
  }
  
  dim(comDNNSTrait)
  i.comDNNSTrait <- (comDNNSTrait)[which(comDNNSTrait[2] == "i"),] ## just invasive species with trait on island
  n.comDNNSTrait <- (comDNNSTrait)[which(comDNNSTrait[2] == "n"),] ## native species with trait on island
  
  # create distance matrix of trait values == absolute value of difference in trait values between each species 
  trait.dist <- as.matrix(dist(comDNNSTrait[4], upper=T))
  rownames(trait.dist) <- comDNNSTrait$Row.names
  colnames(trait.dist) <- comDNNSTrait$Row.names
  dim(trait.dist)
  head(trait.dist)
  
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
functionObsSum <- function(outputTraitDistance){
  outputTraitDistance.naomit <- (outputTraitDistance)[which(!outputTraitDistance[2] == "NA"), ]
  summ <- data.frame()
  if (nrow(outputTraitDistance.naomit)==0){
    summ <- rbind("NA", "NA", "NA", "NA", "NA","NA", "NA", "NA", "NA","NA", "NA", "NA")
    rownames(summ) <- c("meanNNFDnatives", "meanNNFDinvasives", "n.natives", "n.invasives", "t.NNFD.p.value", "t.NNFD.conf.int.Lo","t.NNFD.conf.int.Hi", 
                        "meanMFDnat_nat", "meanMFDinv_nat", "t.MFD.p.value", "t.MFD.conf.int.Lo","t.MFD.conf.int.Hi")  #colnames(summ) <- names(traits[col])
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
                      "meanMFDnat_nat", "meanMFDinv_nat", "t.MFD.p.value", "t.MFD.conf.int.Lo","t.MFD.conf.int.Hi")  #colnames(summ) <- names(traits[col])
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
    
    tmp <- c(All.Sim.com.DistsummaryDNNS["n.natives",], All.Sim.com.DistsummaryDNNS["n.invasives",], All.Sim.com.DistsummaryDNNS["meanDNNSinvasives",], All.Sim.com.DistsummaryDNNS["meanDNNSnatives",], All.Sim.com.DistsummaryDNNS["meanMDNSinv_nat", ], All.Sim.com.DistsummaryDNNS["meanMDNSnat_nat", ])
    sim.mean <- rbind(sim.mean, tmp)
  }
  colnames(sim.mean) <- c("n.native.tips", "n.invasive.tips", "meanDNNSinvasives", "meanDNNSnatives", "meanMDNSinv_nat", "meanMDNSnat_nat")
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
  obs.MDNS.i <- as.numeric(as.character(outputSummary["meanMDNSinv_nat",]))
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
      All.Sim.com.Funct.summaryNNFD <- functionObsSum(All.Sim.com.Functional.Dist) #apply summary funciton across all communities...summary.DNNS
      
      tmp <- c(All.Sim.com.Funct.summaryNNFD["n.natives",], All.Sim.com.Funct.summaryNNFD["n.invasives",], 
               All.Sim.com.Funct.summaryNNFD["meanNNFDinvasives",], All.Sim.com.Funct.summaryNNFD["meanNNFDnatives",], 
               All.Sim.com.Funct.summaryNNFD["meanMFDinv_nat", ], All.Sim.com.Funct.summaryNNFD["meanMFDnat_nat", ])
    }
    
    sim.mean <- rbind(sim.mean, tmp)
  }
  colnames(sim.mean) <- c("n.native.tips", "n.invasive.tips", "meanNNFDinvasives", "meanNNFDnatives", "meanMFDinv_nat", "meanMFDnat_nat")
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
    # observed mean NNFD /MFD for natives, invasives
    outputSummary <- functionObsSum(output)
    
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
      colnames(newlist) <- c("n.native.tips", "n.invasive.tips", "meanNNFDinvasives", "meanNNFDnatives", "meanMFDinv_nat", "meanMFDnat_nat")
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
sum.sesFunctionDist <- function(plottype=c("NullInvOcc", "ses.allIslands", "summary.Bar"), sim.output, islands.sim, phyloObs, 
                                traits, traitname, metadata){
  ## Observed Values
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
  sim.null.distrib.melt <- sim.null.distrib.melt[which(!sim.null.distrib.melt$value == "NA"),]
  #head(sim.null.distrib.melt)
  if(plottype[1] == "NullInvOcc"){
    pdf(file=paste("figs/plots/functionDiv/ses/NullInvOcc.NNFD", traitname, "pdf", sep="."), width=20, height=10)
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
    print(p1)
    dev.off()
    
    pdf(paste("figs/plots/functionDiv/ses/NullInvOcc.MFD", traitname, "pdf", sep="."), width=20, height=10)
    p2 <- ggplot(sim.null.distrib.melt, aes(x=reorder((L1), Area.m2), y=as.numeric(as.character(meanMFDinv_nat))), 
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
    print(p2)
    dev.off()
    
  }
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
  #head(ses.SJ.NNFD)
  #dim(ses.SJ.NNFD) #144
  ses.SJ.NNFD <- na.omit(ses.SJ.NNFD)
  ses.SJ.NNFD$p.value.ranks <- as.numeric(as.character(ses.SJ.NNFD$p.value.ranks)) 
  ses.SJ.NNFD.sig <- subset(ses.SJ.NNFD, p.value.ranks <= 0.05)
  sig = (ses.SJ.NNFD[,"p.value.ranks"] <= 0.05)
  ses.SJ.NNFD <- cbind(ses.SJ.NNFD, sig)
  
  ses.SJ.NNFD[,"Significance"] <- ifelse(ses.SJ.NNFD[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.SJ.NNFD[,"obs.z"])) <= 0, 1,
                                          ifelse(ses.SJ.NNFD[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.SJ.NNFD[,"obs.z"])) >= 0, 2, 
                                                 ifelse(ses.SJ.NNFD[,"p.value.ranks"] >= 0.05 & as.numeric(as.character(ses.SJ.NNFD[,"obs.z"])) <= 0, 3, 4)))
  
  if(plottype[1] == "ses.allIslands"){
    x1 <- length(which(ses.SJ.NNFD$Size.cat == "sm"))/2 + .5
    x2 <- x1 + length(which(ses.SJ.NNFD$Size.cat == "med"))/2
    pdf(paste("figs/plots/functionDiv/ses/ses.SJ.NNFD", traitname, "pdf", sep="."), width=20, height=10)
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
    print(p3)
    dev.off()
    
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
    
    pdf(paste("figs/plots/functionDiv/ses/", traitname, "Functional.SummaryBar.pdf", sep=""))
    p4 <- ggplot(sesFunctional2, aes(x=Significance, fill=factor(Size.cat))) + 
      geom_bar(stat="bin") +
      scale_fill_manual(name="Size Category",
                        breaks=c("sm", "med", "lg"),
                        values=c("sm"="dodgerblue4", "med"="orangered3", "lg"="gold1")) +
      theme_bw() +
      scale_x_discrete(name="", labels=c("0"="NS", "1"="Significant")) +
      facet_wrap(~Metric, ncol = 2) +
      ggtitle(paste(traitname, "SES Functional Difference\n",sep=" "))
    print(p4)
    dev.off()
  }
  
  
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
  
  nrow(sigNNFD_inv)
  nrow(sigMFD_inv)
  
  nrow(sigNNFD_inv) / nrow(MFD_inv_nat)
  nrow(sigMFD_inv) / nrow(MFD_inv_nat)
  
  sum <- rbind(nrow(sigNNFD_inv), nrow(sigNNFD_inv) / nrow(MFD_inv_nat), sm.NNFD, med.NNFD, lg.NNFD,
               nrow(sigMFD_inv),  nrow(sigMFD_inv) / nrow(MFD_inv_nat), sm.MFD, med.MFD, lg.MFD,
               nrow(MFD_inv_nat))
  rownames(sum) <- c("n sig NNFD", "% sig NNFD", "n sig small NNFD", "n sig medium NNFD", "n sig large NNFD",
                     "n sig MFD", "% sig MFD", "n sig small MFD", "n sig medium MFD", "n sig large MFD", 
                     "Total Islands")
  colnames(sum) <- traitname
  
  return(sum)
  
}








