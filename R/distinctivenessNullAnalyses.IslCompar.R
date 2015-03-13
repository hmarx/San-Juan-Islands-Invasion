
#### Functions specific to Island comparative study

randomizeCommunityDF <- function(com, traits, N){
  community.matrix.list <- list()
  com.tmp <- as.list(com[[1]]) 
  
  inv.com <- rownames(traits[traits$Status=="i",]) # natives in species pool
  
  nat.com <- rownames(traits[traits$Status=="n",]) # natives in species pool
  com.island <- names(com.tmp[com.tmp == "n" ])
  nat.com.new <- (c(com.island, nat.com))
  
  rownames(com.island[com.island$Status == "i",]) %in% inv.com
  
  length(com.tmp)
  names(com.tmp) <- rownames(com)
  exclude.natives <- com.tmp[which(names(com.tmp) %in% nat.com.new)] # taxa that are in larger native community: don't shuffle
  shuffle.tips <- com.tmp[which(!names(com.tmp) %in% nat.com.new)] # get all the taxa which are invasives in larger pool
  ###Simulation 
  rep <- 0
  while (rep < N){
    rep <- rep + 1
    sim.comm <- data.frame(row.names=rownames(com), stringsAsFactors=FALSE)
    u <- unlist(shuffle.tips)
    l2 <- relist(u[sample(length(u))],skeleton=shuffle.tips)
    #tmp <- c(l2, exclude.natives)
    #df <- data.frame(matrix(unlist(tmp), byrow=T), row.names=names(tmp))
    new.com <- merge(com[1], t(as.data.frame(l2, stringsAsFactors=FALSE)), by=0, all.x=T)
    new.com[new.com[ , 2] == "n",][3] <- "n"
    new.com[is.na(new.com)] <- 0
    #new.com.final <- l3[,3, drop=F]
    new.com.final <- as.data.frame(new.com[, 3, drop=F], row.names=new.com[,1])
    #new.com.final <- as.data.frame(cbind(new.com[, 3], rownames= new.com[,1], drop=F))
    sim.comm <- cbind(sim.comm, new.com.final)
    colnames(sim.comm) = colnames(com[1:length(com)]) 
    print(length(sim.comm[sim.comm=="i"]))
    community.matrix.list[[rep]] <- sim.comm
  }
  return(community.matrix.list)
}


sim.meanMNNPD.MPD.IslComp <- function(phy, com, island, traits, N){
  tmp <- treedata(phy, com) 
  new.phy <- tmp$phy ## prune the tree to include just species in larger species pool to draw from
  com.data <- as.data.frame(tmp$data, stringsAsFactors=F)   ## subset the community data, to make sure in the right order
  #dim(com.data)
  #head(com.data)
  ## Create a distribution of ranomized communities, summarize 
  rand.inv.com <- randomizeCommunityDF(com=com.data[island], traits=traits, N=N)
  
  sim.mean <- data.frame()
  
  for (i in 1:length(rand.inv.com)){
    comm.island <- data.frame(lapply(rand.inv.com[[i]], as.character), stringsAsFactors=F)
    #length(which(comm.island[1] == "n"))
    rownames(comm.island) <- rownames(rand.inv.com[[i]])
    All.Sim.com.Dist <- phyloDistinct(phy=new.phy, community=comm.island, col=(names(comm.island)))#apply MNNPD funciton across all communities ...find.NN
    All.Sim.com.DistsummaryMNNPD <- summary.MNNPD.MPD(All.Sim.com.Dist) #apply summary funciton across all communities...summary.MNNPD
    
    tmp <- c(All.Sim.com.DistsummaryMNNPD["n.natives",], All.Sim.com.DistsummaryMNNPD["n.invasives",], All.Sim.com.DistsummaryMNNPD["meanMNNPDinvasives",], All.Sim.com.DistsummaryMNNPD["meanMNNPDnatives",], All.Sim.com.DistsummaryMNNPD["meanMPDinv_nat", ], All.Sim.com.DistsummaryMNNPD["meanMPDnat_nat", ])
    sim.mean <- rbind(sim.mean, tmp)
  }
  colnames(sim.mean) <- c("n.native.tips", "n.invasive.tips", "meanMNNPDinvasives", "meanMNNPDnatives", "meanMPDinv_nat", "meanMPDnat_nat")
  return(sim.mean)
  
}




