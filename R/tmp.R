

remove.islands.sim <- c("All_SanJuanIslands", "Unnamed_west_of_Castle_Island") 
SJcommNewSim  <- SJcommNew[, -which(names(SJcommNew) %in% remove.islands.sim)]
head(SJcommNewSim)
SJ_islands.sim <- colnames(SJcommNewSim)

## Null distribution...this will take a while
sim.null.distrib.Height <- lapply(SJ_islands.sim, function(x) sim.meanMNNFD.MFD(phy=SJfinalTree, com=SJcommNewSim, traits=SJtraitLog, island=x, traitname="maxHeight", N = 1000))
names(sim.null.distrib.Height) <- SJ_islands.sim
head(sim.null.distrib.Height)
head(sim.null.distrib.Height["Willow_Island"])
#write.csv(sim.null.distrib.Height, file="sim.null.distrib.Height.null1000.csv")

# output = .csv of sim.meanMNNFD.MFD()
# names of islands 
read.nullOutput <- function(sim.output, islands.sim){
  ## Load output of sim.meanMNNFD.MFD()
  sdf <- read.csv(file=sim.output, as.is=T, row.names=1)
  head(sdf[, 1:6])
  list.sdf <- list()
  try(
  for (i in 1:ncol(sdf)){
    newlist <- (cbind(sdf[,c(1:6)]))
    colnames(newlist) <- c("n.native.tips", "n.invasive.tips", "meanMNNFDinvasives", "meanMNNFDnatives", "meanMFDinv_nat", "meanMFDnat_nat")
    list.sdf[[i]] <- newlist
    sdf <- sdf[, -c(1:6)]
    
  },   silent=TRUE)
  names(list.sdf) <- islands.sim
  return(list.sdf)
  
}

tmp <- read.nullOutput(sim.output  = "output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.Height.null1000.csv", islands.sim = SJ_islands.sim)


sum.sesFunctionDist(plottype = "ses.allIslands", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.Height.null1000.csv", 
               islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="maxHeight", metadata=metadata)

sum.sesFunctionDist(plottype = "ses.allIslands", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.SeedMass.1000.CSV", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="seedMass", metadata=metadata)

sum.sesFunctionDist(plottype = "ses.allIslands", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.SLA.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="sla", metadata=metadata)


# plottype = "NullInvOcc" (boxplot showing distribution of null means, with observed mean) 
#             OR "ses.allIslands" (points of z-scores for SES of each observed mean to null for each island, significant islands hightlighted)
# sim.output = directory of .csv file for the output of the null distribution 
# islands.sim = list of the island names to use
# phyloObs = output of phyloDistinct() applies across all communities 
# traits = trait file; rownames = species, column[1] = "Status", other columns = trait values
# traitname = name of the trait to analize
# metadata = rownames = island names; one column == "Area.m2
sum.sesFunctionDist <- function(plottype=c("NullInvOcc", "ses.allIslands"), sim.output, islands.sim, phyloObs, 
                         traits, traitname, metadata){
  ## Observed Values
  obs.MNNFD <- lapply(islands.sim, function(x) functionDistinct(output=phyloObs[[x]], traits, traitname)) 
  names(obs.MNNFD) <- islands.sim 
  ## Summary Observed Values
  summ.MNNFD  <- lapply(islands.sim, function(x) functionObsSum(obs.MNNFD[[x]])) 
  names(summ.MNNFD) <- islands.sim 
  
  ## Read in output of means from simulated communites  (list elements = communities )
  simIslands <- read.nullOutput(sim.output, islands.sim)
  
  ## append a column of metadata to each element in list
  list.meta.null.distrib <- list()
  for (i in 1:length(simIslands)){ 
    tmp <- metadata[as.character(names(simIslands[i])), "Area.m2"]
    tmp.MNNFD.i <- rep(summ.MNNFD[as.character(names(summ.MNNFD[i]))][[1]][[2]], times=nrow(simIslands[[i]]))
    tmp.MFD.in <- rep(summ.MNNFD[as.character(names(summ.MNNFD[i]))][[1]][[9]], times=nrow(simIslands[[i]]))
    newlist <- mapply(cbind, "islands"=names(simIslands[i]), simIslands[i], "Area.m2"=tmp, "Obs.meanMNNFDinvasives"=tmp.MNNFD.i, "Obs.meanMPFDinv_nat"= tmp.MFD.in, SIMPLIFY=F) 
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
    pdf(file=paste("figs/plots/functionDiv/ses/NullInvOcc.MNNFD", traitname, "pdf", sep="."), width=20, height=10)
        p1 <- ggplot(sim.null.distrib.melt, aes(x=reorder((L1), Area.m2), y=as.numeric(as.character(meanMNNFDinvasives))), 
                    position=position_dodge(width=1)) +
        geom_boxplot(aes(fill=factor(as.character(variable))), width = 1) +
        scale_x_discrete("Island (increasing size)") +
        scale_y_discrete("Null distribution mean(MNNFD)") +
        theme_bw() +
        scale_fill_manual(values="grey", labels="") + 
        guides(fill=guide_legend(title=" Null distibuiton : random invasive occurrence")) +
        theme(legend.position="top") +
        geom_point(aes(y = Obs.meanMNNFDinvasives), shape=1, color="magenta1") +
        theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
        ggtitle(paste("Null distribution and observed MNNFD of\n", traitname, sep=" ")) + 
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
  ses.SanJuan.MNNFD.MFD <- lapply(islands.sim, function(x) ses.FunctionDist(phy=SJfinalTree, com=SJcommNewSim, island=x, 
                                                                                      simOneIslandOneTrait=simIslands, outputMNNPD=phyloObs[[x]], traits=SJtraitLog, traitname=traitname, N=1000))
  names(ses.SanJuan.MNNFD.MFD) <- islands.sim
  
  listIslands <- ses.SanJuan.MNNFD.MFD
  ses.SJ.MNNFD <- data.frame()
  for (i in 1:length(listIslands)){ 
    if (length(listIslands[[i]]) != 2){
      newlist <- NULL
    } else {
      tmp1 <- metadata[as.character(names(listIslands[i])), "Area.m2"]
      tmp2 <- metadata[as.character(names(listIslands[i])), "Size.cat"]
      newlist <-cbind(t(listIslands[[i]]), "Area.m2"=tmp1, "Size.cat"=tmp2) 
    }
    ses.SJ.MNNFD <- rbind(ses.SJ.MNNFD, newlist)
    
  }
  #head(ses.SJ.MNNFD)
  #dim(ses.SJ.MNNFD) #144
  ses.SJ.MNNFD <- na.omit(ses.SJ.MNNFD)
  ses.SJ.MNNFD$p.value.ranks <- as.numeric(as.character(ses.SJ.MNNFD$p.value.ranks)) 
  ses.SJ.MNNFD.sig <- subset(ses.SJ.MNNFD, p.value.ranks <= 0.05)
  sig = (ses.SJ.MNNFD[,"p.value.ranks"] <= 0.05)
  ses.SJ.MNNFD <- cbind(ses.SJ.MNNFD, sig)
  
  ses.SJ.MNNFD[,"Significance"] <- ifelse(ses.SJ.MNNFD[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.SJ.MNNFD[,"obs.z"])) <= 0, 1,
                                                 ifelse(ses.SJ.MNNFD[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.SJ.MNNFD[,"obs.z"])) >= 0, 2, 
                                                        ifelse(ses.SJ.MNNFD[,"p.value.ranks"] >= 0.05 & as.numeric(as.character(ses.SJ.MNNFD[,"obs.z"])) <= 0, 3, 4)))

  if(plottype[1] == "ses.allIslands"){
    pdf(paste("figs/plots/functionDiv/ses/ses.SJ.MNNFD", traitname, "pdf", sep="."), width=20, height=10)
    p3 <- ggplot(ses.SJ.MNNFD, aes(x=reorder(factor(island),as.numeric(as.character(Area.m2))), 
                                    y=as.numeric(as.character((obs.z))), color=factor(sig), shape=Metric, fill=factor(sig))) +
      geom_point(size=10) +  
      coord_cartesian(ylim=c(-5, 5)) + 
      scale_y_continuous("standardized effect size", breaks=seq(-8, 8, 2)) +
      scale_x_discrete("Island (increasing size)") +
      scale_fill_manual(values=alpha(c("FALSE"="white", "TRUE"= "black"), .3), guide="none") +
      scale_color_manual(name =" ",values=c("FALSE"= "grey", "TRUE"="black"), guide="none") +
      scale_shape_manual(name =" ",values=c("MNNFD_inv"= 21, "MFD_inv_nat"= 24), labels=c("MNNFD_inv"="MNNFD i  ", "MFD_inv_nat"="MFD_inv_nat")) +
      theme_bw() +
      theme(legend.position="top") +
      geom_abline(intercept = 0, slope = 0, colour = "grey", size = .5) +
      geom_vline(xintercept = c(18.5, 53.5)) +
      theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
      ggtitle(paste("Significane of", traitname, "difference for\n invasive species to nearest native (NNFDi),\n and native community (MFDNi)\n", sep=" ")) +
      theme(plot.title=element_text(size=rel(1.5))) +
      theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
    print(p3)
    dev.off()
    
  }
  
  #### summarize significance data by size categories
  ses.SJ.MNNFD$Size.cat <- as.character(ses.SJ.MNNFD$Size.cat)
  ses.SJ.MNNFD$Size.cat <- factor(ses.SJ.MNNFD$Size.cat, levels=c("sm", "med", "lg"))
  
  ##### ses MMNPDi 
  MNNFD_inv <- na.omit(ses.SJ.MNNFD[which(ses.SJ.MNNFD$Metric == "MNNFD_inv"), ] )
  sigMNNFD_inv <- (MNNFD_inv[which(as.numeric(as.character(MNNFD_inv$p.value.ranks)) <= 0.05), ]) #7
  ## Significant by size category
  sm.MNNFD <- length(which(sigMNNFD_inv$Size.cat == "sm")) #3/71  0.04225352
  med.MNNFD <- length(which(sigMNNFD_inv$Size.cat == "med")) #3/71   0.04225352
  lg.MNNFD <- length(which(sigMNNFD_inv$Size.cat == "lg")) #1/71   0.01408451
  
  
  #### ses MPDin 
  MFD_inv_nat <- ses.SJ.MNNFD[which(ses.SJ.MNNFD$Metric == "MFD_inv_nat"), ] 
  sigMFD_inv <- (MFD_inv_nat[which(as.numeric(as.character(MFD_inv_nat$p.value.ranks)) <= 0.05), ]) #13
  ## Significant by size category
  sm.MFD <- length(which(sigMFD_inv$Size.cat == "sm")) #2/71  0.02816901
  med.MFD <- length(which(sigMFD_inv$Size.cat == "med")) #7/71   0.09859155
  lg.MFD <- length(which(sigMFD_inv$Size.cat == "lg")) #4/71   0.05633803
  
  nrow(sigMNNFD_inv)
  nrow(sigMFD_inv)
  
  nrow(sigMNNFD_inv) / nrow(MFD_inv_nat)
  nrow(sigMFD_inv) / nrow(MFD_inv_nat)
  
  sum <- rbind(nrow(sigMNNFD_inv), nrow(sigMNNFD_inv) / nrow(MFD_inv_nat), sm.MNNFD, med.MNNFD, lg.MNNFD,
               nrow(sigMFD_inv),  nrow(sigMFD_inv) / nrow(MFD_inv_nat), sm.MFD, med.MFD, lg.MFD,
                nrow(MFD_inv_nat))
  rownames(sum) <- c("n sig MNNFD", "% sig MNNFD", "n sig small MNNFD", "n sig medium MNNFD", "n sig large MNNFD",
                     "n sig MFD", "% sig MFD", "n sig small MFD", "n sig medium MFD", "n sig large MFD", 
                     "Total Islands")
  colnames(sum) <- traitname
  
  return(sum)
 
}








