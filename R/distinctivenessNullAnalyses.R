################### Final Analysis of San Juans Dataset

## Significance of phylogenetic and functional distinctiveness between invasive and native species within each island community 
# compared to a null expectation, randomizing the invasive communtiy in each community with other invasives from the larger species pool

###### 6 Nov 2014  ##########
###### Hannah E. Marx #######


################################################# Read in Final Datasets #####################################################@
source("analysis.R")

SJfinalTree
head(SJtraitLog)
head(SJcommNew)
head(metadata)
dim(metadata)


###########################################  Significance of Phylogenetic Distances  #############################################################

################# RANDOMIZATION == shuffle prese/abs of invasive within each community (ses.MNNPD.MPD.R) ##############
################# Compare observed MMNPDi / MPDinv-nat to random distribution of invasive spevies in each community
## mean obs MMNPDi < mean null MMNPDi == intorduced species is closer to native relative than a random invasive from species pool
## mean obs MMNPDi = mean null MMNPDi == any introduced species would show the same pattern
## mean obs MMNPDi > mean null MMNPDi == intorduced species is more distant to native relative than random; DNH

################################ MNNPD i / MPD inv-nat ################################ 
######################## Simulate means, compare null distribution to the observed mean 

## Simulate MNNPD/MPD by randomizing invasive communities...NOTE: will take a while to run
SJ_islands <- names(SJcommNew) # names of the islands 
sim.null.distrib.SJ <- lapply(SJ_islands, function(x) sim.meanMNNPD.MPD(phy=SJfinalTree, com=SJcommNew, traits=SJtraitLog, island=x, N = 1000))
names(sim.null.distrib.SJ) <- SJ_islands
head(sim.null.distrib.SJ)
head(sim.null.distrib.SJ["All_SanJuanIslands"])
#write.csv(sim.null.distrib.SJ, file="SanJuan.MNNPD.MPD.null1000.csv")

## Read in simulation file
sdf <- read.csv(file="output/10_Analyses/PhylogeneticDiversity/Null/SanJuan.MNNPD.MPD.null1000.csv", as.is=T, row.names=1)
head(sdf[, 1:6])
list.sdf <- list()
for (i in 1:ncol(sdf)){
  newlist <- (cbind(sdf[, c(1:6)]))
  head(newlist)
  colnames(newlist) <- c("n.native.tips", "n.invasive.tips", "meanMNNPDinvasives", "meanMNNPDnatives", "meanMPDinv_nat", "meanMPDnat_nat")
  list.sdf[[i]] <- newlist
  sdf <- sdf[, -c(1:6)]
  
}
names(list.sdf) <- SJ_islands

## Observed difference in mean MNNPD i / MPD inv-nat
obs.MNNPD.SJ  <- lapply(SJ_islands, function(x) phyloDistinct(phy=SJfinalTree, community=SJcommNew, col=x))
names(obs.MNNPD.SJ) <- SJ_islands # = phyloObs
head(obs.MNNPD.SJ)
summ.MNNPD.SJ  <- lapply(SJ_islands, function(x) summary.MNNPD.MPD(obs.MNNPD.SJ[[x]])) #apply summary funciton across all communities
names(summ.MNNPD.SJ) <- SJ_islands

## remove islands that don't make sense for the null hyp.
remove.islands.sim <- c("All_SanJuanIslands", "Unnamed_west_of_Castle_Island") ## no null species pool, only one invasive species
list.sdf.new  <- list.sdf[-which(names(list.sdf) %in% remove.islands.sim)]
length(list.sdf.new) # 72

## Append metadata, observed means to each element in list of communities 
head(metadata)
listIslands <- list.sdf.new
list.meta.null.distrib.SJ <- list()
for (i in 1:length(listIslands)){ 
  tmp <- metadata[as.character(names(listIslands[i])), "Area.m2"]
  tmp.MNNPD.i <- rep(summ.MNNPD.SJ[as.character(names(summ.MNNPD.SJ[i]))][[1]][[2]], times=nrow(listIslands[[i]])) # observed MNNPD.i
  tmp.MPD.in <- rep(summ.MNNPD.SJ[as.character(names(summ.MNNPD.SJ[i]))][[1]][[9]], times=nrow(listIslands[[i]])) # observed MPD.i to mean native community
  newlist <- mapply(cbind, "islands"=names(listIslands[i]), listIslands[i], "Area.m2"=tmp, "Obs.meanMNNPDinvasives"=tmp.MNNPD.i, "Obs.meanMPDinv_nat"= tmp.MPD.in, SIMPLIFY=F)   
  list.meta.null.distrib.SJ[i] <- newlist[1]
}
names(list.meta.null.distrib.SJ) <- names(list.sdf.new)
head(list.meta.null.distrib.SJ[[1]])

## Melt list of simulated mean MNNPD/MPD and observed means into dataframe, and remove NA
sim.null.distrib.melt <- melt(list.meta.null.distrib.SJ, measure.vars="islands")
sim.null.distrib.melt <- sim.null.distrib.melt[which(!sim.null.distrib.melt$value == "NA"),] 
head(sim.null.distrib.melt)

## Plot distribution of null means, with observed mean, for each island 
pdf("figs/plots/phyloDiv/ses/logNullInvOcc.MNNPD.islIncSize.pdf", width=20, height=10)
p <- ggplot(sim.null.distrib.melt, aes(x=reorder((L1), Area.m2), y=as.numeric(as.character(meanMNNPDinvasives))), 
            position=position_dodge(width=1)) 
p <- p + geom_boxplot(aes(fill=factor(as.character(variable))), width = 1)
p <- p + scale_x_discrete("Island (increasing size)") 
p <- p + scale_y_log10("log Null distribution mean(MNNPD)")
p <- p + theme_bw() 
p <- p + scale_fill_manual(values="grey", labels="")
p <- p + guides(fill=guide_legend(title="Null distibuiton : random invasive occurrence"))
p <- p + theme(legend.position="top")
p <- p + geom_point(aes(y = Obs.meanMNNPDinvasives), shape=1, color="magenta1") 
p <- p + theme(axis.text.x = element_text(angle = -45, hjust = 0))
p <- p + ggtitle("Null Distribution and observed mean(MNNPD)\n Increasing Island Size") + theme(plot.title=element_text(size=rel(1.5)))
p <- p + theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
p
dev.off()

pdf("figs/plots/phyloDiv/ses/logNullInvOcc.MPD.islIncSize.pdf", width=20, height=10)
p <- ggplot(sim.null.distrib.melt, aes(x=reorder((L1), Area.m2), y=as.numeric(as.character(meanMPDinv_nat))), 
            position=position_dodge(width=1)) 
p <- p + geom_boxplot(aes(fill=factor(as.character(variable))), width = 1)
p <- p + scale_x_discrete("Island (increasing size)") #, breaks=seq(0, 80, 10)) 
p <- p + scale_y_log10("log Null distribution mean(MPD inv-nat)")
p <- p + theme_bw() 
p <- p + scale_fill_manual(values="grey", labels="")
p <- p + guides(fill=guide_legend(title=" Null distibuiton : random invasive occurrence"))
p <- p + theme(legend.position="top")
p <- p + geom_point(aes(y = Obs.meanMPDinv_nat), shape=1, color="magenta1") 
p <- p + theme(axis.text.x = element_text(angle = -45, hjust = 0))
p <- p + ggtitle("Null Distribution and observed mean(MPD inv-nat)\n Increasing Island Size") + theme(plot.title=element_text(size=rel(1.5)))
p <- p + theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
p
dev.off()


#### Summarize observed mean compared to simualted means == standardized effect size 
## p.rank.DNNS.inv <- min(DNNS.inv.rankLo, DNNS.inv.rankHi) / (N + 1) 
## proportion of simulated means as or more extreme than observed
ses.SanJuan.MNNPD.MPD <- lapply(names(list.sdf.new), function(x) ses.PhyloDist(phy=SJfinalTree, com=SJcommNew, island=x, simOneIsland=list.sdf, N=1000))
names(ses.SanJuan.MNNPD.MPD) <- names(list.sdf.new)
length(ses.SanJuan.MNNPD.MPD)

head(metadata)
listIslands <- ses.SanJuan.MNNPD.MPD
ses.MNNPD.MPD <- data.frame()
for (i in 1:length(listIslands)){ 
  tmp1 <- metadata[as.character(names(listIslands[i])), "Area.m2"]
  tmp2 <- metadata[as.character(names(listIslands[i])), "Size.cat"]
  newlist <-cbind(t(listIslands[[i]]), "Area.m2"=tmp1, "Size.cat"=tmp2) 
  ses.MNNPD.MPD <- rbind(ses.MNNPD.MPD, newlist)
}
head(ses.MNNPD.MPD)
#write.csv(ses.MNNPD.MPD, file="ses.MNNPD.MPD.MPDNew2.csv")
ses.MNNPD.MPD$p.value.ranks <- as.numeric(as.character(ses.MNNPD.MPD$p.value.ranks)) 
ses.MNNPD.MPD <- na.omit(ses.MNNPD.MPD)

## add a column indicating significant islands for plotting
sig = (ses.MNNPD.MPD[,"p.value.ranks"] <= 0.05)
neg = (as.numeric(as.character(ses.MNNPD.MPD[,"obs.z"])) <= 0)
ses.MNNPD.MPD <- cbind(ses.MNNPD.MPD, sig)

## code for signficance of observed means, and for direction of pattern (clusering / overdispersed)
ses.MNNPD.MPD[,"Significance"] <- ifelse(ses.MNNPD.MPD[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.MNNPD.MPD[,"obs.z"])) <= 0, 1,
                                        ifelse(ses.MNNPD.MPD[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.MNNPD.MPD[,"obs.z"])) >= 0, 2, 
                                               ifelse(ses.MNNPD.MPD[,"p.value.ranks"] >= 0.05 & as.numeric(as.character(ses.MNNPD.MPD[,"obs.z"])) <= 0, 3, 4)))

pdf("figs/plots/phyloDiv/ses/ses.MNNPD.MPD.IncSize.Final.pdf", width=20, height=10)
ggplot(ses.MNNPD.MPD, aes(x=reorder(factor(island),as.numeric(as.character(Area.m2))), 
                         y=as.numeric(as.character((obs.z))), color=factor(sig), shape=Metric, fill=factor(sig))) +
  geom_point(size=10) +
  coord_cartesian(ylim=c(-5, 5)) + 
  scale_y_continuous("standardized effect size", breaks=seq(-8, 8, 2)) +
  scale_x_discrete("Island (increasing size)") +
  scale_fill_manual(values=alpha(c("FALSE"="white", "TRUE"= "black"), .3), guide="none") +
  scale_color_manual(name =" ",values=c("FALSE"= "grey", "TRUE"="black"), guide="none") +
  scale_shape_manual(name =" ",values=c("MNNPD_inv"= 21, "MPD_inv_nat"= 24), labels=c("MNNPD_inv"="MNNPD i  ", "MPD_inv_nat"="MPD inv_nat")) +
  theme_bw() +
  theme(legend.position="top") +
  geom_abline(intercept = 0, slope = 0, colour = "grey", size = .5) +
  geom_vline(xintercept = c(18.5, 54.5)) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Significane of Phylogenetic Distances for\n invasive species to nearest native (MNNPD),\n and native community (MPDN)\n") +
  theme(plot.title=element_text(size=rel(1.5))) +
  theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
dev.off()

################################## ses PHYLO summary final
sesMNNPD <- ses.MNNPD.MPD[ses.MNNPD.MPD[,"Metric"] =="MNNPD_inv",]
sesMPD <- ses.MNNPD.MPD[ses.MNNPD.MPD[,"Metric"] !="MNNPD_inv",]

sigMNNPDpos = cbind("Island"=rownames(sesMNNPD), "Size.cat"=as.character(sesMNNPD[,"Size.cat"]), 
                    "Metric"=rep(x = "MNNPD inv > 0", times = nrow(sesMNNPD)), 
                    "Significance"=ifelse(sesMNNPD[,"p.value.ranks"] <= 0.05 & 
                                            as.numeric(as.character(sesMNNPD[,"obs.z"])) >= 0, 1,0))
sigMNNPDneg = cbind("Island"=rownames(sesMNNPD), "Size.cat"=as.character(sesMNNPD[,"Size.cat"]), 
                    "Metric"=rep(x = "MNNPD inv < 0", times = nrow(sesMNNPD)), 
                    "Significance"=ifelse(sesMNNPD[,"p.value.ranks"] <= 0.05 & 
                                            as.numeric(as.character(sesMNNPD[,"obs.z"])) <= 0, 1,0))
sigMPDpos = cbind("Island"=rownames(sesMPD), "Size.cat"=as.character(sesMPD[,"Size.cat"]), 
                  "Metric"=rep(x = "MPD inv > 0", times = nrow(sesMPD)), 
                  "Significance"=ifelse(sesMPD[,"p.value.ranks"] <= 0.05 & 
                                          as.numeric(as.character(sesMPD[,"obs.z"])) >= 0, 1,0))
sigMPDneg = cbind("Island"=rownames(sesMPD), "Size.cat"=as.character(sesMPD[,"Size.cat"]), 
                  "Metric"=rep(x = "MPD inv < 0", times = nrow(sesMPD)), 
                  "Significance"=ifelse(sesMPD[,"p.value.ranks"] <= 0.05 & 
                                          as.numeric(as.character(sesMPD[,"obs.z"])) <= 0, 1,0))
sesPhylo <- as.data.frame(rbind(sigMNNPDpos, sigMNNPDneg, sigMPDpos, sigMPDneg))

neworder <- c("MNNPD inv > 0", "MPD inv > 0", "MNNPD inv < 0", "MPD inv < 0")
sesPhylo2 <- arrange(transform(sesPhylo, Metric=factor(Metric,levels=neworder)),Metric)

pdf("figs/plots/phyloDiv/ses/ses.phylo.SummaryBar.pdf")
ggplot(sesPhylo2, aes(x=Significance, fill=factor(Size.cat))) + 
  geom_bar(stat="bin") +
  scale_fill_manual(name="Size Category",
                    breaks=c("sm", "med", "lg"),
                    values=c("sm"="dodgerblue4", "med"="orangered3", "lg"="gold1")) +
  theme_bw() +
  scale_x_discrete(name="", labels=c("0"="NS", "1"="Significant")) +
  facet_wrap(~Metric, ncol = 2) +
  ggtitle("SES Phylogenetic Distances \n") 
dev.off()

#### summarize significance data by size categories

### ses MMNPDi 
MNNPD_inv <- ses.MNNPD.MPD[which(ses.MNNPD.MPD$Metric == "MNNPD_inv"), ] 
nrow(MNNPD_inv[which(as.numeric(as.character(MNNPD_inv$p.value.ranks)) <= 0.05), ]) # 22
nrow(MNNPD_inv[which(as.numeric(as.character(MNNPD_inv$p.value.ranks)) >= 0.05), ]) #50 NS
22/72  # 0.3055556

sigMNNPD_inv <- (MNNPD_inv[which(as.numeric(as.character(MNNPD_inv$p.value.ranks)) <= 0.05), ])
nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigMNNPD_inv, by=0)) # 10/72 small islands sig 0.1388889
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigMNNPD_inv, by=0)) # 6/72 med islands sig  0.08333333
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigMNNPD_inv, by=0)) # 6/72 large isl sig  0.08333333

## ses MPDin 
MPD_inv_nat <- ses.MNNPD.MPD[which(ses.MNNPD.MPD$Metric == "MPD_inv_nat"), ] 
nrow(MPD_inv_nat[which(as.numeric(as.character(MPD_inv_nat$p.value.ranks)) <= 0.05), ]) #27 significant
nrow(MPD_inv_nat[which(as.numeric(as.character(MPD_inv_nat$p.value.ranks)) >= 0.05), ]) #45 NS
27/72  # 0.375

length(which(as.numeric(as.character(MPD_inv_nat$p.value.ranks)) <= 0.05 &
               as.numeric(as.character(MPD_inv_nat[,"obs.z"])) < 0 )) ## 5 MPD clustered, significant
length(which(as.numeric(as.character(MPD_inv_nat$p.value.ranks)) <= 0.05 &
               as.numeric(as.character(MPD_inv_nat[,"obs.z"])) > 0 )) ## 22 MPD overdispersed, significant
5/72 #0.06944444
22/72 #0.3055556

sigMPD_inv <- (MPD_inv_nat[which(as.numeric(as.character(MPD_inv_nat$p.value.ranks)) <= 0.05), ]) #45 NS
nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigMPD_inv, by=0)) # 8/72 small islands sig 0.1111111
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigMPD_inv, by=0)) # 30/72 med islands sig  0.4166667
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigMPD_inv, by=0)) # 12/72 large isl sig  0.1666667
length(which(sigMPD_inv$Size.cat == "sm")) #7/72  0.09722222
length(which(sigMPD_inv$Size.cat == "med")) #16/72   0.2222222
length(which(sigMPD_inv$Size.cat == "lg")) #4/72   0.05555556



###########################################  Significance of Functional Distances  #############################################################
remove.islands.sim <- c("All_SanJuanIslands", "Unnamed_west_of_Castle_Island") 
SJcommNewSim  <- SJcommNew[, -which(names(SJcommNew) %in% remove.islands.sim)]
head(SJcommNewSim)
SJ_islands.sim <- colnames(SJcommNewSim)

####################  Height  #####################

## Null distribution...this will take a while
sim.null.distrib.Height <- lapply(SJ_islands.sim, function(x) sim.meanMNNFD.MFD(phy=SJfinalTree, com=SJcommNewSim, traits=SJtraitLog, island=x, traitname="maxHeight", N = 1000))
names(sim.null.distrib.Height) <- SJ_islands.sim
head(sim.null.distrib.Height)
head(sim.null.distrib.Height["Willow_Island"])
#write.csv(sim.null.distrib.Height, file="sim.null.distrib.Height.null1000.csv")

## Load output of sim.meanMNNFD.MFD()
sdf.height <- read.csv(file="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.Height.null1000.csv", as.is=T, row.names=1)
head(sdf.height[, 1:6])
list.sdf.height <- list()
for (i in 1:ncol(sdf.height)){
  newlist <- (cbind(sdf.height[, c(1:6)]))
  head(newlist)
  colnames(newlist) <- c("n.native.tips", "n.invasive.tips", "meanMNNFDinvasives", "meanMNNFDnatives", "meanMFDinv_nat", "meanMFDnat_nat")
  list.sdf.height[[i]] <- newlist
  sdf.height <- sdf.height[, -c(1:6)]
  
}
names(list.sdf.height) <- SJ_islands.sim
length(list.sdf.height)
names(list.sdf.height[1])
dim(list.sdf.height[[1]])

## Observed Values
obs.MNNFD.maxheight <- lapply(SJ_islands.sim, function(x) functionDistinct(output=phyloObs[[x]], traits=SJtraitLog, traitname="maxHeight")) 
names(obs.MNNFD.maxheight) <- SJ_islands.sim 
head(obs.MNNFD.maxheight)
length(obs.MNNFD.maxheight)
## Summary Observed Values
summ.MNNFD.maxheight  <- lapply(SJ_islands.sim, function(x) functionObsSum(obs.MNNFD.maxheight[[x]])) 
names(summ.MNNFD.maxheight) <- SJ_islands.sim # = phyloObsSum.tmp 
head(summ.MNNFD.maxheight)

metadata
listIslands <- list.sdf.height
list.meta.null.distrib.SJ <- list()
for (i in 1:length(listIslands)){ 
  tmp <- metadata[as.character(names(listIslands[i])), "Area.m2"]
  #tmp <- rep(x=tmp, times=nrow(listIslands[[i]]))
  tmp.MNNFD.i <- rep(summ.MNNFD.maxheight[as.character(names(summ.MNNFD.maxheight[i]))][[1]][[2]], times=nrow(listIslands[[i]]))
  tmp.MFD.in <- rep(summ.MNNFD.maxheight[as.character(names(summ.MNNFD.maxheight[i]))][[1]][[9]], times=nrow(listIslands[[i]]))
  newlist <- mapply(cbind, "islands"=names(listIslands[i]), listIslands[i], "Area.m2"=tmp, "Obs.meanMNNFDinvasives"=tmp.MNNFD.i, "Obs.meanMPFDinv_nat"= tmp.MFD.in, SIMPLIFY=F) 
  list.meta.null.distrib.SJ[i] <- newlist[1]
}
names(list.meta.null.distrib.SJ) <- SJ_islands.sim
summary(list.meta.null.distrib.SJ[1])
head(list.meta.null.distrib.SJ[[71]])

sim.null.distrib.melt <- melt.list(list.meta.null.distrib.SJ, measure.vars="islands")
sim.null.distrib.melt <- sim.null.distrib.melt[which(!sim.null.distrib.melt$value == "NA"),]
head(sim.null.distrib.melt)
pdf("figs/plots/functionDiv/ses/NullInvOcc.MNNFD.Height.pdf", width=20, height=10)
p <- ggplot(sim.null.distrib.melt, aes(x=reorder((L1), Area.m2), y=as.numeric(as.character(meanMNNFDinvasives))), 
            position=position_dodge(width=1)) 
p <- p + geom_boxplot(aes(fill=factor(as.character(variable))), width = 1)
p <- p + scale_x_discrete("Island (increasing size)") #, breaks=seq(0, 80, 10)) 
p <- p + scale_y_discrete("Null distribution mean(MNNFD)")#, breaks=seq(0, 500, 10))  #scale_y_discrete
p <- p + theme_bw() 
#p <- p + theme(axis.ticks = element_blank(), axis.text.y = element_blank())
p <- p + scale_fill_manual(values="grey", labels="")#c("i"= "#F21A00", "n"="#E1AF00") c("i"= "magenta1", "n"="green3")
p <- p + guides(fill=guide_legend(title=" Null distibuiton : random invasive occurrence"))
p <- p + theme(legend.position="top")
#p <- p + ggtitle("Null") + theme(plot.title=element_text(size=rel(1.5)))
p <- p + geom_point(aes(y = Obs.meanMNNFDinvasives), shape=1, color="magenta1") 
p <- p + theme(axis.text.x = element_text(angle = -45, hjust = 0))
p <- p + ggtitle("Null Distribution and observed mean(MNNFD)\n Maximum Height") + theme(plot.title=element_text(size=rel(1.5)))
p <- p + theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
p
dev.off()

pdf("figs/plots/functionDiv/ses/NullInvOcc.MFD.Height.pdf", width=20, height=10)
p <- ggplot(sim.null.distrib.melt, aes(x=reorder((L1), Area.m2), y=as.numeric(as.character(meanMFDinv_nat))), 
            position=position_dodge(width=1)) 
p <- p + geom_boxplot(aes(fill=factor(as.character(variable))), width = 1)
p <- p + scale_x_discrete("Island (increasing size)") #, breaks=seq(0, 80, 10)) 
p <- p + scale_y_log10("Null distribution mean(MFD inv-nat)") # scale_y_discrete, breaks=seq(0, 500, 10)) 
p <- p + theme_bw() 
#p <- p + theme(axis.ticks = element_blank(), axis.text.y = element_blank())
p <- p + scale_fill_manual(values="grey", labels="")#c("i"= "#F21A00", "n"="#E1AF00") c("i"= "magenta1", "n"="green3")
p <- p + guides(fill=guide_legend(title=" Null distibuiton : random invasive occurrence"))
p <- p + theme(legend.position="top")
#p <- p + ggtitle("Null") + theme(plot.title=element_text(size=rel(1.5)))
p <- p + geom_point(aes(y = Obs.meanMPFDinv_nat), shape=1, color="magenta1") 
p <- p + theme(axis.text.x = element_text(angle = -45, hjust = 0))
p <- p + ggtitle("Null Distribution and observed mean(MFD inv-nat)\n  Maximum Height") + theme(plot.title=element_text(size=rel(1.5)))
p <- p + theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
p
dev.off()

#### Summarize simualted means, standardized effect size 
ses.SanJuan.MNNFD.MFD.height <- lapply(SJ_islands.sim, function(x) ses.FunctionDist(phy=SJfinalTree, com=SJcommNewSim, island=x, 
                                                                                    simOneIslandOneTrait=list.sdf.height, outputMNNPD=phyloObs[[x]], traits=SJtraitLog, traitname="maxHeight", N=1000))
names(ses.SanJuan.MNNFD.MFD.height) <- SJ_islands.sim

metadata
listIslands <- ses.SanJuan.MNNFD.MFD.height
ses.SJ.MNNFD.height <- data.frame()
for (i in 1:length(listIslands)){ 
  if (length(listIslands[[i]]) != 2){
    newlist <- NULL
  } else {
    tmp1 <- metadata[as.character(names(listIslands[i])), "Area.m2"]
    tmp2 <- metadata[as.character(names(listIslands[i])), "Size.cat"]
    newlist <-cbind(t(listIslands[[i]]), "Area.m2"=tmp1, "Size.cat"=tmp2) 
  }
  ses.SJ.MNNFD.height <- rbind(ses.SJ.MNNFD.height, newlist)
  
}
head(ses.SJ.MNNFD.height)
dim(ses.SJ.MNNFD.height) #144
ses.SJ.MNNFD.height <- na.omit(ses.SJ.MNNFD.height)
ses.SJ.MNNFD.height$p.value.ranks <- as.numeric(as.character(ses.SJ.MNNFD.height$p.value.ranks)) 
ses.SJ.MNNFD.height.sig <- subset(ses.SJ.MNNFD.height, p.value.ranks <= 0.05)
dim(ses.SJ.MNNFD.height.sig)
sigH = (ses.SJ.MNNFD.height[,"p.value.ranks"] <= 0.05)
ses.SJ.MNNFD.height <- cbind(ses.SJ.MNNFD.height, sigH)

ses.SJ.MNNFD.height[,"Significance"] <- ifelse(ses.SJ.MNNFD.height[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.SJ.MNNFD.height[,"obs.z"])) <= 0, 1,
                                               ifelse(ses.SJ.MNNFD.height[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.SJ.MNNFD.height[,"obs.z"])) >= 0, 2, 
                                                      ifelse(ses.SJ.MNNFD.height[,"p.value.ranks"] >= 0.05 & as.numeric(as.character(ses.SJ.MNNFD.height[,"obs.z"])) <= 0, 3, 4)))

pdf("figs/plots/functionDiv/ses/ses.SJ.MNNFD.height.FINAL.pdf", width=20, height=10)
ggplot(ses.SJ.MNNFD.height, aes(x=reorder(factor(island),as.numeric(as.character(Area.m2))), 
                                y=as.numeric(as.character((obs.z))), color=factor(sigH), shape=Metric, fill=factor(sigH))) +
  geom_point(size=10) + #, aes(shape = as.numeric(as.character(p.value.ranks >= 0.05)))) +  
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
  ggtitle("Significane of Maximum Height difference for\n invasive species to nearest native (NNFDi),\n and native community (MFDNi)\n") +
  theme(plot.title=element_text(size=rel(1.5))) +
  theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
dev.off()

#### summarize significance data by size categories
ses.SJ.MNNFD.height$Size.cat <- as.character(ses.SJ.MNNFD.height$Size.cat)
ses.SJ.MNNFD.height$Size.cat <- factor(ses.SJ.MNNFD.height$Size.cat, levels=c("sm", "med", "lg"))

##### ses MMNPDi 
MNNFD_inv <- na.omit(ses.SJ.MNNFD.height[which(ses.SJ.MNNFD.height$Metric == "MNNFD_inv"), ] )
nrow(MNNFD_inv[which(as.numeric(as.character(MNNFD_inv$p.value.ranks)) <= 0.05), ]) #7
nrow(MNNFD_inv[which(as.numeric(as.character(MNNFD_inv$p.value.ranks)) >= 0.05), ]) #64 NS
7/71 #0.09859155
#### ses MPDin 
MFD_inv_nat <- ses.SJ.MNNFD.height[which(ses.SJ.MNNFD.height$Metric == "MFD_inv_nat"), ] 
nrow(MFD_inv_nat[which(as.numeric(as.character(MFD_inv_nat$p.value.ranks)) <= 0.05), ]) #13
nrow(MFD_inv_nat[which(as.numeric(as.character(MFD_inv_nat$p.value.ranks)) >= 0.05), ]) #58 NS
13/71 #0.1830986

sigheightMNNFD_inv <- (MNNFD_inv[which(as.numeric(as.character(MNNFD_inv$p.value.ranks)) <= 0.05), ]) #45 NS
sigheightMPD_inv <- (MFD_inv_nat[which(as.numeric(as.character(MFD_inv_nat$p.value.ranks)) <= 0.05), ]) #45 NS

length(which(sigheightMNNFD_inv$Size.cat == "sm")) #3/71  0.04225352
length(which(sigheightMNNFD_inv$Size.cat == "med")) #3/71   0.04225352
length(which(sigheightMNNFD_inv$Size.cat == "lg")) #1/71   0.01408451

length(which(sigheightMPD_inv$Size.cat == "sm")) #2/71  0.02816901
length(which(sigheightMPD_inv$Size.cat == "med")) #7/71   0.09859155
length(which(sigheightMPD_inv$Size.cat == "lg")) #4/71   0.05633803



####### SLA
## Null distribution 
sim.null.distrib.SLA <- lapply(SJ_islands.sim, function(x) sim.meanMNNFD.MFD(phy=SJfinalTree, com=SJcommNewSim, traits=SJtraitLog, island=x, traitname="sla", N = 5))
names(sim.null.distrib.SLA) <- SJ_islands.sim
head(sim.null.distrib.SLA)
head(sim.null.distrib.SLA["Willow_Island"])
#write.csv(sim.null.distrib.SLA, file="sim.null.distrib.SLA.1000.csv")
sdf.SLA <- read.csv(file="~/Dropbox/Work/TankLab/Projects/SanJuans/Manuscript/Drafts/Figs.v4/sim.null.distrib.SLA.1000.csv", as.is=T, row.names=1)
head(sdf.SLA[, 1:6])
list.sdf.SLA <- list()
for (i in 1:ncol(sdf.SLA)){
  newlist <- (cbind(sdf.SLA[, c(1:6)]))
  head(newlist)
  colnames(newlist) <- c("n.native.tips", "n.invasive.tips", "meanMNNFDinvasives", "meanMNNFDnatives", "meanMFDinv_nat", "meanMFDnat_nat")
  list.sdf.SLA[[i]] <- newlist
  sdf.SLA <- sdf.SLA[, -c(1:6)]
  
}
names(list.sdf.SLA) <- SJ_islands.sim
names(list.sdf.SLA[1])
dim(list.sdf.SLA[[1]])
## SES
ses.SanJuan.MNNFD.MFD.SLA <- lapply(SJ_islands.sim, function(x) ses.FunctionDist(phy=SJfinalTree, com=SJcommNewSim, island=x, 
                                                                                 simOneIslandOneTrait=list.sdf.SLA, outputMNNPD=phyloObs[[x]], traits=SJtraitLog, traitname="sla", N=1000))
names(ses.SanJuan.MNNFD.MFD.SLA) <- SJ_islands.sim
metadata
listIslands <- ses.SanJuan.MNNFD.MFD.SLA
ses.SJ.MNNFD.SLA <- data.frame()
for (i in 1:length(listIslands)){ 
  if (length(listIslands[[i]]) != 2){
    newlist <- NULL
  } else {
    tmp1 <- metadata[as.character(names(listIslands[i])), "Area.m2"]
    tmp2 <- metadata[as.character(names(listIslands[i])), "Size.cat"]
    newlist <-cbind(t(listIslands[[i]]), "Area.m2"=tmp1, "Size.cat"=tmp2) 
  }
  ses.SJ.MNNFD.SLA <- rbind(ses.SJ.MNNFD.SLA, newlist)
}
head(ses.SJ.MNNFD.SLA)
dim(ses.SJ.MNNFD.SLA)
ses.SJ.MNNFD.SLA <- na.omit(ses.SJ.MNNFD.SLA)
ses.SJ.MNNFD.SLA$p.value.ranks <- as.numeric(as.character(ses.SJ.MNNFD.SLA$p.value.ranks)) 
ses.SJ.MNNFD.SLA.sig <- subset(ses.SJ.MNNFD.SLA, p.value.ranks <= 0.05)

sigSLA = (ses.SJ.MNNFD.SLA[,"p.value.ranks"] <= 0.05)
ses.SJ.MNNFD.SLA <- cbind(ses.SJ.MNNFD.SLA, sigSLA)
ses.SJ.MNNFD.SLA[,"Significance"] <- ifelse(ses.SJ.MNNFD.SLA[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.SJ.MNNFD.SLA[,"obs.z"])) <= 0, 1,
                                            ifelse(ses.SJ.MNNFD.SLA[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.SJ.MNNFD.SLA[,"obs.z"])) >= 0, 2, 
                                                   ifelse(ses.SJ.MNNFD.SLA[,"p.value.ranks"] >= 0.05 & as.numeric(as.character(ses.SJ.MNNFD.SLA[,"obs.z"])) <= 0, 3, 4)))



pdf("ses.SJ.MNNFD.SLA.pdf", width=20, height=10)
ggplot(ses.SJ.MNNFD.SLA, aes(x=reorder(factor(island),as.numeric(as.character(Area.m2))), 
                             y=as.numeric(as.character((obs.z))), color=Metric, shape=Metric, size=Size.cat)) +
  geom_point() + #, aes(shape = as.numeric(as.character(p.value.ranks >= 0.05)))) +  
  geom_point(colour=alpha("red", 0.5), shape="*", data = ses.SJ.MNNFD.SLA.sig, size = 5) +
  #geom_point(colour=alpha("blue", 0.5), shape="*", data = ses.SJ.MNNFD.SLA.sig.oneminus, size = 10) +
  coord_cartesian(ylim=c(-5, 5)) + scale_y_continuous("standardized effect size", breaks=seq(-8, 8, 2)) +
  scale_x_discrete("Island (increasing size)") +
  scale_color_manual(name =" ",values=c("MNNFD_inv"= "grey", "MFD_inv_nat"="black"), labels=c("MNNFD_inv"= "NNFDi  ", "MFD_inv_nat"="MFDNi")) +
  scale_shape_manual(name =" ",values=c("MNNFD_inv"= 1, "MFD_inv_nat"= 2), labels=c("MNNFD_inv"= "NNFDi  ", "MFD_inv_nat"="MFDNi")) +
  scale_size_manual(values=c("sm"= 2, "med"= 4, "lg"=10), labels=c("sm"= 'Small', "med"= 'Medium', "lg"='Large'), breaks=c("sm", "med", "lg")) +
  theme_bw() +
  theme(legend.position="top") +
  geom_abline(intercept = 0, slope = 0, colour = "grey", size = .5) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Significane of SLA differnence for\n invasive species to nearest native (NNFDi),\n and native community (MFDNi)") +
  theme(plot.title=element_text(size=rel(1.5))) +
  theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
dev.off()

ses.SJ.MNNFD.SLA.four <- ses.SJ.MNNFD.SLA[ses.SJ.MNNFD.SLA$island %in% four.islands, ]
ses.SJ.MNNFD.SLA.sig.four <- ses.SJ.MNNFD.SLA.sig[ses.SJ.MNNFD.SLA.sig$island %in% four.islands, ]
pdf("four.islands.ses.SJ.MNNFD.SLA.pdf")
ggplot(ses.SJ.MNNFD.SLA.four, aes(x=reorder(factor(island),as.numeric(as.character(Area.m2))), 
                                  y=as.numeric(as.character((obs.z))), color=Metric, shape=Metric)) +
  geom_point(size=5) + #, aes(shape = as.numeric(as.character(p.value.ranks >= 0.05)))) +  
  #geom_point(colour=alpha("red", 0.5), shape="*", data = ses.SJ.MNNFD.SLA.sig.four, size = 10) +
  #geom_point(colour=alpha("blue", 0.5), shape="*", data = ses.SJ.MNNFD.height.sig.oneminus, size = 10) +
  coord_cartesian(ylim=c(-5, 5)) + scale_y_continuous("standardized effect size", breaks=seq(-8, 8, 2)) +
  scale_x_discrete("Island (increasing size)") +
  scale_color_manual(name =" ",values=c("MNNFD_inv"= "grey", "MFD_inv_nat"="black"), labels=c("MNNFD_inv"="MNNFD i  ", "MFD_inv_nat"="MFD inv_nat")) +
  scale_shape_manual(name =" ",values=c("MNNFD_inv"= 1, "MFD_inv_nat"= 2), labels=c("MNNFD_inv"="MNNFD i  ", "MFD_inv_nat"="MFD inv_nat")) +
  theme_bw() +
  theme(legend.position="top") +
  geom_abline(intercept = 0, slope = 0, colour = "grey", size = .5) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Significane of SLA Distances for\n invasive species to nearest native(MNNFD),\n and native community (MFD)") +
  theme(plot.title=element_text(size=rel(1.5))) +
  theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
dev.off()

pdf("ses.SJ.MNNFD.SLA.FINAL.pdf", width=20, height=10)
ggplot(ses.SJ.MNNFD.SLA, aes(x=reorder(factor(island),as.numeric(as.character(Area.m2))), 
                             y=as.numeric(as.character((obs.z))), color=factor(sigSLA), shape=Metric, fill=factor(sigSLA))) +
  geom_point(size=10) + #, aes(shape = as.numeric(as.character(p.value.ranks >= 0.05)))) +  
  coord_cartesian(ylim=c(-5, 5)) + 
  scale_y_continuous("standardized effect size", breaks=seq(-8, 8, 2)) +
  scale_x_discrete("Island (increasing size)") +
  scale_fill_manual(values=alpha(c("FALSE"="white", "TRUE"= "black"), .3), legend=F) +
  scale_color_manual(name =" ",values=c("FALSE"= "grey", "TRUE"="black"), legend=F) +
  scale_shape_manual(name =" ",values=c("MNNFD_inv"= 21, "MFD_inv_nat"= 24), labels=c("MNNFD_inv"="MNNFD i  ", "MFD_inv_nat"="MFD_inv_nat")) +
  theme_bw() +
  theme(legend.position="top") +
  geom_abline(intercept = 0, slope = 0, colour = "grey", size = .5) +
  geom_vline(xintercept = c(9.5, 42.5)) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Significane of SLA difference for\n invasive species to nearest native (NNFDi),\n and native community (MFDNi)\n") +
  theme(plot.title=element_text(size=rel(1.5))) +
  theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
dev.off()

#### summarize significance data by size categories
ses.SJ.MNNFD.SLA$Size.cat <- as.character(ses.SJ.MNNFD.SLA$Size.cat)
ses.SJ.MNNFD.SLA$Size.cat <- factor(ses.SJ.MNNFD.SLA$Size.cat, levels=c("sm", "med", "lg"))

pdf("ses.MNNFD.SLA.SummaryBar.pdf")
ggplot(ses.SJ.MNNFD.SLA, aes(x=Size.cat, fill=factor(Significance))) + 
  geom_bar(stat="bin") +
  facet_wrap(~ Metric) +
  scale_fill_manual(name="Significance",
                    breaks=c("1", "2", "3", "4"),
                    labels=c("Significantly Clustered", "Significantly Overdispersed", 
                             "NS Clustered", "NS Overdispersed"),
                    values=c("1"="black", "2"="orange", "3"="skyblue", "4"="springgreen4")) +
  theme_bw() +
  scale_x_discrete(name="") +
  #theme(legend.position="top") +
  ggtitle("Significance of SLA Difference for\n invasive species to nearest native (NNFD_inv),\n and native community (MFDN_inv)\n") 
dev.off()



##### ses MMNPDi 
MNNFD_inv <- na.omit(ses.SJ.MNNFD.SLA[which(ses.SJ.MNNFD.SLA$Metric == "MNNFD_inv"), ] )
nrow(MNNFD_inv[which(as.numeric(as.character(MNNFD_inv$p.value.ranks)) <= 0.05), ]) #0
nrow(MNNFD_inv[which(as.numeric(as.character(MNNFD_inv$p.value.ranks)) >= 0.05), ]) #60 NS
#### ses MPDin 
MFD_inv_nat <- ses.SJ.MNNFD.SLA[which(ses.SJ.MNNFD.SLA$Metric == "MFD_inv_nat"), ] 
nrow(MFD_inv_nat[which(as.numeric(as.character(MFD_inv_nat$p.value.ranks)) <= 0.05), ]) #1
nrow(MFD_inv_nat[which(as.numeric(as.character(MFD_inv_nat$p.value.ranks)) >= 0.05), ]) #59 NS
1/59

sigSLAMNNFD_inv <- (MNNFD_inv[which(as.numeric(as.character(MNNFD_inv$p.value.ranks)) <= 0.05), ]) #45 NS
sigSLAMPD_inv <- (MFD_inv_nat[which(as.numeric(as.character(MFD_inv_nat$p.value.ranks)) <= 0.05), ]) #45 NS

length(which(sigSLAMNNFD_inv$Size.cat == "sm")) #0
length(which(sigSLAMNNFD_inv$Size.cat == "med")) #0
length(which(sigSLAMNNFD_inv$Size.cat == "lg")) #0

length(which(sigSLAMPD_inv$Size.cat == "sm")) #0
length(which(sigSLAMPD_inv$Size.cat == "med")) #1/60 0.01666667
length(which(sigSLAMPD_inv$Size.cat == "lg")) #0


####### leaflet
## Null distribution 
sim.null.distrib.leafletSize <- lapply(SJ_islands.sim, function(x) sim.meanMNNFD.MFD(phy=SJfinalTree, com=SJcommNewSim, traits=SJtraitLog, island=x, traitname="leafletSize", N = 5))
names(sim.null.distrib.leafletSize) <- SJ_islands.sim
head(sim.null.distrib.leafletSize)
head(sim.null.distrib.leafletSize["Willow_Island"])
#write.csv(sim.null.distrib.leafletSize, file="sim.null.distrib.leafletSize.1000.csv")
sdf.leaflet <- read.csv(file="~/Dropbox/Work/TankLab/Projects/SanJuans/Manuscript/Drafts/Figs.v4/sim.null.distrib.leafletSize.1000.csv", as.is=T, row.names=1)
head(sdf.leaflet[, 1:6])
list.sdf.leaflet <- list()
for (i in 1:ncol(sdf.leaflet)){
  newlist <- (cbind(sdf.leaflet[, c(1:6)]))
  head(newlist)
  colnames(newlist) <- c("n.native.tips", "n.invasive.tips", "meanMNNFDinvasives", "meanMNNFDnatives", "meanMFDinv_nat", "meanMFDnat_nat")
  list.sdf.leaflet[[i]] <- newlist
  sdf.leaflet <- sdf.leaflet[, -c(1:6)]
  
}
names(list.sdf.leaflet) <- SJ_islands.sim
names(list.sdf.leaflet[1])
dim(list.sdf.leaflet[[1]])
## SES
ses.SanJuan.MNNFD.MFD.leaflet <- lapply(SJ_islands.sim, function(x) ses.FunctionDist(phy=SJfinalTree, com=SJcommNewSim, island=x, 
                                                                                     simOneIslandOneTrait=list.sdf.leaflet, outputMNNPD=phyloObs[[x]], traits=SJtraitLog, traitname="leafletSize", N=1000))
names(ses.SanJuan.MNNFD.MFD.leaflet) <- SJ_islands.sim
metadata
listIslands <- ses.SanJuan.MNNFD.MFD.leaflet
ses.SJ.MNNFD.leaflet <- data.frame()
for (i in 1:length(listIslands)){ 
  if (length(listIslands[[i]]) != 2){
    newlist <- NULL
  } else {
    tmp1 <- metadata[as.character(names(listIslands[i])), "Area.m2"]
    tmp2 <- metadata[as.character(names(listIslands[i])), "Size.cat"]
    newlist <-cbind(t(listIslands[[i]]), "Area.m2"=tmp1, "Size.cat"=tmp2) 
  }
  ses.SJ.MNNFD.leaflet <- rbind(ses.SJ.MNNFD.leaflet, newlist)
}
head(ses.SJ.MNNFD.leaflet)
dim(ses.SJ.MNNFD.leaflet)
ses.SJ.MNNFD.leaflet <- na.omit(ses.SJ.MNNFD.leaflet)
ses.SJ.MNNFD.leaflet$p.value.ranks <- as.numeric(as.character(ses.SJ.MNNFD.leaflet$p.value.ranks)) 
ses.SJ.MNNFD.leaflet.sig <- subset(ses.SJ.MNNFD.leaflet, p.value.ranks <= 0.05)

sigLeaflet = (ses.SJ.MNNFD.leaflet[,"p.value.ranks"] <= 0.05)
ses.SJ.MNNFD.leaflet <- cbind(ses.SJ.MNNFD.leaflet, sigLeaflet)

ses.SJ.MNNFD.leaflet[,"Significance"] <- ifelse(ses.SJ.MNNFD.leaflet[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.SJ.MNNFD.leaflet[,"obs.z"])) <= 0, 1,
                                                ifelse(ses.SJ.MNNFD.leaflet[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.SJ.MNNFD.leaflet[,"obs.z"])) >= 0, 2, 
                                                       ifelse(ses.SJ.MNNFD.leaflet[,"p.value.ranks"] >= 0.05 & as.numeric(as.character(ses.SJ.MNNFD.leaflet[,"obs.z"])) <= 0, 3, 4)))

pdf("ses.SJ.MNNFD.leaflet.pdf", width=20, height=10)
ggplot(ses.SJ.MNNFD.leaflet, aes(x=reorder(factor(island),as.numeric(as.character(Area.m2))), 
                                 y=as.numeric(as.character((obs.z))), color=Metric, shape=Metric, size=Size.cat)) +
  geom_point() + #, aes(shape = as.numeric(as.character(p.value.ranks >= 0.05)))) +  
  geom_point(colour=alpha("red", 0.5), shape="*", data = ses.SJ.MNNFD.leaflet.sig, size = 5) +
  #geom_point(colour=alpha("blue", 0.5), shape="*", data = ses.SJ.MNNFD.leaflet.sig.oneminus, size = 10) +
  coord_cartesian(ylim=c(-5, 5)) + scale_y_continuous("standardized effect size", breaks=seq(-8, 8, 2)) +
  scale_x_discrete("Island (increasing size)") +
  scale_color_manual(name =" ",values=c("MNNFD_inv"= "grey", "MFD_inv_nat"="black"), labels=c("MNNFD_inv"="NNFDi  ", "MFD_inv_nat"="MFDNi")) +
  scale_shape_manual(name =" ",values=c("MNNFD_inv"= 1, "MFD_inv_nat"= 2), labels=c("MNNFD_inv"="NNFDi  ", "MFD_inv_nat"="MFDNi")) +
  scale_size_manual(values=c("sm"= 2, "med"= 4, "lg"=10), labels=c("sm"= 'Small', "med"= 'Medium', "lg"='Large'), breaks=c("sm", "med", "lg")) +
  theme_bw() +
  theme(legend.position="top") +
  geom_abline(intercept = 0, slope = 0, colour = "grey", size = .5) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Significane of Leaf Size difference for\n invasive species to nearest native (NNFDi),\n and native community (MFDNi)") +
  theme(plot.title=element_text(size=rel(1.5))) +
  theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
dev.off()

ses.SJ.MNNFD.leaflet.four <- ses.SJ.MNNFD.leaflet[ses.SJ.MNNFD.leaflet$island %in% four.islands, ]
ses.SJ.MNNFD.leaflet.sig.four <- ses.SJ.MNNFD.leaflet.sig[ses.SJ.MNNFD.leaflet.sig$island %in% four.islands, ]
pdf("four.islands.ses.SJ.MNNFD.leaflet.four.pdf")
ggplot(ses.SJ.MNNFD.leaflet.four, aes(x=reorder(factor(island),as.numeric(as.character(Area.m2))), 
                                      y=as.numeric(as.character((obs.z))), color=Metric, shape=Metric)) +
  geom_point(size=5) + #, aes(shape = as.numeric(as.character(p.value.ranks >= 0.05)))) +  
  #geom_point(colour=alpha("red", 0.5), shape="*", data = ses.SJ.MNNFD.SLA.sig.four, size = 10) +
  #geom_point(colour=alpha("blue", 0.5), shape="*", data = ses.SJ.MNNFD.height.sig.oneminus, size = 10) +
  coord_cartesian(ylim=c(-5, 5)) + scale_y_continuous("standardized effect size", breaks=seq(-8, 8, 2)) +
  scale_x_discrete("Island (increasing size)") +
  scale_color_manual(name =" ",values=c("MNNFD_inv"= "grey", "MFD_inv_nat"="black"), labels=c("MNNFD_inv"="MNNFD i  ", "MFD_inv_nat"="MFD inv_nat")) +
  scale_shape_manual(name =" ",values=c("MNNFD_inv"= 1, "MFD_inv_nat"= 2), labels=c("MNNFD_inv"="MNNFD i  ", "MFD_inv_nat"="MFD inv_nat")) +
  theme_bw() +
  theme(legend.position="top") +
  geom_abline(intercept = 0, slope = 0, colour = "grey", size = .5) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Significane of leaf size Distances for\n invasive species to nearest native(MNNFD),\n and native community (MFD)") +
  theme(plot.title=element_text(size=rel(1.5))) +
  theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
dev.off()

pdf("ses.SJ.MNNFD.leaflet.FINAL.pdf", width=20, height=10)
ggplot(ses.SJ.MNNFD.leaflet, aes(x=reorder(factor(island),as.numeric(as.character(Area.m2))), 
                                 y=as.numeric(as.character((obs.z))), color=factor(sigLeaflet), shape=Metric, fill=factor(sigLeaflet))) +
  geom_point(size=10) + #, aes(shape = as.numeric(as.character(p.value.ranks >= 0.05)))) +  
  coord_cartesian(ylim=c(-5, 5)) + 
  scale_y_continuous("standardized effect size", breaks=seq(-8, 8, 2)) +
  scale_x_discrete("Island (increasing size)") +
  scale_fill_manual(values=alpha(c("FALSE"="white", "TRUE"= "black"), .3), legend=F) +
  scale_color_manual(name =" ",values=c("FALSE"= "grey", "TRUE"="black"), legend=F) +
  scale_shape_manual(name =" ",values=c("MNNFD_inv"= 21, "MFD_inv_nat"= 24), labels=c("MNNFD_inv"="MNNFD i  ", "MFD_inv_nat"="MFD_inv_nat")) +
  theme_bw() +
  theme(legend.position="top") +
  geom_abline(intercept = 0, slope = 0, colour = "grey", size = .5) +
  geom_vline(xintercept = c(8.5, 38.5)) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Significane of Leaflet Size difference for\n invasive species to nearest native (NNFDi),\n and native community (MFDNi)\n") +
  theme(plot.title=element_text(size=rel(1.5))) +
  theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
dev.off()

#### summarize significance data by size categories
ses.SJ.MNNFD.leaflet$Size.cat <- as.character(ses.SJ.MNNFD.leaflet$Size.cat)
ses.SJ.MNNFD.leaflet$Size.cat <- factor(ses.SJ.MNNFD.leaflet$Size.cat, levels=c("sm", "med", "lg"))

pdf("ses.MNNFD.leaflet.SummaryBar.pdf")
ggplot(ses.SJ.MNNFD.leaflet, aes(x=Size.cat, fill=factor(Significance))) + 
  geom_bar(stat="bin") +
  facet_wrap(~ Metric) +
  scale_fill_manual(name="Significance",
                    breaks=c("1", "2", "3", "4"),
                    labels=c("Significantly Clustered", "Significantly Overdispersed", 
                             "NS Clustered", "NS Overdispersed"),
                    values=c("1"="black", "2"="orange", "3"="skyblue", "4"="springgreen4")) +
  theme_bw() +
  scale_x_discrete(name="") +
  #theme(legend.position="top") +
  ggtitle("Significance of Leaf Difference for\n invasive species to nearest native (NNFD_inv),\n and native community (MFDN_inv)\n") 
dev.off()

##### ses MMNPDi 
MNNFD_inv <- na.omit(ses.SJ.MNNFD.leaflet[which(ses.SJ.MNNFD.leaflet$Metric == "MNNFD_inv"), ] )
nrow(MNNFD_inv[which(as.numeric(as.character(MNNFD_inv$p.value.ranks)) <= 0.05), ]) #6
nrow(MNNFD_inv[which(as.numeric(as.character(MNNFD_inv$p.value.ranks)) >= 0.05), ]) #50 NS
6/56
#### ses MPDin 
MFD_inv_nat <- ses.SJ.MNNFD.leaflet[which(ses.SJ.MNNFD.leaflet$Metric == "MFD_inv_nat"), ] 
nrow(MFD_inv_nat[which(as.numeric(as.character(MFD_inv_nat$p.value.ranks)) <= 0.05), ]) #5
nrow(MFD_inv_nat[which(as.numeric(as.character(MFD_inv_nat$p.value.ranks)) >= 0.05), ]) #51 NS
5/56

sigLeafletMNNFD_inv <- (MNNFD_inv[which(as.numeric(as.character(MNNFD_inv$p.value.ranks)) <= 0.05), ]) #45 NS
sigLeafletMPD_inv <- (MFD_inv_nat[which(as.numeric(as.character(MFD_inv_nat$p.value.ranks)) <= 0.05), ]) #45 NS

length(which(sigLeafletMNNFD_inv$Size.cat == "sm")) #2/56  0.03571429
length(which(sigLeafletMNNFD_inv$Size.cat == "med")) #4/56   0.07142857
length(which(sigLeafletMNNFD_inv$Size.cat == "lg")) #0   

length(which(sigLeafletMPD_inv$Size.cat == "sm")) #1/56  0.01785714
length(which(sigLeafletMPD_inv$Size.cat == "med")) #1/56   0.01785714
length(which(sigLeafletMPD_inv$Size.cat == "lg")) #3/56   0.05357143


####### leaf N
## Null distribution 
sim.null.distrib.leafN <- lapply(SJ_islands.sim, function(x) sim.meanMNNFD.MFD(phy=SJfinalTree, com=SJcommNewSim, traits=SJtraitLog, island=x, traitname="leafN", N = 5))
names(sim.null.distrib.leafN) <- SJ_islands.sim
head(sim.null.distrib.leafN)
head(sim.null.distrib.leafN["Willow_Island"])
#write.csv(sim.null.distrib.Height, file="sim.null.distrib.leafN.1000.csv")
sdf.leafN <- read.csv(file="~/Dropbox/Work/TankLab/Projects/SanJuans/Manuscript/Drafts/Figs.v4/sim.null.distrib.leafN.1000.csv", as.is=T, row.names=1)
head(sdf.leafN[, 1:6])
list.sdf.leafN <- list()
for (i in 1:ncol(sdf.leafN)){
  newlist <- (cbind(sdf.leafN[, c(1:6)]))
  head(newlist)
  colnames(newlist) <- c("n.native.tips", "n.invasive.tips", "meanMNNFDinvasives", "meanMNNFDnatives", "meanMFDinv_nat", "meanMFDnat_nat")
  list.sdf.leafN[[i]] <- newlist
  sdf.leafN <- sdf.leafN[, -c(1:6)]
  
}
names(list.sdf.leafN) <- SJ_islands.sim
names(list.sdf.leafN[1])
head(list.sdf.leafN[[1]])
## SES
ses.SanJuan.MNNFD.MFD.leafN <- lapply(SJ_islands.sim, function(x) ses.FunctionDist(phy=SJfinalTree, com=SJcommNewSim, island=x,
                                                                                   simOneIslandOneTrait=list.sdf.leafN, outputMNNPD=phyloObs[[x]], traits=SJtraitLog, traitname="leafN", N=1000))
names(ses.SanJuan.MNNFD.MFD.leafN) <- SJ_islands.sim
metadata
listIslands <- ses.SanJuan.MNNFD.MFD.leafN
ses.SJ.MNNFD.leafN <- data.frame()
for (i in 1:length(listIslands)){ 
  if (length(listIslands[[i]]) != 2){
    newlist <- NULL
  } else {
    tmp1 <- metadata[as.character(names(listIslands[i])), "Area.m2"]
    tmp2 <- metadata[as.character(names(listIslands[i])), "Size.cat"]
    newlist <-cbind(t(listIslands[[i]]), "Area.m2"=tmp1, "Size.cat"=tmp2) 
  }
  ses.SJ.MNNFD.leafN <- rbind(ses.SJ.MNNFD.leafN, newlist)
}
head(ses.SJ.MNNFD.leafN)
dim(ses.SJ.MNNFD.leafN)
ses.SJ.MNNFD.leafN <- na.omit(ses.SJ.MNNFD.leafN)
ses.SJ.MNNFD.leafN$p.value.ranks <- as.numeric(as.character(ses.SJ.MNNFD.leafN$p.value.ranks)) 
ses.SJ.MNNFD.leafN.sig <- subset(ses.SJ.MNNFD.leafN, p.value.ranks <= 0.05)

sigLN = (ses.SJ.MNNFD.leafN[,"p.value.ranks"] <= 0.05)
ses.SJ.MNNFD.leafN <- cbind(ses.SJ.MNNFD.leafN, sigLN)

ses.SJ.MNNFD.leafN[,"Significance"] <- ifelse(ses.SJ.MNNFD.leafN[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.SJ.MNNFD.leafN[,"obs.z"])) <= 0, 1,
                                              ifelse(ses.SJ.MNNFD.leafN[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.SJ.MNNFD.leafN[,"obs.z"])) >= 0, 2, 
                                                     ifelse(ses.SJ.MNNFD.leafN[,"p.value.ranks"] >= 0.05 & as.numeric(as.character(ses.SJ.MNNFD.leafN[,"obs.z"])) <= 0, 3, 4)))

pdf("ses.SJ.MNNFD.leafN.pdf", width=20, height=10)
ggplot(ses.SJ.MNNFD.leafN, aes(x=reorder(factor(island),as.numeric(as.character(Area.m2))), 
                               y=as.numeric(as.character((obs.z))), color=Metric, shape=Metric, size=Size.cat)) +
  geom_point() + #, aes(shape = as.numeric(as.character(p.value.ranks >= 0.05)))) +  
  geom_point(colour=alpha("red", 0.5), shape="*", data = ses.SJ.MNNFD.leafN.sig, size = 5) +
  coord_cartesian(ylim=c(-5, 5)) + scale_y_continuous("standardized effect size", breaks=seq(-8, 8, 2)) +
  scale_x_discrete("Island (increasing size)") +
  scale_color_manual(name =" ",values=c("MNNFD_inv"= "grey", "MFD_inv_nat"="black"), labels=c("MNNFD_inv"="NNFDi  ", "MFD_inv_nat"="MFDNi")) +
  scale_shape_manual(name =" ",values=c("MNNFD_inv"= 1, "MFD_inv_nat"= 2), labels=c("MNNFD_inv"="NNFDi  ", "MFD_inv_nat"="MFDNi")) +
  scale_size_manual(values=c("sm"= 2, "med"= 4, "lg"=10), labels=c("sm"= 'Small', "med"= 'Medium', "lg"='Large'), breaks=c("sm", "med", "lg")) +
  theme_bw() +
  theme(legend.position="top") +
  geom_abline(intercept = 0, slope = 0, colour = "grey", size = .5) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Significane of Leaf Nitrogen differnece for\n invasive species to nearest native (NNFDi),\n and native community (MFDNi)") +
  theme(plot.title=element_text(size=rel(1.5))) +
  theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
dev.off()

ses.SJ.MNNFD.leafN.four <- ses.SJ.MNNFD.leafN[ses.SJ.MNNFD.leafN$island %in% four.islands, ]
ses.SJ.MNNFD.leafN.sig.four <- ses.SJ.MNNFD.leafN.sig[ses.SJ.MNNFD.leafN.sig$island %in% four.islands, ]
pdf("four.islands.ses.SJ.MNNFD.leafN.four.pdf")
ggplot(ses.SJ.MNNFD.leafN.four, aes(x=reorder(factor(island),as.numeric(as.character(Area.m2))), 
                                    y=as.numeric(as.character((obs.z))), color=Metric, shape=Metric)) +
  geom_point(size=5) + #, aes(shape = as.numeric(as.character(p.value.ranks >= 0.05)))) +  
  #geom_point(colour=alpha("red", 0.5), shape="*", data = ses.SJ.MNNFD.SLA.sig.four, size = 10) +
  #geom_point(colour=alpha("blue", 0.5), shape="*", data = ses.SJ.MNNFD.height.sig.oneminus, size = 10) +
  coord_cartesian(ylim=c(-5, 5)) + scale_y_continuous("standardized effect size", breaks=seq(-8, 8, 2)) +
  scale_x_discrete("Island (increasing size)") +
  scale_color_manual(name =" ",values=c("MNNFD_inv"= "grey", "MFD_inv_nat"="black"), labels=c("MNNFD_inv"="MNNFD i  ", "MFD_inv_nat"="MFD inv_nat")) +
  scale_shape_manual(name =" ",values=c("MNNFD_inv"= 1, "MFD_inv_nat"= 2), labels=c("MNNFD_inv"="MNNFD i  ", "MFD_inv_nat"="MFD inv_nat")) +
  theme_bw() +
  theme(legend.position="top") +
  geom_abline(intercept = 0, slope = 0, colour = "grey", size = .5) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Significane of leaf size Distances for\n invasive species to nearest native(MNNFD),\n and native community (MFD)") +
  theme(plot.title=element_text(size=rel(1.5))) +
  theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
dev.off()

pdf("ses.SJ.MNNFD.leafN.FINAL.pdf", width=20, height=10)
ggplot(ses.SJ.MNNFD.leafN, aes(x=reorder(factor(island),as.numeric(as.character(Area.m2))), 
                               y=as.numeric(as.character((obs.z))), color=factor(sigLN), shape=Metric, fill=factor(sigLN))) +
  geom_point(size=10) + #, aes(shape = as.numeric(as.character(p.value.ranks >= 0.05)))) +  
  coord_cartesian(ylim=c(-5, 5)) + 
  scale_y_continuous("standardized effect size", breaks=seq(-8, 8, 2)) +
  scale_x_discrete("Island (increasing size)") +
  scale_fill_manual(values=alpha(c("FALSE"="white", "TRUE"= "black"), .3), legend=F) +
  scale_color_manual(name =" ",values=c("FALSE"= "grey", "TRUE"="black"), legend=F) +
  scale_shape_manual(name =" ",values=c("MNNFD_inv"= 21, "MFD_inv_nat"= 24), labels=c("MNNFD_inv"="MNNFD i  ", "MFD_inv_nat"="MFD_inv_nat")) +
  theme_bw() +
  theme(legend.position="top") +
  geom_abline(intercept = 0, slope = 0, colour = "grey", size = .5) +
  geom_vline(xintercept = c(4.5, 31.5)) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Significane of Leaf Nitrogen difference for\n invasive species to nearest native (NNFDi),\n and native community (MFDNi)\n") +
  theme(plot.title=element_text(size=rel(1.5))) +
  theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
dev.off()

#### summarize significance data by size categories
ses.SJ.MNNFD.leafN$Size.cat <- as.character(ses.SJ.MNNFD.leafN$Size.cat)
ses.SJ.MNNFD.leafN$Size.cat <- factor(ses.SJ.MNNFD.leafN$Size.cat, levels=c("sm", "med", "lg"))

pdf("ses.MNNFD.leafN.SummaryBar.pdf")
ggplot(ses.SJ.MNNFD.leafN, aes(x=Size.cat, fill=factor(Significance))) + 
  geom_bar(stat="bin") +
  facet_wrap(~ Metric) +
  scale_fill_manual(name="Significance",
                    breaks=c("1", "2", "3", "4"),
                    labels=c("Significantly Clustered", "Significantly Overdispersed", 
                             "NS Clustered", "NS Overdispersed"),
                    values=c("1"="black", "2"="orange", "3"="skyblue", "4"="springgreen4")) +
  theme_bw() +
  scale_x_discrete(name="") +
  #theme(legend.position="top") +
  ggtitle("Significance of Leaf N difference for\n invasive species to nearest native (NNFD_inv),\n and native community (MFDN_inv)\n") 
dev.off()

##### ses MMNPDi 
MNNFD_inv <- na.omit(ses.SJ.MNNFD.leafN[which(ses.SJ.MNNFD.leafN$Metric == "MNNFD_inv"), ] )
nrow(MNNFD_inv[which(as.numeric(as.character(MNNFD_inv$p.value.ranks)) <= 0.05), ]) #4
nrow(MNNFD_inv[which(as.numeric(as.character(MNNFD_inv$p.value.ranks)) >= 0.05), ]) #44 NS
4/48
#### ses MPDin 
MFD_inv_nat <- ses.SJ.MNNFD.leafN[which(ses.SJ.MNNFD.leafN$Metric == "MFD_inv_nat"), ] 
nrow(MFD_inv_nat[which(as.numeric(as.character(MFD_inv_nat$p.value.ranks)) <= 0.05), ]) #11
nrow(MFD_inv_nat[which(as.numeric(as.character(MFD_inv_nat$p.value.ranks)) >= 0.05), ]) #37 NS
11/48

sigLeafN_MNNFD_inv <- (MNNFD_inv[which(as.numeric(as.character(MNNFD_inv$p.value.ranks)) <= 0.05), ]) #45 NS
sigLeafN_MPD_inv <- (MFD_inv_nat[which(as.numeric(as.character(MFD_inv_nat$p.value.ranks)) <= 0.05), ]) #45 NS

length(which(sigLeafN_MNNFD_inv$Size.cat == "sm")) #1/48  0.02083333
length(which(sigLeafN_MNNFD_inv$Size.cat == "med")) #3/48   0.0625
length(which(sigLeafN_MNNFD_inv$Size.cat == "lg")) #0   

length(which(sigLeafN_MPD_inv$Size.cat == "sm")) #1/48  0.02083333
length(which(sigLeafN_MPD_inv$Size.cat == "med")) #6/48   0.125
length(which(sigLeafN_MPD_inv$Size.cat == "lg")) #4/48   0.08333333

####### seed mass
## Null distribution 
sim.null.distrib.SeedMass <- lapply(SJ_islands.sim, function(x) sim.meanMNNFD.MFD(phy=SJfinalTree, com=SJcommNewSim, traits=SJtraitLog, island=x, traitname="seedMass", N = 5))
names(sim.null.distrib.SeedMass) <- SJ_islands.sim
head(sim.null.distrib.SeedMass)
head(sim.null.distrib.SeedMass["Willow_Island"])
#write.csv(sim.null.distrib.Height, file="sim.null.distrib.SeedMass.1000.CSV")
sdf.SeedMass <- read.csv(file="~/Dropbox/Work/TankLab/Projects/SanJuans/Manuscript/Drafts/Figs.v4/sim.null.distrib.SeedMass.1000.csv", as.is=T, row.names=1)
head(sdf.SeedMass[, 1:6])
list.sdf.SeedMass <- list()
for (i in 1:ncol(sdf.SeedMass)){
  newlist <- (cbind(sdf.SeedMass[, c(1:6)]))
  head(newlist)
  colnames(newlist) <- c("n.native.tips", "n.invasive.tips", "meanMNNFDinvasives", "meanMNNFDnatives", "meanMFDinv_nat", "meanMFDnat_nat")
  list.sdf.SeedMass[[i]] <- newlist
  sdf.SeedMass <- sdf.SeedMass[, -c(1:6)]
  
}
names(list.sdf.SeedMass) <- SJ_islands.sim
names(list.sdf.SeedMass[1])
head(list.sdf.SeedMass[[1]])
## SES
ses.SanJuan.MNNFD.MFD.SeedMass <- lapply(SJ_islands.sim, function(x) ses.FunctionDist(phy=SJfinalTree, com=SJcommNewSim, island=x,
                                                                                      simOneIslandOneTrait=list.sdf.SeedMass, outputMNNPD=phyloObs[[x]], traits=SJtraitLog, traitname="seedMass", N=1000))
names(ses.SanJuan.MNNFD.MFD.SeedMass) <- SJ_islands.sim
metadata
listIslands <- ses.SanJuan.MNNFD.MFD.SeedMass
ses.SJ.MNNFD.SeedMass <- data.frame()
for (i in 1:length(listIslands)){ 
  if (length(listIslands[[i]]) != 2){
    newlist <- NULL
  } else {
    tmp1 <- metadata[as.character(names(listIslands[i])), "Area.m2"]
    tmp2 <- metadata[as.character(names(listIslands[i])), "Size.cat"]
    newlist <-cbind(t(listIslands[[i]]), "Area.m2"=tmp1, "Size.cat"=tmp2) 
  }
  ses.SJ.MNNFD.SeedMass <- rbind(ses.SJ.MNNFD.SeedMass, newlist)
}
head(ses.SJ.MNNFD.SeedMass)
dim(ses.SJ.MNNFD.SeedMass)
ses.SJ.MNNFD.SeedMass <- na.omit(ses.SJ.MNNFD.SeedMass)
ses.SJ.MNNFD.SeedMass$p.value.ranks <- as.numeric(as.character(ses.SJ.MNNFD.SeedMass$p.value.ranks)) 
ses.SJ.MNNFD.SeedMass.sig <- subset(ses.SJ.MNNFD.SeedMass, p.value.ranks <= 0.05)

sigSM = (ses.SJ.MNNFD.SeedMass[,"p.value.ranks"] <= 0.05)
ses.SJ.MNNFD.SeedMass <- cbind(ses.SJ.MNNFD.SeedMass, sigSM)

ses.SJ.MNNFD.SeedMass[,"Significance"] <- ifelse(ses.SJ.MNNFD.SeedMass[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.SJ.MNNFD.SeedMass[,"obs.z"])) <= 0, 1,
                                                 ifelse(ses.SJ.MNNFD.SeedMass[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.SJ.MNNFD.SeedMass[,"obs.z"])) >= 0, 2, 
                                                        ifelse(ses.SJ.MNNFD.SeedMass[,"p.value.ranks"] >= 0.05 & as.numeric(as.character(ses.SJ.MNNFD.SeedMass[,"obs.z"])) <= 0, 3, 4)))

pdf("ses.SJ.MNNFD.SeedMass.pdf", width=20, height=10)
ggplot(ses.SJ.MNNFD.SeedMass, aes(x=reorder(factor(island),as.numeric(as.character(Area.m2))), 
                                  y=as.numeric(as.character((obs.z))), color=Metric, shape=Metric, size=Size.cat)) +
  geom_point() + #, aes(shape = as.numeric(as.character(p.value.ranks >= 0.05)))) +  
  geom_point(colour=alpha("red", 0.5), shape="*", data = ses.SJ.MNNFD.SeedMass.sig, size = 10) +
  #geom_point(colour=alpha("blue", 0.5), shape="*", data = ses.SJ.MNNFD.SeedMass.sig.oneminus, size = 10) +
  coord_cartesian(ylim=c(-5, 5)) + scale_y_continuous("standardized effect size", breaks=seq(-8, 8, 2)) +
  scale_x_discrete("Island (increasing size)") +
  scale_color_manual(name =" ", values=c("MNNFD_inv"= "grey", "MFD_inv_nat"="black"), labels=c("MNNFD_inv"= "NNFDi  ", "MFD_inv_nat"="MFDNi")) +
  scale_shape_manual(name =" ", values=c("MNNFD_inv"= 1, "MFD_inv_nat"= 2), labels=c("MNNFD_inv"= "NNFDi  ", "MFD_inv_nat"="MFDNi")) + #circle, triangle =2
  scale_size_manual(values=c("sm"= 2, "med"= 4, "lg"=10), labels=c("sm"= 'Small', "med"= 'Medium', "lg"='Large'), breaks=c("sm", "med", "lg")) +
  theme_bw() +
  theme(legend.position="top") +
  geom_abline(intercept = 0, slope = 0, colour = "grey", size = .5) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Significane of Seed Mass difference for\n invasive species to nearest native (NNFDi),\n and native community (MFDNi)") +
  theme(plot.title=element_text(size=rel(1.5))) +
  theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
dev.off()

ses.SJ.MNNFD.SeedMass.four <- ses.SJ.MNNFD.SeedMass[ses.SJ.MNNFD.SeedMass$island %in% four.islands, ]
ses.SJ.MNNFD.SeedMass.sig.four <- ses.SJ.MNNFD.SeedMass.sig[ses.SJ.MNNFD.SeedMass.sig$island %in% four.islands, ]
pdf("four.islands.ses.SJ.MNNFD.SeedMass.four.pdf")
ggplot(ses.SJ.MNNFD.SeedMass.four, aes(x=reorder(factor(island),as.numeric(as.character(Area.m2))), 
                                       y=as.numeric(as.character((obs.z))), color=Metric, shape=Metric)) +
  geom_point(size=5) + #, aes(shape = as.numeric(as.character(p.value.ranks >= 0.05)))) +  
  geom_point(colour=alpha("red", 0.5), shape="*", data = ses.SJ.MNNFD.SeedMass.sig.four, size = 10) +
  #geom_point(colour=alpha("blue", 0.5), shape="*", data = ses.SJ.MNNFD.height.sig.oneminus, size = 10) +
  coord_cartesian(ylim=c(-5, 5)) + scale_y_continuous("standardized effect size", breaks=seq(-8, 8, 2)) +
  scale_x_discrete("Island (increasing size)") +
  scale_color_manual(name =" ",values=c("MNNFD_inv"= "grey", "MFD_inv_nat"="black"), labels=c("MNNFD_inv"="MNNFD i  ", "MFD_inv_nat"="MFD inv_nat")) +
  scale_shape_manual(name =" ",values=c("MNNFD_inv"= 1, "MFD_inv_nat"= 2), labels=c("MNNFD_inv"="MNNFD i  ", "MFD_inv_nat"="MFD inv_nat")) +
  theme_bw() +
  theme(legend.position="top") +
  geom_abline(intercept = 0, slope = 0, colour = "grey", size = .5) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Significane of Seed Mass Distances for\n invasive species to nearest native(MNNFD),\n and native community (MFD)") +
  theme(plot.title=element_text(size=rel(1.5))) +
  theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
dev.off()

pdf("ses.SJ.MNNFD.seedmass.FINAL.pdf", width=20, height=10)
ggplot(ses.SJ.MNNFD.SeedMass, aes(x=reorder(factor(island),as.numeric(as.character(Area.m2))), 
                                  y=as.numeric(as.character((obs.z))), color=factor(sigSM), shape=Metric, fill=factor(sigSM))) +
  geom_point(size=10) + #, aes(shape = as.numeric(as.character(p.value.ranks >= 0.05)))) +  
  coord_cartesian(ylim=c(-5, 5)) + 
  scale_y_continuous("standardized effect size", breaks=seq(-8, 8, 2)) +
  scale_x_discrete("Island (increasing size)") +
  scale_fill_manual(values=alpha(c("FALSE"="white", "TRUE"= "black"), .3), legend=F) +
  scale_color_manual(name =" ",values=c("FALSE"= "grey", "TRUE"="black"), legend=F) +
  scale_shape_manual(name =" ",values=c("MNNFD_inv"= 21, "MFD_inv_nat"= 24), labels=c("MNNFD_inv"="MNNFD i  ", "MFD_inv_nat"="MFD_inv_nat")) +
  theme_bw() +
  theme(legend.position="top") +
  geom_abline(intercept = 0, slope = 0, colour = "grey", size = .5) +
  geom_vline(xintercept = c(18.5, 54.5)) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Significane of Seed Mass difference for\n invasive species to nearest native (NNFDi),\n and native community (MFDNi)\n") +
  theme(plot.title=element_text(size=rel(1.5))) +
  theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
dev.off()

#### summarize significance data by size categories
ses.SJ.MNNFD.SeedMass$Size.cat <- as.character(ses.SJ.MNNFD.SeedMass$Size.cat)
ses.SJ.MNNFD.SeedMass$Size.cat <- factor(ses.SJ.MNNFD.SeedMass$Size.cat, levels=c("sm", "med", "lg"))

pdf("ses.MNNFD.SeedMass.SummaryBar.pdf")
ggplot(ses.SJ.MNNFD.SeedMass, aes(x=Size.cat, fill=factor(Significance))) + 
  geom_bar(stat="bin") +
  facet_wrap(~ Metric) +
  scale_fill_manual(name="Significance",
                    breaks=c("1", "2", "3", "4"),
                    labels=c("Significantly Clustered", "Significantly Overdispersed", 
                             "NS Clustered", "NS Overdispersed"),
                    values=c("1"="black", "2"="orange", "3"="skyblue", "4"="springgreen4")) +
  theme_bw() +
  scale_x_discrete(name="") +
  #theme(legend.position="top") +
  ggtitle("Significance of Seed Mass difference for\n invasive species to nearest native (NNFD_inv),\n and native community (MFDN_inv)\n") 
dev.off()

##### ses MMNPDi 
MNNFD_inv <- na.omit(ses.SJ.MNNFD.SeedMass[which(ses.SJ.MNNFD.SeedMass$Metric == "MNNFD_inv"), ] )
nrow(MNNFD_inv[which(as.numeric(as.character(MNNFD_inv$p.value.ranks)) <= 0.05), ]) #13
nrow(MNNFD_inv[which(as.numeric(as.character(MNNFD_inv$p.value.ranks)) >= 0.05), ]) #59 NS
13/72 #0.1805556
#### ses MPDin 
MFD_inv_nat <- ses.SJ.MNNFD.SeedMass[which(ses.SJ.MNNFD.SeedMass$Metric == "MFD_inv_nat"), ] 
nrow(MFD_inv_nat[which(as.numeric(as.character(MFD_inv_nat$p.value.ranks)) <= 0.05), ]) #18
nrow(MFD_inv_nat[which(as.numeric(as.character(MFD_inv_nat$p.value.ranks)) >= 0.05), ]) #54 NS
18/72 # 0.25

sigseedmassMNNFD_inv <- (MNNFD_inv[which(as.numeric(as.character(MNNFD_inv$p.value.ranks)) <= 0.05), ]) #45 NS
sigseedmassMPD_inv <- (MFD_inv_nat[which(as.numeric(as.character(MFD_inv_nat$p.value.ranks)) <= 0.05), ]) #45 NS

length(which(sigseedmassMNNFD_inv$Size.cat == "sm")) #4/72  0.05555556
length(which(sigseedmassMNNFD_inv$Size.cat == "med")) #6/72   0.08333333
length(which(sigseedmassMNNFD_inv$Size.cat == "lg")) #3/72   0.04166667

length(which(sigseedmassMPD_inv$Size.cat == "sm")) #6/72  0.08333333
length(which(sigseedmassMPD_inv$Size.cat == "med")) #9/72   0.125
length(which(sigseedmassMPD_inv$Size.cat == "lg")) #3/72   0.04166667




################################## ses FUNCTIONAL summary final #################################################################### 
sesMNNPDseed <- ses.SJ.MNNFD.SeedMass[ses.SJ.MNNFD.SeedMass[,"Metric"] =="MNNFD_inv",]
sesMPDseed <- ses.SJ.MNNFD.SeedMass[ses.SJ.MNNFD.SeedMass[,"Metric"] =="MFD_inv_nat",]

sesMNNPDseedPos = cbind("Island"=rownames(sesMNNPDseed), "Size.cat"=as.character(sesMNNPDseed[,"Size.cat"]), 
                       "Metric"=rep(x = "MNNFD inv > 0", times = nrow(sesMNNPDseed)), 
                       "Significance"=ifelse(sesMNNPDseed[,"p.value.ranks"] <= 0.05 & 
                                               as.numeric(as.character(sesMNNPDseed[,"obs.z"])) >= 0, 1,0))
sesMNNPDseedNeg = cbind("Island"=rownames(sesMNNPDseed), "Size.cat"=as.character(sesMNNPDseed[,"Size.cat"]), 
                       "Metric"=rep(x = "MNNFD inv < 0", times = nrow(sesMNNPDseed)), 
                       "Significance"=ifelse(sesMNNPDseed[,"p.value.ranks"] <= 0.05 & 
                                               as.numeric(as.character(sesMNNPDseed[,"obs.z"])) <= 0, 1,0))
sesMPDseedPos = cbind("Island"=rownames(sesMPDseed), "Size.cat"=as.character(sesMPDseed[,"Size.cat"]), 
                      "Metric"=rep(x = "MFD inv > 0", times = nrow(sesMPDseed)), 
                      "Significance"=ifelse(sesMPDseed[,"p.value.ranks"] <= 0.05 & 
                                              as.numeric(as.character(sesMPDseed[,"obs.z"])) >= 0, 1,0))
sesMPDseedNeg = cbind("Island"=rownames(sesMPDseed), "Size.cat"=as.character(sesMPDseed[,"Size.cat"]), 
                      "Metric"=rep(x = "MFD inv < 0", times = nrow(sesMPDseed)), 
                      "Significance"=ifelse(sesMPDseed[,"p.value.ranks"] <= 0.05 & 
                                              as.numeric(as.character(sesMPDseed[,"obs.z"])) <= 0, 1,0))
sesFunctionalSeed <- as.data.frame(rbind(sesMNNPDseedPos, sesMNNPDseedNeg, sesMPDseedPos, sesMPDseedNeg))

length(which(sesMPDseed[,"p.value.ranks"] <= 0.05 & 
               as.numeric(as.character(sesMPDseed[,"obs.z"])) <= 0)) #17/72 ..0.2361111

neworder <- c("MNNFD inv > 0", "MFD inv > 0", "MNNFD inv < 0", "MFD inv < 0")
sesFunctionalSeed2 <- arrange(transform(sesFunctionalSeed, Metric=factor(Metric,levels=neworder)),Metric)

pdf("ses.Functional.Seed.SummaryBar.pdf")
ggplot(sesFunctionalSeed2, aes(x=Significance, fill=factor(Size.cat))) + 
  geom_bar(stat="bin") +
  scale_fill_manual(name="Size Category",
                    breaks=c("sm", "med", "lg"),
                    values=c("sm"="dodgerblue4", "med"="orangered3", "lg"="gold1")) +
  theme_bw() +
  scale_x_discrete(name="", labels=c("0"="NS", "1"="Significant")) +
  facet_wrap(~Metric, ncol = 2) +
  ggtitle("SES Functional Difference\n Seed Mass \n") 
dev.off()

########
sesMNNPDheight <- ses.SJ.MNNFD.height[ses.SJ.MNNFD.height[,"Metric"] =="MNNFD_inv",]
sesMPDheight <- ses.SJ.MNNFD.height[ses.SJ.MNNFD.height[,"Metric"] =="MFD_inv_nat",]

sesMNNPDheightPos = cbind("Island"=rownames(sesMNNPDheight), "Size.cat"=as.character(sesMNNPDheight[,"Size.cat"]), 
                         "Metric"=rep(x = "MNNFD inv > 0", times = nrow(sesMNNPDheight)), 
                         "Significance"=ifelse(sesMNNPDheight[,"p.value.ranks"] <= 0.05 & 
                                                 as.numeric(as.character(sesMNNPDheight[,"obs.z"])) >= 0, 1,0))
sesMNNPDheightNeg = cbind("Island"=rownames(sesMNNPDheight), "Size.cat"=as.character(sesMNNPDheight[,"Size.cat"]), 
                         "Metric"=rep(x = "MNNFD inv < 0", times = nrow(sesMNNPDheight)), 
                         "Significance"=ifelse(sesMNNPDheight[,"p.value.ranks"] <= 0.05 & 
                                                 as.numeric(as.character(sesMNNPDheight[,"obs.z"])) <= 0, 1,0))
sesMPDheightPos = cbind("Island"=rownames(sesMPDheight), "Size.cat"=as.character(sesMPDheight[,"Size.cat"]), 
                        "Metric"=rep(x = "MFD inv > 0", times = nrow(sesMPDheight)), 
                        "Significance"=ifelse(sesMPDheight[,"p.value.ranks"] <= 0.05 & 
                                                as.numeric(as.character(sesMPDheight[,"obs.z"])) >= 0, 1,0))
sesMPDheightNeg = cbind("Island"=rownames(sesMPDheight), "Size.cat"=as.character(sesMPDheight[,"Size.cat"]), 
                        "Metric"=rep(x = "MFD inv < 0", times = nrow(sesMPDheight)), 
                        "Significance"=ifelse(sesMPDheight[,"p.value.ranks"] <= 0.05 & 
                                                as.numeric(as.character(sesMPDheight[,"obs.z"])) <= 0, 1,0))

length(which(sesMPDheight[,"p.value.ranks"] <= 0.05 & 
               as.numeric(as.character(sesMPDheight[,"obs.z"])) <= 0)) # ...0.1690141

sesFunctionalheight <- as.data.frame(rbind(sesMNNPDheightPos, sesMNNPDheightNeg, sesMPDheightPos, sesMPDheightNeg))

neworder <- c("MNNFD inv > 0", "MFD inv > 0", "MNNFD inv < 0", "MFD inv < 0")
sesFunctionalHeight2 <- arrange(transform(sesFunctionalheight, Metric=factor(Metric,levels=neworder)),Metric)

pdf("ses.Functional.Height.SummaryBar.pdf")
ggplot(sesFunctionalHeight2, aes(x=Significance, fill=factor(Size.cat))) + 
  geom_bar(stat="bin") +
  scale_fill_manual(name="Size Category",
                    breaks=c("sm", "med", "lg"),
                    values=c("sm"="dodgerblue4", "med"="orangered3", "lg"="gold1")) +
  theme_bw() +
  scale_x_discrete(name="", labels=c("0"="NS", "1"="Significant")) +
  facet_wrap(~Metric, ncol = 2) +
  ggtitle("SES Functional Difference\n Height \n") 
dev.off()




#######
sesMNNPDSLA <- ses.SJ.MNNFD.SLA[ses.SJ.MNNFD.SLA[,"Metric"] =="MNNFD_inv",]
sesMPDSLA <- ses.SJ.MNNFD.SLA[ses.SJ.MNNFD.SLA[,"Metric"] =="MFD_inv_nat",]

sesMNNPDSLAPos = cbind("Island"=rownames(sesMNNPDSLA), "Size.cat"=as.character(sesMNNPDSLA[,"Size.cat"]), 
                      "Metric"=rep(x = "MNNFD inv > 0", times = nrow(sesMNNPDSLA)), 
                      "Significance"=ifelse(sesMNNPDSLA[,"p.value.ranks"] <= 0.05 & 
                                              as.numeric(as.character(sesMNNPDSLA[,"obs.z"])) >= 0, 1,0))
sesMNNPDSLANeg = cbind("Island"=rownames(sesMNNPDSLA), "Size.cat"=as.character(sesMNNPDSLA[,"Size.cat"]), 
                      "Metric"=rep(x = "MNNFD inv < 0", times = nrow(sesMNNPDSLA)), 
                      "Significance"=ifelse(sesMNNPDSLA[,"p.value.ranks"] <= 0.05 & 
                                              as.numeric(as.character(sesMPDSLA[,"obs.z"])) <= 0, 1,0))
sesMPDSLAPos = cbind("Island"=rownames(sesMPDSLA), "Size.cat"=as.character(sesMPDSLA[,"Size.cat"]), 
                     "Metric"=rep(x = "MFD inv > 0", times = nrow(sesMPDSLA)), 
                     "Significance"=ifelse(sesMPDSLA[,"p.value.ranks"] <= 0.05 & 
                                             as.numeric(as.character(sesMPDSLA[,"obs.z"])) >= 0, 1,0))
sesMPDSLANeg = cbind("Island"=rownames(sesMPDSLA), "Size.cat"=as.character(sesMPDSLA[,"Size.cat"]), 
                     "Metric"=rep(x = "MFD inv < 0", times = nrow(sesMPDSLA)), 
                     "Significance"=ifelse(sesMPDSLA[,"p.value.ranks"] <= 0.05 & 
                                             as.numeric(as.character(sesMPDSLA[,"obs.z"])) <= 0, 1,0))
sesFunctionalSLA <- as.data.frame(rbind(sesMNNPDSLAPos, sesMNNPDSLANeg, sesMPDSLAPos, sesMPDSLANeg))

neworder <- c("MNNFD inv > 0", "MFD inv > 0", "MNNFD inv < 0", "MFD inv < 0")
sesFunctionalSLA2 <- arrange(transform(sesFunctionalSLA, Metric=factor(Metric,levels=neworder)),Metric)

pdf("ses.Functional.SLA.SummaryBar.pdf")
ggplot(sesFunctionalSLA2, aes(x=Significance, fill=factor(Size.cat))) + 
  geom_bar(stat="bin") +
  scale_fill_manual(name="Size Category",
                    breaks=c("sm", "med", "lg"),
                    values=c("sm"="dodgerblue4", "med"="orangered3", "lg"="gold1")) +
  theme_bw() +
  scale_x_discrete(name="", labels=c("0"="NS", "1"="Significant"), limits=c("0","1")) +
  facet_wrap(~Metric, ncol = 2) +
  ggtitle("SES Functional Difference\n SLA \n") 
dev.off()



#######
sesMNNPDleaflet <- ses.SJ.MNNFD.leaflet[ses.SJ.MNNFD.leaflet[,"Metric"] =="MNNFD_inv",]
sesMPDleaflet <- ses.SJ.MNNFD.leaflet[ses.SJ.MNNFD.leaflet[,"Metric"] =="MFD_inv_nat",]

sesMNNPDleafletPos = cbind("Island"=rownames(sesMNNPDleaflet), "Size.cat"=as.character(sesMNNPDleaflet[,"Size.cat"]), 
                          "Metric"=rep(x = "MNNFD inv > 0", times = nrow(sesMNNPDleaflet)), 
                          "Significance"=ifelse(sesMNNPDleaflet[,"p.value.ranks"] <= 0.05 & 
                                                  as.numeric(as.character(sesMNNPDleaflet[,"obs.z"])) >= 0, 1,0))
sesMNNPDleafletNeg = cbind("Island"=rownames(sesMNNPDleaflet), "Size.cat"=as.character(sesMNNPDleaflet[,"Size.cat"]), 
                          "Metric"=rep(x = "MNNFD inv < 0", times = nrow(sesMNNPDleaflet)), 
                          "Significance"=ifelse(sesMNNPDleaflet[,"p.value.ranks"] <= 0.05 & 
                                                  as.numeric(as.character(sesMPDleaflet[,"obs.z"])) <= 0, 1,0))
sesMPDleafletPos = cbind("Island"=rownames(sesMPDleaflet), "Size.cat"=as.character(sesMPDleaflet[,"Size.cat"]), 
                         "Metric"=rep(x = "MFD inv > 0", times = nrow(sesMPDleaflet)), 
                         "Significance"=ifelse(sesMPDleaflet[,"p.value.ranks"] <= 0.05 & 
                                                 as.numeric(as.character(sesMPDleaflet[,"obs.z"])) >= 0, 1,0))
sesMPDleafletNeg = cbind("Island"=rownames(sesMPDleaflet), "Size.cat"=as.character(sesMPDleaflet[,"Size.cat"]), 
                         "Metric"=rep(x = "MFD inv < 0", times = nrow(sesMPDleaflet)), 
                         "Significance"=ifelse(sesMPDleaflet[,"p.value.ranks"] <= 0.05 & 
                                                 as.numeric(as.character(sesMPDleaflet[,"obs.z"])) <= 0, 1,0))
sesFunctionalleaflet <- as.data.frame(rbind(sesMNNPDleafletPos, sesMNNPDleafletNeg, sesMPDleafletPos, sesMPDleafletNeg))

neworder <- c("MNNFD inv > 0", "MFD inv > 0", "MNNFD inv < 0", "MFD inv < 0")
sesFunctionalleaflet2 <- arrange(transform(sesFunctionalleaflet, Metric=factor(Metric,levels=neworder)),Metric)

pdf("ses.Functional.leaflet.SummaryBar.pdf")
ggplot(sesFunctionalleaflet2, aes(x=Significance, fill=factor(Size.cat))) + 
  geom_bar(stat="bin") +
  scale_fill_manual(name="Size Category",
                    breaks=c("sm", "med", "lg"),
                    values=c("sm"="dodgerblue4", "med"="orangered3", "lg"="gold1")) +
  theme_bw() +
  scale_x_discrete(name="", labels=c("0"="NS", "1"="Significant"), limits=c("0","1")) +
  facet_wrap(~Metric, ncol = 2) +
  ggtitle("SES Functional Difference\n Leaf Size \n") 
dev.off()


#######
ses.SJ.MNNFD.leafN.rem <- ses.SJ.MNNFD.leafN[which(ses.SJ.MNNFD.leafN[,'obs.z'] != Inf),]
sesMNNPDleafN <- ses.SJ.MNNFD.leafN.rem[ses.SJ.MNNFD.leafN.rem[,"Metric"] =="MNNFD_inv",]
sesMPDleafN <- ses.SJ.MNNFD.leafN.rem[ses.SJ.MNNFD.leafN.rem[,"Metric"] =="MFD_inv_nat",]

sesMNNPDleafNPos = cbind("Island"=rownames(sesMNNPDleafN), "Size.cat"=as.character(sesMNNPDleafN[,"Size.cat"]), 
                        "Metric"=rep(x = "MNNFD inv > 0", times = nrow(sesMNNPDleafN)), 
                        "Significance"=ifelse(sesMNNPDleafN[,"p.value.ranks"] <= 0.05 & 
                                                as.numeric(as.character(sesMNNPDleafN[,"obs.z"])) >= 0, 1,0))
sesMNNPDleafNNeg = cbind("Island"=rownames(sesMNNPDleafN), "Size.cat"=as.character(sesMNNPDleafN[,"Size.cat"]), 
                        "Metric"=rep(x = "MNNFD inv < 0", times = nrow(sesMNNPDleafN)), 
                        "Significance"=ifelse(sesMNNPDleafN[,"p.value.ranks"] <= 0.05 & 
                                                as.numeric(as.character(sesMPDleafN[,"obs.z"])) <= 0, 1,0))
sesMPDleafNPos = cbind("Island"=rownames(sesMPDleafN), "Size.cat"=as.character(sesMPDleafN[,"Size.cat"]), 
                       "Metric"=rep(x = "MFD inv > 0", times = nrow(sesMPDleafN)), 
                       "Significance"=ifelse(sesMPDleafN[,"p.value.ranks"] <= 0.05 & 
                                               as.numeric(as.character(sesMPDleafN[,"obs.z"])) >= 0, 1,0))
sesMPDleafNNeg = cbind("Island"=rownames(sesMPDleafN), "Size.cat"=as.character(sesMPDleafN[,"Size.cat"]), 
                       "Metric"=rep(x = "MFD inv < 0", times = nrow(sesMPDleafN)), 
                       "Significance"=ifelse(sesMPDleafN[,"p.value.ranks"] <= 0.05 & 
                                               as.numeric(as.character(sesMPDleafN[,"obs.z"])) <= 0, 1,0))
sesFunctionalleafN <- as.data.frame(rbind(sesMNNPDleafNPos, sesMNNPDleafNNeg, sesMPDleafNPos, sesMPDleafNNeg))


neworder <- c("MNNFD inv > 0", "MFD inv > 0", "MNNFD inv < 0", "MFD inv < 0")
sesFunctionalleafN2 <- arrange(transform(sesFunctionalleafN, Metric=factor(Metric,levels=neworder)),Metric)

pdf("ses.Functional.leafN.SummaryBar.pdf")
ggplot(sesFunctionalleafN2, aes(x=Significance, fill=factor(Size.cat))) + 
  geom_bar(stat="bin") +
  scale_fill_manual(name="Size Category",
                    breaks=c("sm", "med", "lg"),
                    values=c("sm"="dodgerblue4", "med"="orangered3", "lg"="gold1")) +
  theme_bw() +
  scale_x_discrete(name="", labels=c("0"="NS", "1"="Significant"), limits=c("0","1")) +
  facet_wrap(~Metric, ncol = 2) +
  ggtitle("SES Functional Difference\n Leaf N \n") 
dev.off()




