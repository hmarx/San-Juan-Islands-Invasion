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


############################ Null distributions for each trait...this will take a while 

####### seed mass
## Null distribution 
sim.null.distrib.SeedMass <- lapply(SJ_islands.sim, function(x) sim.meanMNNFD.MFD(phy=SJfinalTree, com=SJcommNewSim, traits=SJtraitLog, island=x, traitname="seedMass", N = 1000))
names(sim.null.distrib.SeedMass) <- SJ_islands.sim
head(sim.null.distrib.SeedMass)
head(sim.null.distrib.SeedMass["Willow_Island"])
#write.csv(sim.null.distrib.Height, file="sim.null.distrib.SeedMass.1000.CSV")

####  Height
sim.null.distrib.Height <- lapply(SJ_islands.sim, function(x) sim.meanMNNFD.MFD(phy=SJfinalTree, com=SJcommNewSim, traits=SJtraitLog, island=x, traitname="maxHeight", N = 1000))
names(sim.null.distrib.Height) <- SJ_islands.sim
head(sim.null.distrib.Height)
head(sim.null.distrib.Height["Willow_Island"])
#write.csv(sim.null.distrib.Height, file="sim.null.distrib.Height.null1000.csv")

#### SLA
## Null distribution 
sim.null.distrib.SLA <- lapply(SJ_islands.sim, function(x) sim.meanMNNFD.MFD(phy=SJfinalTree, com=SJcommNewSim, traits=SJtraitLog, island=x, traitname="sla", N = 1000))
names(sim.null.distrib.SLA) <- SJ_islands.sim
#write.csv(sim.null.distrib.SLA, file="sim.null.distrib.SLA.1000.csv")

####### leaflet
## Null distribution 
sim.null.distrib.leafletSize <- lapply(SJ_islands.sim, function(x) sim.meanMNNFD.MFD(phy=SJfinalTree, com=SJcommNewSim, traits=SJtraitLog, island=x, traitname="leafletSize", N = 1000))
names(sim.null.distrib.leafletSize) <- SJ_islands.sim
head(sim.null.distrib.leafletSize)
head(sim.null.distrib.leafletSize["Willow_Island"])
#write.csv(sim.null.distrib.leafletSize, file="sim.null.distrib.leafletSize.1000.csv")

####### leaf N
## Null distribution 
sim.null.distrib.leafN <- lapply(SJ_islands.sim, function(x) sim.meanMNNFD.MFD(phy=SJfinalTree, com=SJcommNewSim, traits=SJtraitLog, island=x, traitname="leafN", N = 1000))
names(sim.null.distrib.leafN) <- SJ_islands.sim
head(sim.null.distrib.leafN)
head(sim.null.distrib.leafN["Willow_Island"])
#write.csv(sim.null.distrib.Height, file="sim.null.distrib.leafN.1000.csv")


##################### Summarize results with plots

#### Seed mass
sum.sesFunctionDist(plottype = "NullInvOcc", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.SeedMass.1000.CSV", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="seedMass", metadata=metadata)

sum.sesFunctionDist(plottype = "ses.allIslands", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.SeedMass.1000.CSV", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="seedMass", metadata=metadata)

sum.sesFunctionDist(plottype = "summary.Bar", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.SeedMass.1000.CSV", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="seedMass", metadata=metadata)

#### Max Height
sum.sesFunctionDist(plottype = "NullInvOcc", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.Height.null1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="maxHeight", metadata=metadata)

sum.sesFunctionDist(plottype = "ses.allIslands", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.Height.null1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="maxHeight", metadata=metadata)

sum.sesFunctionDist(plottype = "summary.Bar", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.Height.null1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="maxHeight", metadata=metadata)


#### sla
sum.sesFunctionDist(plottype = "NullInvOcc", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.SLA.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="sla", metadata=metadata)

sum.sesFunctionDist(plottype = "ses.allIslands", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.SLA.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="sla", metadata=metadata)

sum.sesFunctionDist(plottype = "summary.Bar", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.SLA.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="sla", metadata=metadata)


#### leaf size
sum.sesFunctionDist(plottype = "NullInvOcc", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.leafletSize.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="leafletSize", metadata=metadata)

sum.sesFunctionDist(plottype = "ses.allIslands", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.leafletSize.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="leafletSize", metadata=metadata)

sum.sesFunctionDist(plottype = "summary.Bar", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.leafletSize.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="leafletSize", metadata=metadata)

#### leaf N
sum.sesFunctionDist(plottype = "NullInvOcc", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.leafN.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="leafN", metadata=metadata)

sum.sesFunctionDist(plottype = "ses.allIslands", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.leafN.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="leafN", metadata=metadata)

sum.sesFunctionDist(plottype = "summary.Bar", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.leafN.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="leafN", metadata=metadata)




