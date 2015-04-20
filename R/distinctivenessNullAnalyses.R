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

################# RANDOMIZATION == shuffle prese/abs of invasive within each community (ses.DNNS.MDNS.R) ##############
################# Compare observed MMNPDi / MDNSinv-nat to random distribution of invasive spevies in each community
## mean obs MMNPDi < mean null MMNPDi == intorduced species is closer to native relative than a random invasive from species pool
## mean obs MMNPDi = mean null MMNPDi == any introduced species would show the same pattern
## mean obs MMNPDi > mean null MMNPDi == intorduced species is more distant to native relative than random; DNH

################################ DNNS i / MDNS inv-nat ################################ 
######################## Simulate means, compare null distribution to the observed mean 

## Simulate DNNS/MDNS by randomizing invasive communities...NOTE: will take a while to run
SJ_islands <- names(SJcommNew) # names of the islands 
sim.null.distrib.SJ <- lapply(SJ_islands, function(x) sim.meanDNNS.MDNS(phy=SJfinalTree, com=SJcommNew, traits=SJtraitLog, island=x, N = 1000))
names(sim.null.distrib.SJ) <- SJ_islands
head(sim.null.distrib.SJ)
head(sim.null.distrib.SJ["All_SanJuanIslands"])
#write.csv(sim.null.distrib.SJ, file="SanJuan.DNNS.MDNS.null1000.csv")

## Read in simulation file
sdf <- read.csv(file="output/10_Analyses/PhylogeneticDiversity/Null/SanJuan.DNNS.MDNS.null1000.csv", as.is=T, row.names=1)
head(sdf[, 1:6])
list.sdf <- list()
for (i in 1:ncol(sdf)){
  newlist <- (cbind(sdf[, c(1:6)]))
  head(newlist)
  colnames(newlist) <- c("n.native.tips", "n.invasive.tips", "meanDNNSinvasives", "meanDNNSnatives", "meanMDNSinv_nat", "meanMDNSnat_nat")
  list.sdf[[i]] <- newlist
  sdf <- sdf[, -c(1:6)]
  
}
names(list.sdf) <- SJ_islands

## Observed difference in mean DNNS i / MDNS inv-nat
obs.DNNS.SJ  <- lapply(SJ_islands, function(x) phyloDistinct(phy=SJfinalTree, community=SJcommNew, col=x))
names(obs.DNNS.SJ) <- SJ_islands # = phyloObs
head(obs.DNNS.SJ)
summ.DNNS.SJ  <- lapply(SJ_islands, function(x) summary.DNNS.MDNS(obs.DNNS.SJ[[x]])) #apply summary funciton across all communities
names(summ.DNNS.SJ) <- SJ_islands

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
  tmp.DNNS.i <- rep(summ.DNNS.SJ[as.character(names(summ.DNNS.SJ[i]))][[1]][[2]], times=nrow(listIslands[[i]])) # observed DNNS.i
  tmp.MDNS.in <- rep(summ.DNNS.SJ[as.character(names(summ.DNNS.SJ[i]))][[1]][[9]], times=nrow(listIslands[[i]])) # observed MDNS.i to mean native community
  newlist <- mapply(cbind, "islands"=names(listIslands[i]), listIslands[i], "Area.m2"=tmp, "Obs.meanDNNSinvasives"=tmp.DNNS.i, "Obs.meanMDNSinv_nat"= tmp.MDNS.in, SIMPLIFY=F)   
  list.meta.null.distrib.SJ[i] <- newlist[1]
}
names(list.meta.null.distrib.SJ) <- names(list.sdf.new)
head(list.meta.null.distrib.SJ[[1]])

## Melt list of simulated mean DNNS/MDNS and observed means into dataframe, and remove NA
sim.null.distrib.melt <- melt(list.meta.null.distrib.SJ, measure.vars="islands")
sim.null.distrib.melt <- sim.null.distrib.melt[which(!sim.null.distrib.melt$value == "NA"),] 
head(sim.null.distrib.melt)

## Plot distribution of null means, with observed mean, for each island 
pdf("figs/plots/phyloDiv/ses/logNullInvOcc.DNNS.islIncSize.pdf", width=20, height=10)
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
p
dev.off()

pdf("figs/plots/phyloDiv/ses/logNullInvOcc.MDNS.islIncSize.pdf", width=20, height=10)
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
p
dev.off()


#### Summarize observed mean compared to simualted means == standardized effect size 
## p.rank.DNNS.inv <- min(DNNS.inv.rankLo, DNNS.inv.rankHi) / (N + 1) 
## proportion of simulated means as or more extreme than observed
ses.SanJuan.DNNS.MDNS <- lapply(names(list.sdf.new), function(x) ses.PhyloDist(phy=SJfinalTree, com=SJcommNew, island=x, simOneIsland=list.sdf, N=1000))
names(ses.SanJuan.DNNS.MDNS) <- names(list.sdf.new)
length(ses.SanJuan.DNNS.MDNS)

head(metadata)
listIslands <- ses.SanJuan.DNNS.MDNS
ses.DNNS.MDNS <- data.frame()
for (i in 1:length(listIslands)){ 
  tmp1 <- metadata[as.character(names(listIslands[i])), "Area.m2"]
  tmp2 <- metadata[as.character(names(listIslands[i])), "Size.cat"]
  newlist <-cbind(t(listIslands[[i]]), "Area.m2"=tmp1, "Size.cat"=tmp2) 
  ses.DNNS.MDNS <- rbind(ses.DNNS.MDNS, newlist)
}
head(ses.DNNS.MDNS)
#write.csv(ses.DNNS.MDNS, file="ses.DNNS.MDNS.MDNSNew2.csv")
ses.DNNS.MDNS$p.value.ranks <- as.numeric(as.character(ses.DNNS.MDNS$p.value.ranks)) 
ses.DNNS.MDNS <- na.omit(ses.DNNS.MDNS)

## add a column indicating significant islands for plotting
sig = (ses.DNNS.MDNS[,"p.value.ranks"] <= 0.05)
neg = (as.numeric(as.character(ses.DNNS.MDNS[,"obs.z"])) <= 0)
ses.DNNS.MDNS <- cbind(ses.DNNS.MDNS, sig)

## code for signficance of observed means, and for direction of pattern (clusering / overdispersed)
ses.DNNS.MDNS[,"Significance"] <- ifelse(ses.DNNS.MDNS[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.DNNS.MDNS[,"obs.z"])) <= 0, 1,
                                        ifelse(ses.DNNS.MDNS[,"p.value.ranks"] <= 0.05 & as.numeric(as.character(ses.DNNS.MDNS[,"obs.z"])) >= 0, 2, 
                                               ifelse(ses.DNNS.MDNS[,"p.value.ranks"] >= 0.05 & as.numeric(as.character(ses.DNNS.MDNS[,"obs.z"])) <= 0, 3, 4)))

pdf("figs/plots/phyloDiv/ses/ses.DNNS.MDNS.IncSize.Final.pdf", width=20, height=10)
ggplot(ses.DNNS.MDNS, aes(x=reorder(factor(island),as.numeric(as.character(Area.m2))), 
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
dev.off()

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
DNNS_inv <- ses.DNNS.MDNS[which(ses.DNNS.MDNS$Metric == "DNNS_inv"), ] 
nrow(DNNS_inv[which(as.numeric(as.character(DNNS_inv$p.value.ranks)) <= 0.05), ]) # 22
nrow(DNNS_inv[which(as.numeric(as.character(DNNS_inv$p.value.ranks)) >= 0.05), ]) #50 NS
22/72  # 0.3055556

sigDNNS_inv <- (DNNS_inv[which(as.numeric(as.character(DNNS_inv$p.value.ranks)) <= 0.05), ])
nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigDNNS_inv, by=0)) # 10/72 small islands sig 0.1388889
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigDNNS_inv, by=0)) # 6/72 med islands sig  0.08333333
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigDNNS_inv, by=0)) # 6/72 large isl sig  0.08333333

## ses MDNSin 
MDNS_inv_nat <- ses.DNNS.MDNS[which(ses.DNNS.MDNS$Metric == "MDNS_inv_nat"), ] 
nrow(MDNS_inv_nat[which(as.numeric(as.character(MDNS_inv_nat$p.value.ranks)) <= 0.05), ]) #27 significant
nrow(MDNS_inv_nat[which(as.numeric(as.character(MDNS_inv_nat$p.value.ranks)) >= 0.05), ]) #45 NS
27/72  # 0.375

length(which(as.numeric(as.character(MDNS_inv_nat$p.value.ranks)) <= 0.05 &
               as.numeric(as.character(MDNS_inv_nat[,"obs.z"])) < 0 )) ## 5 MDNS clustered, significant
length(which(as.numeric(as.character(MDNS_inv_nat$p.value.ranks)) <= 0.05 &
               as.numeric(as.character(MDNS_inv_nat[,"obs.z"])) > 0 )) ## 22 MDNS overdispersed, significant
5/72 #0.06944444
22/72 #0.3055556

sigMDNS_inv <- (MDNS_inv_nat[which(as.numeric(as.character(MDNS_inv_nat$p.value.ranks)) <= 0.05), ]) #45 NS
nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigMDNS_inv, by=0)) # 8/72 small islands sig 0.1111111
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigMDNS_inv, by=0)) # 30/72 med islands sig  0.4166667
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigMDNS_inv, by=0)) # 12/72 large isl sig  0.1666667
length(which(sigMDNS_inv$Size.cat == "sm")) #7/72  0.09722222
length(which(sigMDNS_inv$Size.cat == "med")) #16/72   0.2222222
length(which(sigMDNS_inv$Size.cat == "lg")) #4/72   0.05555556



###########################################  Significance of Functional Distances  #############################################################
remove.islands.sim <- c("All_SanJuanIslands", "Unnamed_west_of_Castle_Island") 
SJcommNewSim  <- SJcommNew[, -which(names(SJcommNew) %in% remove.islands.sim)]
head(SJcommNewSim)
SJ_islands.sim <- colnames(SJcommNewSim)


############################ Null distributions for each trait...this will take a while 

####### seed mass
## Null distribution 
sim.null.distrib.SeedMass <- lapply(SJ_islands.sim, function(x) sim.meanNNFD.MFD(phy=SJfinalTree, com=SJcommNewSim, traits=SJtraitLog, island=x, traitname="seedMass", N = 1000))
names(sim.null.distrib.SeedMass) <- SJ_islands.sim
head(sim.null.distrib.SeedMass)
head(sim.null.distrib.SeedMass["Willow_Island"])
#write.csv(sim.null.distrib.SeedMass, file="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.SeedMass.1000.CSV")

####  Height
sim.null.distrib.Height <- lapply(SJ_islands.sim, function(x) sim.meanNNFD.MFD(phy=SJfinalTree, com=SJcommNewSim, traits=SJtraitLog, island=x, traitname="maxHeight", N = 1000))
names(sim.null.distrib.Height) <- SJ_islands.sim
head(sim.null.distrib.Height)
head(sim.null.distrib.Height["Willow_Island"])
#write.csv(sim.null.distrib.Height, file="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.Height.null1000.csv")

#### SLA
## Null distribution 
sim.null.distrib.SLA <- lapply(SJ_islands.sim, function(x) sim.meanNNFD.MFD(phy=SJfinalTree, com=SJcommNewSim, traits=SJtraitLog, island=x, traitname="sla", N = 1000))
names(sim.null.distrib.SLA) <- SJ_islands.sim
#write.csv(sim.null.distrib.SLA, file="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.SLA.1000.csv")

####### leaflet
## Null distribution 
sim.null.distrib.leafletSize <- lapply(SJ_islands.sim, function(x) sim.meanNNFD.MFD(phy=SJfinalTree, com=SJcommNewSim, traits=SJtraitLog, island=x, traitname="leafletSize", N = 1000))
names(sim.null.distrib.leafletSize) <- SJ_islands.sim
head(sim.null.distrib.leafletSize)
head(sim.null.distrib.leafletSize["Willow_Island"])
#write.csv(sim.null.distrib.leafletSize, file="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.leafletSize.1000.csv")

####### leaf N
## Null distribution 
sim.null.distrib.leafN <- lapply(SJ_islands.sim, function(x) sim.meanNNFD.MFD(phy=SJfinalTree, com=SJcommNewSim, traits=SJtraitLog, island=x, traitname="leafN", N = 1000))
names(sim.null.distrib.leafN) <- SJ_islands.sim
head(sim.null.distrib.leafN)
head(sim.null.distrib.leafN["Willow_Island"])
#write.csv(sim.null.distrib.leafN, file="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.leafN.1000.csv")


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




