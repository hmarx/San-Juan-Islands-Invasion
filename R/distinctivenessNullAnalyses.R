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

## Observed difference in mean DNNS i / MDNS inv-nat
obs.DNNS.SJ  <- lapply(SJ_islands, function(x) phyloDistinct(phy=SJfinalTree, community=SJcommNew, col=x))
names(obs.DNNS.SJ) <- SJ_islands # = phyloObs
head(obs.DNNS.SJ)
summ.DNNS.SJ  <- lapply(SJ_islands, function(x) summary.DNNS.MDNS(obs.DNNS.SJ[[x]])) #apply summary funciton across all communities
names(summ.DNNS.SJ) <- SJ_islands


#### Summarize observed mean compared to simualted means == standardized effect size 
## p.rank.DNNS.inv <- min(DNNS.inv.rankLo, DNNS.inv.rankHi) / (N + 1) 
## proportion of simulated means as or more extreme than observed
ses.SanJuan.DNNS.MDNS <- lapply(names(list.sdf), function(x) ses.PhyloDist(phy=SJfinalTree, com=SJcommNew, island=x, simOneIsland=list.sdf, N=1000))
names(ses.SanJuan.DNNS.MDNS) <- names(list.sdf)
length(ses.SanJuan.DNNS.MDNS)


#pdf("figs/plots/phyloDiv/ses/logNullInvOcc.DNNS.islIncSize.pdf", width=20, height=10)
sum.sesPhyloDist(plottype = "NullObsIntervalDNNS", simPhyloOut = "output/10_Analyses/PhylogeneticDiversity/Null/SanJuan.DNNS.MDNS.null1000.csv",
                 metadata = metadata)
#dev.off()

#pdf("figs/plots/phyloDiv/ses/logNullInvOcc.MDNS.islIncSize.pdf", width=20, height=10)
sum.sesPhyloDist(plottype = "NullObsIntervalMDNS", simPhyloOut = "output/10_Analyses/PhylogeneticDiversity/Null/SanJuan.DNNS.MDNS.null1000.csv",
                 metadata = metadata)
#dev.off()

#pdf("figs/plots/phyloDiv/ses/ses.DNNS.MDNS.IncSize.Final.pdf", width=20, height=10)
sum.sesPhyloDist(plottype = "ses", simPhyloOut = "output/10_Analyses/PhylogeneticDiversity/Null/SanJuan.DNNS.MDNS.null1000.csv",
                 metadata = metadata)
#dev.off()

#pdf("figs/plots/phyloDiv/ses/ses.phylo.SummaryBar.pdf")
sum.sesPhyloDist(plottype = "summary.Bar", simPhyloOut = "output/10_Analyses/PhylogeneticDiversity/Null/SanJuan.DNNS.MDNS.null1000.csv",
                 metadata = metadata)
#dev.off()


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

########### Seed mass ########### 
#pdf(file="figs/plots/functionDiv/ses/NullInvOcc.NNFD.seedMass")
sum.sesFunctionDist(plottype = "NullObsIntervalNNFD", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.SeedMass.1000.CSV", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="seedMass", metadata=metadata)
#dev.off()

#pdf(paste("figs/plots/functionDiv/ses/NullInvOcc.MFD", traitname, "pdf", sep="."), width=20, height=10)
sum.sesFunctionDist(plottype = "NullObsIntervalMFD", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.SeedMass.1000.CSV", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="seedMass", metadata=metadata)
#dev.off()

#pdf(paste("figs/plots/functionDiv/ses/ses.SJ.NNFD", traitname, "pdf", sep="."), width=20, height=10)
sum.sesFunctionDist(plottype = "ses", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.SeedMass.1000.CSV", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="seedMass", metadata=metadata)
#dev.off()

#pdf(paste("figs/plots/functionDiv/ses/", traitname, "Functional.SummaryBar.pdf", sep=""))
sum.sesFunctionDist(plottype = "summary.Bar", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.SeedMass.1000.CSV", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="seedMass", metadata=metadata)
#dev.off()

########### Max Height ########### 
sum.sesFunctionDist(plottype = "NullObsIntervalNNFD", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.Height.null1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="maxHeight", metadata=metadata)

sum.sesFunctionDist(plottype = "NullObsIntervalMFD", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.Height.null1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="maxHeight", metadata=metadata)

sum.sesFunctionDist(plottype = "ses", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.Height.null1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="maxHeight", metadata=metadata)

sum.sesFunctionDist(plottype = "summary.Bar", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.Height.null1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="maxHeight", metadata=metadata)


########### sla ########### 
sum.sesFunctionDist(plottype = "NullObsIntervalNNFD", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.SLA.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="sla", metadata=metadata)

sum.sesFunctionDist(plottype = "NullObsIntervalMFD", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.SLA.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="sla", metadata=metadata)

sum.sesFunctionDist(plottype = "ses", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.SLA.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="sla", metadata=metadata)

sum.sesFunctionDist(plottype = "summary.Bar", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.SLA.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="sla", metadata=metadata)


########### leaf size ########### 
sum.sesFunctionDist(plottype = "NullObsIntervalNNFD", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.leafletSize.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="leafletSize", metadata=metadata)

sum.sesFunctionDist(plottype = "NullObsIntervalMFD", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.leafletSize.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="leafletSize", metadata=metadata)

sum.sesFunctionDist(plottype = "ses", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.leafletSize.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="leafletSize", metadata=metadata)

sum.sesFunctionDist(plottype = "summary.Bar", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.leafletSize.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="leafletSize", metadata=metadata)

########### leaf N ########### 
sum.sesFunctionDist(plottype = "NullObsIntervalNNFD", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.leafN.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="leafN", metadata=metadata)

sum.sesFunctionDist(plottype = "NullObsIntervalMFD", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.leafN.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="leafN", metadata=metadata)

sum.sesFunctionDist(plottype = "ses", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.leafN.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="leafN", metadata=metadata)

sum.sesFunctionDist(plottype = "summary.Bar", sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.leafN.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="leafN", metadata=metadata)




