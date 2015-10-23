 
#### New plots of observed and expected phylogenetic and functional distinctiveness ###
#### Following reviewer's suggestions ####
#### 4 September 2015 ####
#### Hannah E. Marx ####

source("R/DistinctivenessFunctions.R")

################################ Phylogenetic Distinctiveness  ################################ 
########  Observed
head(phyloObsSum)

phyloObsSum.DNNS = cbind("Island"=phyloObsSum[,1], "Size.cat"=as.character(phyloObsSum[,"Size.cat"]), 
                         "Metric"=rep(x = "DNNS", times = nrow(phyloObsSum)), "p.value"=round(as.numeric(phyloObsSum[,"t.DNNS.p.value"]), 
                                                                        digits = 2), "Area" = phyloObsSum$Area.m2)
phyloObsSum.MDNS = cbind("Island"=phyloObsSum[,1], "Size.cat"=as.character(phyloObsSum[,"Size.cat"]), 
                         "Metric"=rep(x = "MDNS", times = nrow(phyloObsSum)), "p.value"=round(as.numeric(phyloObsSum[,"t.MDNS.p.value"]), digits=2), "Area" =phyloObsSum$Area.m2)
phyloObsSum.DNNS.MDNS <- as.data.frame(rbind(phyloObsSum.DNNS, phyloObsSum.MDNS))
head(phyloObsSum.DNNS.MDNS)
phyloObsSum.DNNS.MDNS <- phyloObsSum.DNNS.MDNS[phyloObsSum.DNNS.MDNS$Island != 'Unnamed_west_of_Castle_Island',]
phyloObsSum.DNNS.MDNS$p.value <- as.numeric(as.character(phyloObsSum.DNNS.MDNS$p.value))
phyloObsSum.DNNS.MDNS$Island <- gsub(phyloObsSum.DNNS.MDNS$Island, pattern = "___", replacement = " - ")
phyloObsSum.DNNS.MDNS$Island <- gsub(phyloObsSum.DNNS.MDNS$Island, pattern = "__", replacement = "_#")
phyloObsSum.DNNS.MDNS$Island <- gsub(phyloObsSum.DNNS.MDNS$Island, pattern = "_", replacement = " ")
phyloObsSum.DNNS.MDNS$sig <- ifelse(phyloObsSum.DNNS.MDNS$p.value <= 0.05, 1, 0)


phyloObsSum.DNNS.MDNS$Y1 <- cut(phyloObsSum.DNNS.MDNS$p.value,breaks = c(0:0.049, 0.05, 1) ,right = FALSE)

p <- ggplot(phyloObsSum.DNNS.MDNS, aes(y=reorder(factor(Island), as.numeric(as.character(Area))), x=Metric))#, col=c("magenta1", "green3"))
p <- p + geom_tile(aes(fill = p.value), colour="white")
p <- p + scale_fill_gradient(low="steelblue", high="white", na.value="black")
base_size <- 9
p <- p + theme_grey(base_size = base_size) + labs(x = "",  y = "") 
p <- p + scale_x_discrete(expand = c(0, 0)) 
p <- p + scale_y_discrete(expand = c(0, 0)) 
p <- p + theme(legend.position = "none", axis.ticks = element_blank(), 
       axis.text.x = element_text(size = base_size * 0.8, angle = -45, hjust = 0, colour = "grey50"))
p <- p + coord_fixed(ratio=1)
p

obs.plot <- ggplot(phyloObsSum.DNNS.MDNS, aes(y=reorder(factor(Island), as.numeric(as.character(Area))), x=Metric))#, col=c("magenta1", "green3"))
obs.plot <- obs.plot + geom_tile(aes(fill = p.value), colour="white")
obs.plot <- obs.plot + scale_fill_gradient(low="steelblue", high="white", na.value="black")
obs.plot <- obs.plot + geom_point(aes(size=ifelse(sig, "dot", "no_dot")))
obs.plot <- obs.plot + scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")
base_size <- 9
obs.plot <- obs.plot + theme_grey(base_size = base_size) + labs(x = "",  y = "") 
obs.plot <- obs.plot + scale_x_discrete(expand = c(0, 0)) 
obs.plot <- obs.plot + scale_y_discrete(expand = c(0, 0)) 
obs.plot <- obs.plot + theme(legend.position = "none", axis.ticks = element_blank(), 
               axis.text.x = element_text(size = base_size * 0.8, angle = -45, hjust = 0, colour = "black"))
obs.plot <- obs.plot + coord_fixed(ratio=1)
#pdf(file="figs/plots/phyloDiv/observed/obs.phylo.tile.pdf")
obs.plot
#dev.off()


########  Expected 
head(ses.DNNS.MDNS)
phylo.exp.DNNS.MDNS <- ses.DNNS.MDNS
phylo.exp.DNNS.MDNS$obs.z <- as.numeric(as.character(phylo.exp.DNNS.MDNS$obs.z))
phylo.exp.DNNS.MDNS$island <- gsub(phylo.exp.DNNS.MDNS$island, pattern = "___", replacement = " - ")
phylo.exp.DNNS.MDNS$island <- gsub(phylo.exp.DNNS.MDNS$island, pattern = "__", replacement = "_#")
phylo.exp.DNNS.MDNS$island <- gsub(phylo.exp.DNNS.MDNS$island, pattern = "_", replacement = " ")

p <- ggplot(phylo.exp.DNNS.MDNS, aes(y=reorder(factor(island), as.numeric(as.character(Area.m2))), x=Metric))#, col=c("magenta1", "green3"))
p <- p + geom_tile(aes(fill = obs.z), colour="white")
p <- p + scale_fill_gradient(low="darkblue", high="darkorange", na.value="black")
base_size <- 9
p <- p + theme_grey(base_size = base_size) + labs(x = "",  y = "") 
p <- p + scale_x_discrete(expand = c(0, 0)) 
p <- p + scale_y_discrete(expand = c(0, 0)) 
p <- p + theme(legend.position = "none", axis.ticks = element_blank(), 
               axis.text.x = element_text(size = base_size * 0.8, angle = -45, hjust = 0, colour = "grey50"))
p <- p + coord_fixed(ratio=1)
p

ses.plot <- ggplot(phylo.exp.DNNS.MDNS, aes(y=reorder(factor(island), as.numeric(as.character(Area.m2))), x=Metric, fill=obs.z))#, col=c("magenta1", "green3"))
ses.plot <- ses.plot + geom_tile(colour="white")
ses.plot <- ses.plot + scale_fill_gradient(low="darkblue", high="darkorange", na.value="black")
ses.plot <- ses.plot + geom_point(aes(size=ifelse(sig, "dot", "no_dot")))
ses.plot <- ses.plot + scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")
base_size <- 9
ses.plot <- ses.plot + theme_grey(base_size = base_size) + labs(x = "",  y = "") 
ses.plot <- ses.plot + scale_x_discrete(expand = c(0, 0)) 
ses.plot <- ses.plot + scale_y_discrete(expand = c(0, 0)) 
ses.plot <- ses.plot + theme(legend.position = "none", axis.ticks = element_blank(), 
               axis.text.x = element_text(size = base_size * 0.8, angle = -45, hjust = 0, colour = "black"))
ses.plot <- ses.plot + coord_fixed(ratio=1)
#pdf(file="figs/plots/phyloDiv/ses/ses.phylo.tile.pdf")
ses.plot
#dev.off()

################################ Functional Distinctiveness################################ 
sigSeed = cbind("Island"=seedmass[,"Row.names"], "Size.cat"= as.character(seedmass[,"Size.cat"]), "Metric"=rep(x = "Seed Mass", times = nrow(seedmass)), "p.value" = seedmass[,"p.value.seedMass"], "Area"=as.character(seedmass[,"Area.m2"]),  "sig"=ifelse(seedmass[,"p.value.seedMass"] <= 0.05, 1,0))
sigNNFDSeed = cbind("Island"=seedmass[,"Row.names"],"Size.cat"= as.character(seedmass[,"Size.cat"]), "Metric"=rep(x = "NNFD Seed Mass", times = nrow(seedmass)),"p.value" = seedmass[,"t.NNFD.p.value"], "Area"=as.character(seedmass[,"Area.m2"]),  "sig"=ifelse(seedmass[,"t.NNFD.p.value"] <= 0.05, 1,0))
sigMFDSeed = cbind("Island"=seedmass[,"Row.names"],"Size.cat"=as.character(seedmass[,"Size.cat"]), "Metric"=rep(x = "MFD Seed Mass", times = nrow(seedmass)), "p.value" = seedmass[,"t.NNFD.p.value"], "Area"=as.character(seedmass[,"Area.m2"]),  "sig"=ifelse(seedmass[,"t.MFD.p.value"] <= 0.05, 1,0))

sigHeight = cbind("Island"=height[,"Row.names"],"Size.cat"=as.character(seedmass[,"Size.cat"]), "Metric"=rep(x = "Height", times = nrow(height)), "p.value" = height[,"p.value.maxHeight"], "Area"=as.character(height[,"Area.m2"]),  "sig"=ifelse(height[,"p.value.maxHeight"] <= 0.05, 1,0))
sigNNFDheight = cbind("Island"=height[,"Row.names"],"Size.cat"=as.character(seedmass[,"Size.cat"]),"Metric"=rep(x = "NNFD Height", times = nrow(height)), "p.value" = height[,"t.NNFD.p.value"], "Area"=as.character(height[,"Area.m2"]),  "sig"=ifelse(height[,"t.NNFD.p.value"] <= 0.05, 1,0))
sigMFDheight = cbind("Island"=height[,"Row.names"],"Size.cat"=as.character(seedmass[,"Size.cat"]), "Metric"=rep(x = "MFD Height", times = nrow(height)), "p.value" = height[,"t.NNFD.p.value"], "Area"=as.character(height[,"Area.m2"]),  "sig"=ifelse(height[,"t.MFD.p.value"] <= 0.05, 1,0))

sigSLA = cbind("Island"=SLA[,"Row.names"], "Size.cat"=as.character(seedmass[,"Size.cat"]),"Metric"=rep(x = "SLA", times = nrow(SLA)), "p.value" = SLA[,"p.value.sla"], "Area"=as.character(SLA[,"Area.m2"]),  "sig"=ifelse(SLA[,"p.value.sla"] <= 0.05, 1,0))
sigNNFDSLA = cbind("Island"=SLA[,"Row.names"],"Size.cat"=as.character(seedmass[,"Size.cat"]), "Metric"=rep(x = "NNFD SLA", times = nrow(SLA)),"p.value" = SLA[,"t.NNFD.p.value"], "Area"=as.character(SLA[,"Area.m2"]),  "sig"=ifelse(SLA[,"t.NNFD.p.value"] <= 0.05, 1,0))
sigMFDSLA = cbind("Island"=SLA[,"Row.names"],"Size.cat"=as.character(seedmass[,"Size.cat"]), "Metric"=rep(x = "MFD SLA", times = nrow(SLA)),"p.value" = SLA[,"t.NNFD.p.value"], "Area"=as.character(SLA[,"Area.m2"]),  "sig"=ifelse(SLA[,"t.MFD.p.value"] <= 0.05, 1,0))

sigleaflet = cbind("Island"=leaflet[,"Row.names"],"Size.cat"=as.character(seedmass[,"Size.cat"]), "Metric"=rep(x = "Leaf Size", times = nrow(leaflet)),  "p.value" = leaflet[,"p.value.leafletSize"], "Area"=as.character(leaflet[,"Area.m2"]), "sig"=ifelse(leaflet[,"p.value.leafletSize"] <= 0.05, 1,0))
sigNNFDleaflet = cbind("Island"=leaflet[,"Row.names"],"Size.cat"=as.character(seedmass[,"Size.cat"]),  "Metric"=rep(x = "NNFD Leaf Size",times = nrow(leaflet)), "p.value" = leaflet[,"t.NNFD.p.value"], "Area"=as.character(leaflet[,"Area.m2"]),  "sig"=ifelse(leaflet[,"t.NNFD.p.value"] <= 0.05, 1,0))
sigMFDleaflet = cbind("Island"=leaflet[,"Row.names"],"Size.cat"=as.character(seedmass[,"Size.cat"]), "Metric"=rep(x = "MFD Leaf Size", times = nrow(leaflet)),"p.value" = leaflet[,"t.NNFD.p.value"], "Area"=as.character(leaflet[,"Area.m2"]),  "sig"=ifelse(leaflet[,"t.MFD.p.value"] <= 0.05, 1,0))

sigleafN = cbind("IleafNnd"=leafN[,"Row.names"], "Size.cat"=as.character(seedmass[,"Size.cat"]),  "Metric"=rep(x = "Leaf N", times = nrow(leafN)),"p.value" = leafN[,"p.value.leafN"], "Area"=as.character(leafN[,"Area.m2"]), "sig"=ifelse(leafN[,"p.value.leafN"] <= 0.05, 1,0))
sigNNFDleafN = cbind("IleafNnd"=leafN[,"Row.names"],"Size.cat"=as.character(seedmass[,"Size.cat"]), "Metric"=rep(x = "NNFD Leaf N", times = nrow(leafN)), "p.value" = leafN[,"t.NNFD.p.value"], "Area"=as.character(leafN[,"Area.m2"]), "sig"=ifelse(leafN[,"t.NNFD.p.value"] <= 0.05, 1,0))
sigMFDleafN = cbind("IleafNnd"=leafN[,"Row.names"],"Size.cat"=as.character(seedmass[,"Size.cat"]), "Metric"=rep(x = "MFD Leaf N", times = nrow(leafN)),"p.value" = leafN[,"t.NNFD.p.value"], "Area"=as.character(leafN[,"Area.m2"]), "sig"=ifelse(leafN[,"t.MFD.p.value"] <= 0.05, 1,0))


obs.function <- as.data.frame(rbind(sigSeed, sigHeight, sigSLA,sigleaflet,sigleafN,
                                    sigNNFDSeed,  sigNNFDheight,  sigNNFDSLA,sigNNFDleaflet,  sigNNFDleafN, 
                                    sigMFDSeed, sigMFDheight,sigMFDSLA,sigMFDleaflet, sigMFDleafN))
                          
neworder <- c('Seed Mass', 'Height','SLA', 'Leaf Size','Leaf N',
              'NNFD Seed Mass', 'NNFD Height',  'NNFD SLA', 'NNFD Leaf Size','NNFD Leaf N',
              'MFD Seed Mass', 'MFD Height','MFD SLA','MFD Leaf Size','MFD Leaf N')

tmpppp <- arrange(transform(obs.function, Metric=factor(Metric,levels=neworder)),Metric)
head(tmpppp)
tmpppp <- tmpppp[tmpppp$Island != 'Unnamed_west_of_Castle_Island',]
tmpppp$p.value <- as.numeric(as.character(tmpppp$p.value))
tmpppp$Island <- gsub(tmpppp$Island, pattern = "___", replacement = " - ")
tmpppp$Island <- gsub(tmpppp$Island, pattern = "__", replacement = "_#")
tmpppp$Island <- gsub(tmpppp$Island, pattern = "_", replacement = " ")
str(tmpppp)

head(phyloObsSum.DNNS.MDNS)
head(tmpppp)
tmpOBS <- as.data.frame(cbind(phyloObsSum.DNNS.MDNS[-7], Measure=rep(x = "Phylogenetic Diversity", times = nrow(phyloObsSum.DNNS.MDNS))))
obsALL <- as.data.frame(rbind(tmpOBS, cbind(tmpppp, Measure=rep(x = "Functional Diversity", times = nrow(tmpppp)))))
obsALL$Size.cat <- factor(obsALL$Size.cat, c("sm", "med", "lg"))

b<- c(0, 0.05, 0.25, 0.5, 0.75, 1)
obs.plot.function <- ggplot(obsALL, aes(y=reorder(factor(Island), as.numeric(as.character(Area))), x=Metric))#, col=c("magenta1", "green3"))
obs.plot.function <- obs.plot.function + geom_tile(aes(fill = p.value), colour="white")
obs.plot.function <- obs.plot.function + scale_fill_gradientn(colours=c("steelblue4", "azure4", "azure3", "azure2", "azure1", "azure"), breaks=b, na.value="white", labels=format(b))
obs.plot.function <- obs.plot.function + geom_point(aes(size=ifelse(as.numeric(as.character(sig)), "dot", "no_dot")))
obs.plot.function <- obs.plot.function + scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")
obs.plot.function <- obs.plot.function + facet_grid(Size.cat~Measure,drop=T,space="free",scales="free", as.table = F)   
obs.plot.function <- obs.plot.function + theme_grey(base_size = 6) + labs(x = "",  y = "") 
obs.plot.function <- obs.plot.function + scale_x_discrete(expand = c(0, 0)) 
obs.plot.function <- obs.plot.function + scale_y_discrete(expand = c(0, 0)) 
obs.plot.function <- obs.plot.function + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  axis.text.x = element_text(size = base_size * 0.8, angle = -45, hjust = 0, colour = "black"))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7,title.position = "top", title.hjust = 0.5))
obs.plot.function <- obs.plot.function + coord_fixed(ratio=1)
obs.plot.function <- obs.plot.function + geom_point(aes(shape=ifelse(is.na(p.value), "is_NA", "not_NA")))
obs.plot.function <- obs.plot.function + scale_shape_manual(values=c(is_NA=4, not_NA=NA), guide="none")
obs.plot.function

#pdf(file="figs/plots/obs.tile.pdf",bg="white", paper="USr")
obs.plot.function
#dev.off()

#################### ses functional 

function.seedmass <- getSesFunctionalDataframe(sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.SeedMass.1000.CSV", 
                          islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="seedMass", metadata=metadata)

function.maxheight <- getSesFunctionalDataframe(sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.Height.null1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="maxHeight", metadata=metadata)

function.sla <- getSesFunctionalDataframe(sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.SLA.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="sla", metadata=metadata)

function.leaflet <- getSesFunctionalDataframe(sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.leafletSize.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="leafletSize", metadata=metadata)

function.leafN <- getSesFunctionalDataframe(sim.output ="output/10_Analyses/FunctionalDiversity/Null/sim.null.distrib.leafN.1000.csv", 
                    islands.sim = SJ_islands.sim, phyloObs = phyloObs, traits=SJtraitLog, traitname="leafN", metadata=metadata)

ses.function <- as.data.frame(rbind(function.seedmass, function.maxheight, function.sla, function.leaflet, function.leafN))
ses.function[ses.function$island =="Secar_Rock",]
names(phylo.exp.DNNS.MDNS)
names(cbind(phylo.exp.DNNS.MDNS[-c(13,14,15,21)], Metric2 = phylo.exp.DNNS.MDNS$Metric))
tmpSES <- as.data.frame(cbind(phylo.exp.DNNS.MDNS[-c(13,14,15,21)], Metric2 = phylo.exp.DNNS.MDNS$Metric, Measure=rep(x = "Phylogenetic Diversity", times = nrow(phylo.exp.DNNS.MDNS))))
sesALL <- as.data.frame(rbind(tmpSES, cbind(ses.function, Measure=rep(x = "Functional Diversity", times = nrow(ses.function)))))

unique(sesALL$Metric2)
sesALL$Metric2 <- revalue(sesALL$Metric2, c('DNNS_inv'= "DNNS", 'MDNS_inv_nat'="MDNS",
                            'NNFD_inv seedMass'="NNFD Seed Mass", 'NNFD_inv maxHeight'="NNFD Height",'NNFD_inv sla'="NNFD SLA", 'NNFD_inv leafletSize'="NNFD Leaf Size", 'NNFD_inv leafN'="NNFD Leaf N",  
                            'MFD_inv_nat seedMass'="MFD Seed Mass", 'MFD_inv_nat maxHeight'="MFD Height",'MFD_inv_nat sla'="MFD SLA", 'MFD_inv_nat leafletSize'="MFD Leaf Size", 'MFD_inv_nat leafN'="MFD Leaf N"))
neworder <- c("DNNS", "MDNS",
              "NNFD Seed Mass", "NNFD Height","NNFD SLA", "NNFD Leaf Size", "NNFD Leaf N",  
              "MFD Seed Mass", "MFD Height","MFD SLA", "MFD Leaf Size", "MFD Leaf N")

sesALL <- arrange(transform(sesALL, Metric2=factor(Metric2,levels=neworder)),Metric2)
head(sesALL)
sesALL$obs.z <- as.numeric(as.character(sesALL$obs.z))
sesALL$island <- gsub(sesALL$island, pattern = "___", replacement = " - ")
sesALL$island <- gsub(sesALL$island, pattern = "__", replacement = "_#")
sesALL$island <- gsub(sesALL$island, pattern = "_", replacement = " ")
str(sesALL)
sesALL$obs.z
sesALL[sesALL$p.value.ranks == 0,]
sesALL.clean <- sesALL
sesALL.clean$obs.z[!is.finite(sesALL.clean$obs.z)] <- NA
is.infinite(sesALL.clean$obs.z)

sesALL.clean[sesALL.clean$island == "Secar Rock",] # complete blank == No trait data for species on that island
sesALL.clean[sesALL.clean$island == "Swirl Rock Central",] # yellow with dot == 
#sesALL.clean <- na.omit(sesALL.clean)
sesALL.clean$Size.cat <- factor(sesALL.clean$Size.cat, c("sm", "med", "lg"))

ses.plot.function <- ggplot(sesALL.clean, aes(y=reorder(factor(island), as.numeric(as.character(Area.m2))), x=Metric2, fill=obs.z))#, col=c("magenta1", "green3"))
ses.plot.function <- ses.plot.function + geom_tile(colour="white", alpha=.75)
ses.plot.function <- ses.plot.function + scale_fill_gradient2(low="cyan2", high="red", mid="beige", na.value="white")
ses.plot.function <- ses.plot.function + geom_point(aes(size=ifelse(sig, "dot", "no_dot")))
ses.plot.function <- ses.plot.function + scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")
ses.plot.function <- ses.plot.function + theme_grey(base_size = 6) + labs(x = "",  y = "") 
ses.plot.function <- ses.plot.function + facet_grid(Size.cat~Measure,space="free",scales="free", as.table = F)   
ses.plot.function <- ses.plot.function + scale_x_discrete(expand = c(0, 0)) 
ses.plot.function <- ses.plot.function + scale_y_discrete(expand = c(0, 0)) 
#ses.plot.function <- ses.plot.function + theme(legend.position = "none", axis.ticks = element_blank(), 
#                             axis.text.x = element_text(size = base_size * 0.8, angle = -45, hjust = 0, colour = "black"))
ses.plot.function <- ses.plot.function + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  axis.text.x = element_text(size = base_size * 0.8, angle = -45, hjust = 0, colour = "black"))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7,title.position = "top", title.hjust = 0.5))
ses.plot.function <- ses.plot.function + coord_fixed(ratio=1)
ses.plot.function <- ses.plot.function + geom_point(aes(shape=ifelse(is.na(obs.z), "is_NA", "not_NA")))
ses.plot.function <- ses.plot.function + scale_shape_manual(values=c(is_NA=4, not_NA=NA), guide="none")
ses.plot.function

#pdf(file="figs/plots/ses.tile.pdf",bg="white", paper="USr")
ses.plot.function
#dev.off()





