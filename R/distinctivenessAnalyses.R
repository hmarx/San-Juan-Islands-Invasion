########### Final Analysis of San Juans Dataset

## Phylogenetic and functional distinctiveness between invasive and native species within each island community 

###### 6 Nov 2014  ##########
###### Hannah E. Marx #######


################################################# Read in Final Datasets #####################################################@
source("analysis.R")
source("R/DistinctivenessFunctions.R")

SJfinalTree
head(SJtraitLog)
head(SJcommNew)
head(metadata)
dim(metadata)

####################################################################################################################################
######################################### PHYLOGENETIC DISTINCTIVENESS #############################################################

######### Calculate observed distances to nearest native (DNNS) and mean distance to native community (MDNS) #######################

### Calculate observed phylo distinctiveness for each species across all islands
com.list <- 1:ncol(SJcommNew) # make a list of the columns in the community data you want the fucntion to itterate across
# Observed find DNNS / MDNS
phyloObs <- lapply(com.list, function(x) phyloDistinct(SJfinalTree, SJcommNew, x)) #apply across all communities
head(phyloObs)
length(phyloObs) #74 islands 
names(phyloObs) <- colnames(SJcommNew)
head(phyloObs[["All_SanJuanIslands"]])
dim(phyloObs[[1]]) #366   7
#write.csv(phyloObs[[1]], file="phyloDistinct.AllSJ.csv")

##### Summarize means for natives and introduced on each island
phyloObsSum.tmp <- lapply(com.list, function(x) summary.DNNS.MDNS(phyloObs[[x]])) #apply summary funciton across all communities...summary.DNNS
names(phyloObsSum.tmp) <- colnames(SJcommNew)
phyloObsSum <- as.data.frame(do.call(cbind, phyloObsSum.tmp), stringsAsFactors=F)
colnames(phyloObsSum) <- colnames(SJcommNew)
phyloObsSum <-  as.data.frame(t(phyloObsSum), stringsAsFactors=F)
#write.csv(phyloObsSum, file="phyloObsSum.MDNS.csv")
head(phyloObsSum)
phyloObsSum <- merge(phyloObsSum, metadata, by=0) # add metadata for islands for plotting and discussing


#pdf("figs/plots/phyloDiv/observed/obs.phylo.SummaryBar.pdf")
sig.obs.phyloDiversty(summaryDF = phyloObsSum, metadataFULL = metadataFULL, plottype = "summary.Bar")
#dev.off()

#pdf("figs/plots/phyloDiv/observed/Observed.DNNS.islIncSize.pdf", width=20, height=10)
sig.obs.phyloDiversty(summaryDF = phyloObsSum, metadataFULL = metadataFULL, plottype = "IslSizeDNNS")
#dev.off()

#pdf("figs/plots/phyloDiv/observed/Observed.MDNS.islIncSize.pdf", width=20, height=10)
sig.obs.phyloDiversty(summaryDF = phyloObsSum, metadataFULL = metadataFULL, plottype = "IslSizeMDNS")
#dev.off()

#pdf("figs/plots/phyloDiv/observed/four.islands.phyloObs.pdf")
sig.obs.phyloDiversty(summaryDF = phyloObsSum, metadataFULL = metadataFULL, plottype = "4islDNNS")
#dev.off()

#pdf("figs/plots/phyloDiv/observed/four.islands.MDNS.islIncSize.pdf")
sig.obs.phyloDiversty(summaryDF = phyloObsSum, metadataFULL = metadataFULL, plottype = "4islMDNS")
#dev.off()

phyloObsDF <- sig.obs.phyloDiversty(summaryDF = phyloObsSum, metadataFULL = metadataFULL, "none")[[2]]


###################### Other summary figures

## Proportion of invasive and native species per island, by increasing island size
#pdf("figs/misc/Proportion.Inv.islIncSize.pdf", width=20, height=10)
ggplot(phyloObsDF,aes(x=reorder(factor(L2),value),fill = Species.Status)) + 
  geom_bar(position = "fill") + 
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values=c("i"= "magenta1", "n"="green3"), labels=c("Native", "Introduced")) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  theme(legend.position="top") +
  scale_x_discrete("Island (increasing size)") +
  theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm")) +
  ggtitle("Proportion of Invasive Species Per Island") + theme(plot.title=element_text(size=rel(1.5)))
#dev.off()


## Richness of invasive and native species per islands, shapes = size categories, 
pdf("figs/misc/SR.percenttips.islsize.pdf", width=20, height=10)
metadataFULL.plot <- ggplot(metadataFULL, aes(x=reorder(factor(Row.names),Area.m2), y=as.numeric(as.character(Total.species)), shape=Size.cat)) +
  geom_point(aes(size=as.numeric(Percent.native.tips)), color="green") +
  geom_point(aes(size=as.numeric(Percent.invasive.tips)), color="magenta") +
  scale_size(range = c(2, 10)) +
  scale_shape_manual(name ="Island Size Categories",values=c("sm"= 1, "med"= 2, "lg"=0), labels=c("sm"="< 1,800", "med"="1,800-34,000", "lg"=">34,000")) +
  theme_bw() +
  theme(legend.position="top") +
  geom_abline(intercept = 0, slope = 0, colour = "grey", size = .5) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  theme(plot.title=element_text(size=rel(1.5))) +
  theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm")) +
  guides(size=guide_legend(title="Percent tips", override.aes = list(colour = "grey")), shape = guide_legend(reverse = TRUE, override.aes = list(colour = "grey"))) +
  ggtitle("") + theme(plot.title=element_text(size=rel(1.5))) +
  scale_y_continuous("species richness") +
  scale_x_discrete("Island (increasing size)") 
metadataFULL.plot
dev.off()


##################################################################################################################################
######################################### FUNCTIONAL DISTINCTIVENESS #############################################################

##################### Calculate the observed difference in trait values between all native and invasive species ##################

######## Number of invasive/natives in spcies list, and phylogeny for each trait
# in the San Juan dataset, all only traits in tree were used, so there is no differnece
head(SJtraitLog)
trait.list <- c(1:ncol(SJtraitLog)) # make a list of the columns in the community data you want the fucntion to itterate across
foo4 <- lapply(trait.list, function(x) traitSummary(SJfinalTree, SJtrait, x)) #summarize all traits
SJtraitSummary <- as.data.frame(do.call(cbind, foo4)) #bring together in a data.frame
SJtraitSummary
#write.csv(SJtraitSummary, file="figs/plots/SJtraitSummary.csv")  ## edited in excel to make Table1_SJtraitSummary.csv


######## Prune trait data to get just species with all trait data for all SJ coommunity (col=1)
SJ_islands <- names(SJcommNew) # names of the islands 
length(SJ_islands) #74
head(SJtraitLog)

################# Summary of each trait, and p-value for difference in mean between native and invasive species

## seed mass
seedmassSummary <- lapply(com.list, function(x) summary.trait.meas(phy=SJfinalTree, community=SJcommNew, traits=SJtraitLog, OneTrait=SJtraitLog[2], col=x))
names(seedmassSummary) <- names(SJcommNew)
seedmassSummary <- as.data.frame(do.call(cbind, seedmassSummary), stringsAsFactors=F)
colnames(seedmassSummary) <- colnames(SJcommNew)
seedmassSummary <-  as.data.frame(t(seedmassSummary), stringsAsFactors=F)
head(seedmassSummary)
#write.csv(seedmassSummary, file="seedmassSummary.csv")


## maximum height
maxheightSummary <- lapply(com.list, function(x) summary.trait.meas(phy=SJfinalTree, community=SJcommNew, traits=SJtraitLog, OneTrait=SJtraitLog[3], col=x))
names(maxheightSummary) <- names(SJcommNew)
maxheightSummary <- as.data.frame(do.call(cbind, maxheightSummary), stringsAsFactors=F)
colnames(maxheightSummary) <- colnames(SJcommNew)
maxheightSummary <-  as.data.frame(t(maxheightSummary), stringsAsFactors=F)
head(maxheightSummary)
#write.csv(maxheightSummary, file="maxheightSummary.csv")

## SLA
SLAtraitSummary <- lapply(com.list, function(x) summary.trait.meas(phy=SJfinalTree, community=SJcommNew, traits=SJtraitLog, OneTrait=SJtraitLog[4], col=x))
names(SLAtraitSummary) <- names(SJcommNew)
SLAtraitSummary <- as.data.frame(do.call(cbind, SLAtraitSummary), stringsAsFactors=F)
colnames(SLAtraitSummary) <- colnames(SJcommNew)
SLAtraitSummary <-  as.data.frame(t(SLAtraitSummary), stringsAsFactors=F)
head(SLAtraitSummary)
#write.csv(SLAtraitSummary, file  ="SLASummary.csv")

## Leaflet size
leafletSummary <- lapply(com.list, function(x) summary.trait.meas(phy=SJfinalTree, community=SJcommNew, traits=SJtraitLog, OneTrait=SJtraitLog[5], col=x))
names(leafletSummary) <- names(SJcommNew)
leafletSummary <- as.data.frame(do.call(cbind, leafletSummary), stringsAsFactors=F)
colnames(leafletSummary) <- colnames(SJcommNew)
leafletSummary <-  as.data.frame(t(leafletSummary), stringsAsFactors=F)
head(leafletSummary)
#write.csv(leafletSummary, file="leafletSummary.csv")

## Leaf N
leafNSummary <- lapply(com.list, function(x) summary.trait.meas(phy=SJfinalTree, community=SJcommNew, traits=SJtraitLog, OneTrait=SJtraitLog[6], col=x))
names(leafNSummary) <- names(SJcommNew)
leafNSummary <- as.data.frame(do.call(cbind, leafNSummary), stringsAsFactors=F)
colnames(leafNSummary) <- colnames(SJcommNew)
leafNSummary <-  as.data.frame(t(leafNSummary), stringsAsFactors=F)
head(leafNSummary)
#write.csv(leafNSummary, file="leafNSummary.csv")

## Plot richness of species with trait data for each island
#### Bring observed trait summaries together all into one dataframe
MyMerge <- function(x, y){
  df <- merge(x, y, by= "row.names", all.x= T, all.y= F)
  rownames(df) <- df$Row.names
  return(df)
}

#colnames(seedmassSummary) <- unlist(lapply(1:length(colnames(seedmassSummary)), function(x) paste(colnames(seedmassSummary)[x], "seedMass", sep = ".")))
#colnames(maxheightSummary) <- unlist(lapply(1:length(colnames(maxheightSummary)), function(x) paste(colnames(maxheightSummary)[x], "maxHeight", sep = ".")))
#colnames(SLAtraitSummary) <- unlist(lapply(1:length(colnames(SLAtraitSummary)), function(x) paste(colnames(SLAtraitSummary)[x], "SLA", sep = ".")))
#colnames(leafletSummary) <- unlist(lapply(1:length(colnames(leafletSummary)), function(x) paste(colnames(leafletSummary)[x], "leafSize", sep = ".")))
#colnames(leafNSummary) <- unlist(lapply(1:length(colnames(leafNSummary)), function(x) paste(colnames(leafNSummary)[x], "leafN", sep = ".")))

masterTraitSummary <- Reduce(MyMerge, list(metadataFULL.rownames, SLAtraitSummary, leafletSummary, seedmassSummary, leafNSummary, maxheightSummary))
masterTraitSummary <- as.data.frame(masterTraitSummary)
head(masterTraitSummary)
names(masterTraitSummary)

#pdf("figs/misc/SR.percenttipsTraits.islsize.pdf", width=20, height=10)
ggplot(masterTraitSummary, aes(x=reorder(factor(Row.names),Area.m2), 
                               y=as.numeric(as.character(Total.species)), shape=Size.cat)) +
  geom_point(aes(x=reorder(factor(Row.names),Area.m2), 
                 y=as.numeric(as.character(Total.tips.with.trait.seedMass)),
                 size=as.numeric(Percent.of.total.tips.with.trait.seedMass)), color="purple") +
  geom_point(data=masterTraitSummary, aes(x=reorder(factor(Row.names),Area.m2), 
                 y=as.numeric(Total.tips.with.trait.maxHeight),
                 size=as.numeric(Percent.of.total.tips.with.trait.maxHeight)), color="red2") + 
  geom_point(aes(x=reorder(factor(Row.names),Area.m2), 
                 y=as.numeric(as.character(Total.tips.with.trait.SLA)),
                 size=as.numeric(Percent.of.total.tips.with.trait.SLA)), color="gold1") +
  geom_point(aes(x=reorder(factor(Row.names),Area.m2), 
                 y=as.numeric(as.character(Total.tips.with.trait.leafSize)),
                 size=as.numeric(Percent.of.total.tips.with.trait.leafSize)), color="green4") +
  geom_point(aes(x=reorder(factor(Row.names),Area.m2), 
                 y=as.numeric(as.character(Total.tips.with.trait.leafN)),
                 size=as.numeric(Percent.of.total.tips.with.trait.leafN)), color="blue") +
  scale_size(range = c(1, 5)) +
  scale_shape_manual(name ="Island Size Categories",values=c("sm"= 1, "med"= 2, "lg"=0), 
                     labels=c("sm"="< 1,800", "med"="1,800-34,000", "lg"=">34,000")) +
  theme_bw() +
  theme(legend.position="top") +
  geom_abline(intercept = 0, slope = 0, colour = "grey", size = .5) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  theme(plot.title=element_text(size=rel(1.5))) +
  theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm")) +
  guides(size=guide_legend(title="Percent tips", override.aes = list(colour = "grey")), 
         shape = guide_legend(reverse = TRUE, override.aes = list(colour = "grey")), 
         colour = "colorbar") +
  ggtitle("") + theme(plot.title=element_text(size=rel(1.5))) +
  scale_y_continuous("species richness") +
  scale_x_discrete("Island (increasing size)") 
#dev.off()


##########  List of species and trait values across each island communtiy, for each trait
com.list <- 1:length(SJ_islands) # make a list of the columns in the community data you want the fucntion to itterate across
seedmass.allSJ <- lapply(SJ_islands, function(x) pruneTrait(phy=SJfinalTree, community=SJcommNew, traits=SJtraitLog, OneTrait=SJtraitLog[2], col=x))
names(seedmass.allSJ) <- names(SJcommNew) 
maxHeight.allSJ <- lapply(SJ_islands, function(x) pruneTrait(phy=SJfinalTree, community=SJcommNew, traits=SJtraitLog, OneTrait=SJtraitLog[3], col=x))
names(maxHeight.allSJ) <- names(SJcommNew) 
sla.allSJ <- lapply(SJ_islands, function(x) pruneTrait(phy=SJfinalTree, community=SJcommNew, traits=SJtraitLog, OneTrait=SJtraitLog[4], col=x))
names(sla.allSJ) <- names(SJcommNew) 
leafletize.allSJ <- lapply(SJ_islands, function(x) pruneTrait(phy=SJfinalTree, community=SJcommNew, traits=SJtraitLog, OneTrait=SJtraitLog[5], col=x))
names(leafletize.allSJ) <- names(SJcommNew) 
leafN.allSJ <- lapply(SJ_islands, function(x) pruneTrait(phy=SJfinalTree, community=SJcommNew, traits=SJtraitLog, OneTrait=SJtraitLog[6], col=x))
names(leafN.allSJ) <- names(SJcommNew) 


######## Plot distribution of observed values for each functional trait on each island, for native and invasive status groups
#SeedMass
#pdf("figs/plots/functionDiv/observed/SeedMass_ObservedValues.pdf", width=20, height=10)
melt.trait.to.meta(list = seedmass.allSJ, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Seed Mass", y.axis.title="Log Seed Mass")
#dev.off()

#maxHeight
pdf("figs/plots/functionDiv/observed/MaxHeight_ObservedValues.pdf", width=20, height=10)
melt.trait.to.meta(list = maxHeight.allSJ, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Maximum Height", y.axis.title="Log Maximum Height")
dev.off()

#SLA
pdf("figs/plots/functionDiv/observed/SLA_ObservedValues.pdf", width=20, height=10)
melt.trait.to.meta(list = sla.allSJ, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="SLA", y.axis.title="SLA")
dev.off()

#leafSize
pdf("figs/plots/functionDiv/observed/LeafletSize_ObservedValues.pdf", width=20, height=10)
melt.trait.to.meta(list = leafletize.allSJ, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Leaflet Size", y.axis.title="Log Leaflet Size")
dev.off()

#leafN
pdf("figs/plots/functionDiv/observed/LeafN_ObservedValues.pdf", width=20, height=10)
melt.trait.to.meta(list = leafN.allSJ, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Leaf Nitrogen", y.axis.title="Log Leaf %N")
dev.off()


######## Just For Four Islands (in each size category) for manuscript figs #######
### Plot distribution of Raw traits 
seedmass.four <- lapply(com.list.4, function(x) pruneTrait(phy=SJfinalTree, community=SJcommNew, traits=SJtraitLog, OneTrait=SJtraitLog[2], col=x))
names(seedmass.four) <- four.islands
maxHeight.four <- lapply(com.list.4, function(x) pruneTrait(phy=SJfinalTree, community=SJcommNew, traits=SJtraitLog, OneTrait=SJtraitLog[3], col=x))
names(maxHeight.four) <- four.islands
sla.four <- lapply(com.list.4, function(x) pruneTrait(phy=SJfinalTree, community=SJcommNew, traits=SJtraitLog, OneTrait=SJtraitLog[4], col=x))
names(sla.four) <- four.islands
leafletize.four <- lapply(com.list.4, function(x) pruneTrait(phy=SJfinalTree, community=SJcommNew, traits=SJtraitLog, OneTrait=SJtraitLog[5], col=x))
names(leafletize.four) <- four.islands
leafN.four <- lapply(com.list.4, function(x) pruneTrait(phy=SJfinalTree, community=SJcommNew, traits=SJtraitLog, OneTrait=SJtraitLog[6], col=x))
names(leafN.four) <-four.islands

#SeedMass
pdf("figs/plots/functionDiv/observed/4Islands_SeedMass_ObservedValues.pdf")
melt.trait.to.meta(list = seedmass.four, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Seed Mass", y.axis.title="Log Seed Mass")
dev.off()
#maxHeight
pdf("figs/plots/functionDiv/observed/4Islands_MaxHeight_ObservedValues.pdf")
melt.trait.to.meta(list = maxHeight.four, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Maximum Height", y.axis.title="Log Maximum Height")
dev.off()
#SLA
pdf("figs/plots/functionDiv/observed/4Islands_SLA_ObservedValues.pdf")
melt.trait.to.meta(list = sla.four, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="SLA", y.axis.title="SLA")
dev.off()
#leafSize
pdf("figs/plots/functionDiv/observed/4Islands_LeafletSize_ObservedValues.pdf")
melt.trait.to.meta(list = leafletize.four, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Leaflet Size", y.axis.title="Log Leaflet Size")
dev.off()
#leafN
pdf("figs/plots/functionDiv/observed/4Islands_LeafN_ObservedValues.pdf")
melt.trait.to.meta(list = leafN.four, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Leaf Nitrogen", y.axis.title="Log Leaf %N")
dev.off()

####################################################################################################################################
######### Calculate observed difference in trait values to nearest native (NNFD) and to mean value odnative community (MFD) #######################

############# Nearest Neighbor Trait Difference (NNFD) ##################################################################

## All Islands (n = 74)
seedmass.distObs <- lapply(com.list, function(x) functionDistinct(output=phyloObs[[x]], traits=SJtraitLog, traitname="seedMass"))
names(seedmass.distObs) <- names(SJcommNew) 

maxheight.distObs <- lapply(com.list, function(x) functionDistinct(output=phyloObs[[x]], traits=SJtraitLog, traitname="maxHeight"))
names(maxheight.distObs) <- names(SJcommNew) 

sla.distObs <- lapply(com.list, function(x) functionDistinct(output=phyloObs[[x]], traits=SJtraitLog, traitname="sla"))
names(sla.distObs) <- names(SJcommNew) 

leaflet.distObs <- lapply(com.list, function(x) functionDistinct(output=phyloObs[[x]], traits=SJtraitLog, traitname="leafletSize"))
names(leaflet.distObs) <- names(SJcommNew) 

leafN.distObs <- lapply(com.list, function(x) functionDistinct(output=phyloObs[[x]], traits=SJtraitLog, traitname="leafN"))
names(leafN.distObs) <- names(SJcommNew) 

#SeedMass
pdf("figs/plots/functionDiv/observed/SeedMass_ObservedNNFD.pdf", width=20, height=10)
melt.NNFD.to.meta(list = seedmass.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Seed Mass NNFD", y.axis.title="NNFD (Seed Mass)")
dev.off()
#maxHeight
pdf("figs/plots/functionDiv/observed/MaxHeight_ObservedNNFD.pdf", width=20, height=10)
melt.NNFD.to.meta(list = maxheight.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Maximum Height NNFD", y.axis.title="NNFD (Maximum Height)")
dev.off()
#SLA
pdf("figs/plots/functionDiv/observed/SLA_ObservedNNFD.pdf", width=20, height=10)
melt.NNFD.to.meta(list = sla.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="SLA NNFD", y.axis.title="NNFD (SLA)")
dev.off()
#leafSize
pdf("figs/plots/functionDiv/observed/Leaflet_ObservedNNFD.pdf", width=20, height=10)
melt.NNFD.to.meta(list = leaflet.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Leaflet Size NNFD", y.axis.title="NNFD (Leaflet Size)")
dev.off()
#leafN
pdf("figs/plots/functionDiv/observed/LeafN_ObservedNNFD.pdf", width=20, height=10)
melt.NNFD.to.meta(list = leafN.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Leaf N NNFD", y.axis.title="NNFD (% Leaf N)")
dev.off()


####### Four Islands
four.islands.seedmass.distObs <- lapply(com.list.4, function(x) functionDistinct(output=phyloObs[[x]], traits=SJtraitLog, traitname="seedMass"))
length(four.islands.seedmass.distObs)
names(four.islands.seedmass.distObs) <- four.islands

four.islands.maxheight.distObs <- lapply(com.list.4, function(x) functionDistinct(output=phyloObs[[x]], traits=SJtraitLog, traitname="maxHeight"))
length(four.islands.maxheight.distObs)
names(four.islands.maxheight.distObs) <-four.islands

four.islands.sla.distObs <- lapply(com.list.4, function(x) functionDistinct(output=phyloObs[[x]], traits=SJtraitLog, traitname="sla"))
length(four.islands.sla.distObs)
names(four.islands.sla.distObs) <- four.islands

four.islands.leaflet.distObs <- lapply(com.list.4, function(x) functionDistinct(output=phyloObs[[x]], traits=SJtraitLog, traitname="leafletSize"))
length(four.islands.leaflet.distObs)
names(four.islands.leaflet.distObs) <-four.islands

four.islands.leafN.distObs <- lapply(com.list.4, function(x) functionDistinct(output=phyloObs[[x]], traits=SJtraitLog, traitname="leafN"))
length(four.islands.leafN.distObs)
names(four.islands.leafN.distObs) <- four.islands


#SLA
pdf("figs/plots/functionDiv/observed/4Islands_SLA_ObservedNNFD.pdf")
melt.NNFD.to.meta(list = four.islands.sla.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="SLA NNFD", y.axis.title="SLA")
dev.off()
#leafSize
pdf("figs/plots/functionDiv/observed/4Islands_Leaflet_ObservedNNFD.pdf")
melt.NNFD.to.meta(list = four.islands.leaflet.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Leaflet Size NNFD", y.axis.title="Leaflet Size")
dev.off()
#SeedMass
pdf("figs/plots/functionDiv/observed/4Islands_SeedMass_ObservedNNFD.pdf")
melt.NNFD.to.meta(list = four.islands.seedmass.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Seed Mass NNFD", y.axis.title="Seed Mass")
dev.off()
#leafN
pdf("figs/plots/functionDiv/observed/4Islands_LeafN_ObservedNNFD.pdf")
melt.NNFD.to.meta(list = four.islands.leafN.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Leaf N NNFD", y.axis.title="% Leaf N")
dev.off()
#maxHeight
pdf("figs/plots/functionDiv/observed/4Islands_MaxHeight_ObservedNNFD.pdf")
melt.NNFD.to.meta(list = four.islands.maxheight.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Maximum Height NNFD", y.axis.title="Maximum Height")
dev.off()


######################################## Mean Functional Distance (MFD) #############################################
####### All Islands
#SeedMass
pdf("figs/plots/functionDiv/observed/SeedMass_ObservedMFD.pdf", width=20, height=10)
melt.MFD.to.meta(list = seedmass.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Seed Mass Mean Functional Distance to Native Community (MFD)", y.axis.title="MFD (Seed Mass)")
dev.off()
#maxHeight
pdf("figs/plots/functionDiv/observed/MaxHeight_ObservedMFD.pdf", width=20, height=10)
melt.MFD.to.meta(list = maxheight.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Maximum Height Mean Functional Distance to Native Community (MFD)", y.axis.title="MFD (Maximum Height)")
dev.off()
#SLA
pdf("figs/plots/functionDiv/observed/SLA_ObservedMFD.pdf", width=20, height=10)
melt.MFD.to.meta(list = sla.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="SLA Mean Functional Distance to Native Community (MFD)", y.axis.title="MFD (SLA)")
dev.off()
#leafSize
pdf("figs/plots/functionDiv/observed/Leaflet_ObservedMFD.pdf", width=20, height=10)
melt.MFD.to.meta(list = leaflet.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Leaflet Size Mean Functional Distance to Native Community (MFD)", y.axis.title="MFD (Leaflet Size)")
dev.off()
#leafN
pdf("figs/plots/functionDiv/observed/LeafN_ObservedMFD.pdf", width=20, height=10)
melt.MFD.to.meta(list = leafN.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Leaf N Mean Functional Distance to Native Community (MFD)", y.axis.title="MFD (% Leaf N)")
dev.off()


####### Four Islands
#SeedMass
pdf("figs/plots/functionDiv/observed/4islands_SeedMass_ObservedMFD.pdf")
melt.MFD.to.meta(list = four.islands.seedmass.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Seed Mass Mean Functional Distance\n to Native Community (MFD)", y.axis.title="Seed Mass")
dev.off()
#maxHeight
pdf("figs/plots/functionDiv/observed/4islands_MaxHeight_ObservedMFD.pdf")
melt.MFD.to.meta(list = four.islands.maxheight.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Maximum Height Mean Functional Distance\n to Native Community (MFD)", y.axis.title="Maximum Height")
dev.off()
#SLA
pdf("figs/plots/functionDiv/observed/4islands_SLA_ObservedMFD.pdf")
melt.MFD.to.meta(list = four.islands.sla.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="SLA Mean Functional Distance\n to Native Community (MFD)", y.axis.title="SLA")
dev.off()
#leafSize
pdf("figs/plots/functionDiv/observed/4islands_Leaflet_ObservedMFD.pdf")
melt.MFD.to.meta(list = four.islands.leaflet.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Leaflet Size Mean Functional Distance\n to Native Community (MFD)", y.axis.title="Leaflet Size")
dev.off()
#leafN
pdf("figs/plots/functionDiv/observed/4islands_LeafN_ObservedMFD.pdf")
melt.MFD.to.meta(list = four.islands.leafN.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Leaf N Mean Functional Distance\n to Native Community (MFD)", y.axis.title="% Leaf N")
dev.off()


################### Summarize means for NNFDi / NNFDn and MFD inv-nat / MFD nat-nat for each trait on each island

### Seed Mass 
SeedMass.distObsSum <- lapply(com.list, function(x) summary.function.diff(seedmass.distObs[[x]])) #apply summary funciton across all communities...summary.DNNS
names(SeedMass.distObsSum) <- colnames(SJcommNew)
SeedMass.distObsSum <- as.data.frame(do.call(cbind, SeedMass.distObsSum), stringsAsFactors=F)
colnames(SeedMass.distObsSum) <- colnames(SJcommNew)
SeedMass.distObsSum <-  as.data.frame(t(SeedMass.distObsSum), stringsAsFactors=F)
#write.csv(SeedMass.distObsSum, file="output/9_Traits/SeedMass.distObsSum.MFD.csv")
head(SeedMass.distObsSum)

sig.obs.functionDiversty(obsdiffSum = SeedMass.distObsSum, obstraitSum = seedmassSummary, metadata = metadata, traitName="Seed Mass")


### Maximum height
Heightsummary.distObsSum <- lapply(com.list, function(x) summary.function.diff(maxheight.distObs[[x]])) #apply summary funciton across all communities...summary.DNNS
names(Heightsummary.distObsSum) <- colnames(SJcommNew)
Heightsummary.distObsSum <- as.data.frame(do.call(cbind, Heightsummary.distObsSum), stringsAsFactors=F)
colnames(Heightsummary.distObsSum) <- colnames(SJcommNew)
Heightsummary.distObsSum <-  as.data.frame(t(Heightsummary.distObsSum), stringsAsFactors=F)
#write.csv(Heightsummary.distObsSum, file="output/9_Traits/Heightsummary.distObsSum.MFD.csv")
head(Heightsummary.distObsSum)

sig.obs.functionDiversty(obsdiffSum = Heightsummary.distObsSum, obstraitSum = maxheightSummary, metadata = metadata, traitName="Max Height")


### SLA
SLAsummary.distObsSum  <- lapply(com.list, function(x) summary.function.diff(sla.distObs[[x]])) #apply summary funciton across all communities...summary.DNNS
names(SLAsummary.distObsSum) <- colnames(SJcommNew)
SLAsummary.distObsSum <- as.data.frame(do.call(cbind, SLAsummary.distObsSum), stringsAsFactors=F)
colnames(SLAsummary.distObsSum) <- colnames(SJcommNew)
SLAsummary.distObsSum <-  as.data.frame(t(SLAsummary.distObsSum), stringsAsFactors=F)
#write.csv(SLAsummary.distObsSum, file="output/9_Traits/SLAsummary.distObsSum.MFD.csv")
head(SLAsummary.distObsSum)

sig.obs.functionDiversty(obsdiffSum = SLAsummary.distObsSum, obstraitSum = SLAtraitSummary, metadata = metadata, traitName="SLA")


### Leaflet 
Leafletsummary.distObsSum <- lapply(com.list, function(x) summary.function.diff(leaflet.distObs[[x]])) #apply summary funciton across all communities...summary.DNNS
names(Leafletsummary.distObsSum) <- colnames(SJcommNew)
Leafletsummary.distObsSum <- as.data.frame(do.call(cbind, Leafletsummary.distObsSum), stringsAsFactors=F)
colnames(Leafletsummary.distObsSum) <- colnames(SJcommNew)
Leafletsummary.distObsSum <-  as.data.frame(t(Leafletsummary.distObsSum), stringsAsFactors=F)
#write.csv(Leafletsummary.distObsSum, file="output/9_Traits/Leafletsummary.distObsSum.MFD.csv")
head(Leafletsummary.distObsSum)

sig.obs.functionDiversty(obsdiffSum = Leafletsummary.distObsSum, obstraitSum = leafletSummary, metadata = metadata, traitName="Leaf Size")


### Leaf N
LeafNsummary.distObsSum <- lapply(com.list, function(x) summary.function.diff(leafN.distObs[[x]])) #apply summary funciton across all communities...summary.DNNS
names(LeafNsummary.distObsSum) <- colnames(SJcommNew)
LeafNsummary.distObsSum <- as.data.frame(do.call(cbind, LeafNsummary.distObsSum), stringsAsFactors=F)
colnames(LeafNsummary.distObsSum) <- colnames(SJcommNew)
LeafNsummary.distObsSum <-  as.data.frame(t(LeafNsummary.distObsSum), stringsAsFactors=F)
#write.csv(LeafNsummary.distObsSum, file="output/9_Traits/LeafNsummary.distObsSum.MFD.csv")
head(LeafNsummary.distObsSum)

sig.obs.functionDiversty(obsdiffSum = LeafNsummary.distObsSum, obstraitSum = leafNSummary, metadata = metadata, traitName="Leaf N")



##################### TOTAL OBSERVED SUMMARY  #########################

aa <- merge(SeedMass.distObsSum.ranm, Heightsummary.distObsSum.ranm, by=0, all.x = T, all.y = T)
bb <- merge(SLAsummary.distObsSum.ranm, Leafletsummary.distObsSum.ranm,by=0, all.x = T, all.y = T)
cc <- merge(aa, bb, by=1, all.x = T, all.y = T)           
dd <- merge(cc, LeafNsummary.distObsSum.ranm, by.x=1, by.y=0, all.x = T, all.y = T)
#write.csv(dd, file="output/9_Traits/masterTraitObs.csv")

##################### TOTAL OBSERVED SUMMARY BAR PLOTS #########################
sigDNNS = cbind("Island"=phyloObsSum[,1], "Size.cat"=as.character(phyloObsSum[,"Size.cat"]), "Metric"=rep(x = "DNNS", times = nrow(phyloObsSum)), "Significance"=ifelse(phyloObsSum[,"t.DNNS.p.value"] <= 0.05, 1,0))
sigMDNS = cbind("Island"=phyloObsSum[,1], "Size.cat"=as.character(phyloObsSum[,"Size.cat"]), "Metric"=rep(x = "MDNS", times = nrow(phyloObsSum)), "Significance"=ifelse(phyloObsSum[,"t.MDNS.p.value"] <= 0.05, 1,0))

sigSeed = cbind("Island"=seedmass[,1], "Size.cat"=as.character(seedmass[,"Size.cat"]), "Metric"=rep(x = "Seed Mass", times = nrow(seedmass)), "Significance"=ifelse(seedmass[,"p.value.seedMass"] <= 0.05, 1,0))
sigNNFDSeed = cbind("Island"=seedmass[,1], "Size.cat"=as.character(seedmass[,"Size.cat"]), "Metric"=rep(x = "NNFD Seed Mass", times = nrow(seedmass)), "Significance"=ifelse(seedmass[,"t.NNFD.p.value"] <= 0.05, 1,0))
sigMFDSeed = cbind("Island"=seedmass[,1], "Size.cat"=as.character(seedmass[,"Size.cat"]), "Metric"=rep(x = "MFD Seed Mass", times = nrow(seedmass)), "Significance"=ifelse(seedmass[,"t.MFD.p.value"] <= 0.05, 1,0))

sigHeight = cbind("Island"=height[,1], "Size.cat"=as.character(height[,"Size.cat"]), "Metric"=rep(x = "Height", times = nrow(height)), "Significance"=ifelse(height[,"p.value.maxHeight"] <= 0.05, 1,0))
sigNNFDheight = cbind("Island"=height[,1], "Size.cat"=as.character(height[,"Size.cat"]), "Metric"=rep(x = "NNFD Height", times = nrow(height)), "Significance"=ifelse(height[,"t.NNFD.p.value"] <= 0.05, 1,0))
sigMFDheight = cbind("Island"=height[,1], "Size.cat"=as.character(height[,"Size.cat"]), "Metric"=rep(x = "MFD Height", times = nrow(height)), "Significance"=ifelse(height[,"t.MFD.p.value"] <= 0.05, 1,0))

sigSLA = cbind("Island"=SLA[,1], "Size.cat"=as.character(SLA[,"Size.cat"]), "Metric"=rep(x = "SLA", times = nrow(SLA)), "Significance"=ifelse(SLA[,"p.value.sla"] <= 0.05, 1,0))
sigNNFDSLA = cbind("Island"=SLA[,1], "Size.cat"=as.character(SLA[,"Size.cat"]), "Metric"=rep(x = "NNFD SLA", times = nrow(SLA)), "Significance"=ifelse(SLA[,"t.NNFD.p.value"] <= 0.05, 1,0))
sigMFDSLA = cbind("Island"=SLA[,1], "Size.cat"=as.character(SLA[,"Size.cat"]), "Metric"=rep(x = "MFD SLA", times = nrow(SLA)), "Significance"=ifelse(SLA[,"t.MFD.p.value"] <= 0.05, 1,0))

sigleaflet = cbind("Island"=leaflet[,1], "Size.cat"=as.character(leaflet[,"Size.cat"]), "Metric"=rep(x = "Leaf Size", times = nrow(leaflet)), "Significance"=ifelse(leaflet[,"p.value.leafletSize"] <= 0.05, 1,0))
sigNNFDleaflet = cbind("Island"=leaflet[,1], "Size.cat"=as.character(leaflet[,"Size.cat"]), "Metric"=rep(x = "NNFD Leaf Size", times = nrow(leaflet)), "Significance"=ifelse(leaflet[,"t.NNFD.p.value"] <= 0.05, 1,0))
sigMFDleaflet = cbind("Island"=leaflet[,1], "Size.cat"=as.character(leaflet[,"Size.cat"]), "Metric"=rep(x = "MFD Leaf Size", times = nrow(leaflet)), "Significance"=ifelse(leaflet[,"t.MFD.p.value"] <= 0.05, 1,0))

sigleafN = cbind("IleafNnd"=leafN[,1], "Size.cat"=as.character(leafN[,"Size.cat"]), "Metric"=rep(x = "Leaf N", times = nrow(leafN)), "Significance"=ifelse(leafN[,"p.value.leafN"] <= 0.05, 1,0))
sigNNFDleafN = cbind("IleafNnd"=leafN[,1], "Size.cat"=as.character(leafN[,"Size.cat"]), "Metric"=rep(x = "NNFD Leaf N", times = nrow(leafN)), "Significance"=ifelse(leafN[,"t.NNFD.p.value"] <= 0.05, 1,0))
sigMFDleafN = cbind("IleafNnd"=leafN[,1], "Size.cat"=as.character(leafN[,"Size.cat"]), "Metric"=rep(x = "MFD Leaf N", times = nrow(leafN)), "Significance"=ifelse(leafN[,"t.MFD.p.value"] <= 0.05, 1,0))

obs <- as.data.frame(rbind(sigSeed, sigNNFDSeed, sigMFDSeed, 
                           sigHeight, sigNNFDheight, sigMFDheight,
                           sigSLA, sigNNFDSLA, sigMFDSLA,
                           sigleaflet, sigNNFDleaflet, sigMFDleaflet,
                          sigleafN, sigNNFDleafN, sigMFDleafN))

neworder <- c('Seed Mass', 'NNFD Seed Mass', 'MFD Seed Mass', 
              'Height', 'NNFD Height', 'MFD Height',
              'SLA', 'NNFD SLA', 'MFD SLA',
              'Leaf Size', 'NNFD Leaf Size', 'MFD Leaf Size',
              'Leaf N', 'NNFD Leaf N', 'MFD Leaf N')

obs2 <- arrange(transform(obs, Metric=factor(Metric,levels=neworder)),Metric)

#pdf("figs/plots/functionDiv/observed/obs.functionDistinct.SummaryBar.pdf")
ggplot(obs2, aes(x=Significance, fill=factor(Size.cat))) + 
  geom_bar(stat="bin") +
  scale_fill_manual(name="Size Category",
                    breaks=c("sm", "med", "lg"),
                    values=c("sm"="dodgerblue4", "med"="orangered3", "lg"="gold1"),
                    labels=c("sm"="small", "med"="medium", "lg"="large")) +
  theme_bw() +
  scale_x_discrete(name="", labels=c("0"="NS", "1"="Significant")) +
  facet_wrap(~Metric, ncol = 3) +
  ggtitle("Observed Functional Distinctiveness \n") 
#dev.off()


