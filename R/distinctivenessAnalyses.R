################### Final Analysis of San Juans Dataset

## Phylogenetic and functional distinctiveness between invasive and native species within each island community 

###### 6 Nov 2014  ##########
###### Hannah E. Marx #######


################################################# Read in Final Datasets #####################################################@
source("analysis.R")

SJfinalTree
head(SJtraitLog)
head(SJcommNew)
head(metadata)
dim(metadata)

####################################################################################################################################
######################################### PHYLOGENETIC DISTINCTIVENESS #############################################################

######### Calculate observed distances to nearest native (MNNPD) and mean distance to native community (MPD) #######################

### Calculate observed phylo distinctiveness for each species across all islands
com.list <- 1:ncol(SJcommNew) # make a list of the columns in the community data you want the fucntion to itterate across
# Observed find MNNPD / MPD
phyloObs <- lapply(com.list, function(x) phyloDistinct(SJfinalTree, SJcommNew, x)) #apply across all communities
head(phyloObs)
length(phyloObs) #74
names(phyloObs) <- colnames(SJcommNew)
names(phyloObs[1]) #"All_SanJuanIslands"
head(phyloObs[[1]])
dim(phyloObs[[1]]) #366   7
#write.csv(phyloObs[[1]], file="phyloDistinct.AllSJ.csv")

##### Summarize means for natives and introduced on each island
phyloObsSum.tmp <- lapply(com.list, function(x) summary.MNNPD.MPD(phyloObs[[x]])) #apply summary funciton across all communities...summary.MNNPD
names(phyloObsSum.tmp) <- colnames(SJcommNew)
phyloObsSum <- as.data.frame(do.call(cbind, phyloObsSum.tmp), stringsAsFactors=F)
colnames(phyloObsSum) <- colnames(SJcommNew)
phyloObsSum <-  as.data.frame(t(phyloObsSum), stringsAsFactors=F)
#write.csv(phyloObsSum, file="phyloObsSum.MPD.csv")
head(phyloObsSum)
phyloObsSum <- merge(phyloObsSum, metadata, by=0) # add metadata for islands for plotting and discussing

# indicate islands with significant difference in observed means between status groups for each island
phyloObsSum[,"Significance"] <- ifelse(phyloObsSum[,"t.MNNPD.p.value"] <= 0.05, 1,
                                            ifelse(phyloObsSum[,"t.MPD.p.value"] <= 0.05, 2, 0)) 

#### summarize significance data by size categories
phyloObsSum$Size.cat <- as.character(phyloObsSum$Size.cat)
phyloObsSum$Size.cat <- factor(phyloObsSum$Size.cat, levels=c("sm", "med", "lg"))


sigDNNS = cbind("Island"=phyloObsSum[,1], "Size.cat"=as.character(phyloObsSum[,"Size.cat"]), "Metric"=rep(x = "MNNPD", times = nrow(phyloObsSum)), "Significance"=ifelse(phyloObsSum[,"t.MNNPD.p.value"] <= 0.05, 1,0))
sigMPD = cbind("Island"=phyloObsSum[,1], "Size.cat"=as.character(phyloObsSum[,"Size.cat"]), "Metric"=rep(x = "MPD", times = nrow(phyloObsSum)), "Significance"=ifelse(phyloObsSum[,"t.MPD.p.value"] <= 0.05, 1,0))
obs.phylo <- as.data.frame(rbind(sigDNNS, sigMPD))

pdf("figs/plots/phyloDiv/observed/obs.phylo.SummaryBar.pdf")
ggplot(obs.phylo, aes(x=Significance, fill=factor(Size.cat))) + 
  geom_bar(stat="bin") +
  scale_fill_manual(name="Size Category",
                    breaks=c("sm", "med", "lg"),
                    values=c("sm"="dodgerblue4", "med"="orangered3", "lg"="gold1"),
                    labels=c("sm"="small", "med"="medium", "lg"="large")) +
  theme_bw() +
  scale_x_discrete(name="", labels=c("0"="NS", "1"="Significant")) +
  facet_wrap(~ Metric) +
  ggtitle("Observed difference of \n Phylogenetic Distances\n") 
dev.off()

##### significance in observed difference bewteen mean MMNPDi / mean MMNPDn 
length(phyloObsSum$t.MNNPD.p.value[which(phyloObsSum$t.MNNPD.p.value >= 0.05)]) # 34 islands NS
length(phyloObsSum$t.MNNPD.p.value[which(phyloObsSum$t.MNNPD.p.value <= 0.05)]) 
# 40 isalnds with significant differnece in mean betweenmean MMNPDi / mean MMNPDn 
phyloObsSum.MNNPD.sig <- phyloObsSum[which(phyloObsSum$t.MNNPD.p.value <= 0.05),]
phyloObsSum[which(phyloObsSum$t.MNNPD.p.value <= 0.005),] # 3
40/74 # 0.5405405 

nrow(merge(metadataFULL[which(metadataFULL$Size.cat == "sm"),], phyloObsSum.MNNPD.sig, by=1)) # 7/74 small islands sig MNNPD 0.09459459
nrow(merge(metadataFULL[which(metadataFULL$Size.cat == "med"),], phyloObsSum.MNNPD.sig, by=1)) # 20/74 medium isl sig MNNPD 0.2702703
nrow(merge(metadataFULL[which(metadataFULL$Size.cat == "lg"),], phyloObsSum.MNNPD.sig, by=1)) # 13/74 large isl sig MNNPD 0.1756757

#### significance in observed difference bewteen mean MPDin / mean MPDnn
length(phyloObsSum$t.MPD.p.value[which(phyloObsSum$t.MPD.p.value >= 0.05)]) # 53 NS
length(phyloObsSum$t.MPD.p.value[which(phyloObsSum$t.MPD.p.value <= 0.05)]) # 21 significant differneces between mean MPD
phyloObsSum.MPD.sig <- phyloObsSum[which(phyloObsSum$t.MPD.p.value <= 0.05),]
phyloObsSum[which(phyloObsSum$t.MPD.p.value <= 0.005),] # 5
21/74 # 0.2837838

nrow(merge(metadataFULL[which(metadataFULL$Size.cat == "sm"),], phyloObsSum.MPD.sig, by=1)) # 2/74 small islands sig MPD 0.02702703
nrow(merge(metadataFULL[which(metadataFULL$Size.cat == "med"),], phyloObsSum.MPD.sig, by=1)) # 4/74 med islands sig MPD 0.05405405
nrow(merge(metadataFULL[which(metadataFULL$Size.cat == "lg"),], phyloObsSum.MPD.sig, by=1)) # 15/74 large isl sig MPD 0.2027027



########################  Plot observed phylogenetic distinctiveness for each island, ordered by increasing island size

### append metadata for island area to column in list of MNNPD/MPD for each island
list <- phyloObs
list.meta <- list()
for (i in 1:length(list)){ 
  tmp <- metadata[names(list[i]), "Area.m2"]
  rep(x=tmp, times=nrow(list[[i]]))
  newlist <- mapply(cbind, list[i], "Area"=tmp, SIMPLIFY=F) 
  list.meta[[i]] <- newlist
}
head(list.meta)

### melt into datafrome for plotting
phyloObs_melt  <- melt(data=list.meta)
head(phyloObs_melt)
dim(phyloObs_melt) # 4835   11

## Regression of differnence in means to island size
model <- lm(as.numeric(as.character(MinDist.Nearest.native)) ~ value , data=phyloObs_melt)
summary(model)
anova(model)
r2.MNNPD <- paste("R^2 = ", signif(summary(model)$r.squared, 3), sep="") #Explained variation / Total variation
p.MNNPD <- paste("p-value = ", signif(anova(model)[[5]][1], 3), "***", sep="")

head(phyloObs_melt)
pdf("figs/plots/phyloDiv/observed/Observed.MNNPD.islIncSize.pdf", width=20, height=10)
MNNPD <- ggplot(phyloObs_melt, aes(x=reorder(factor(L2),value), y=as.numeric(as.character(MinDist.Nearest.native))), position=position_dodge(width=1))#, col=c("magenta1", "green3"))
MNNPD <- MNNPD + geom_boxplot(aes(fill = Species.Status), width = 1)
MNNPD <- MNNPD + geom_smooth(method="lm", se=T, color="black", aes(group=1))
MNNPD <- MNNPD + coord_cartesian(ylim=c(.5, 5000)) 
MNNPD <- MNNPD + scale_x_discrete("Island (increasing size)") #, breaks=seq(0, 80, 10)) 
MNNPD <- MNNPD + scale_y_log10("Log MNNPD") #+ coord_fixed(ratio=4) 
MNNPD <- MNNPD + theme_bw() 
MNNPD <- MNNPD + scale_fill_manual(values=c("i"= "magenta1", "n"="green3"), labels=c("i"= "MMNPDi", "n"="MMNPDn")) ##breaks=rev(factor(SJnew$status)),
MNNPD <- MNNPD + guides(fill=guide_legend(title=""))
MNNPD <- MNNPD + theme(legend.position="top")
MNNPD <- MNNPD + theme(axis.text.x = element_text(angle = -45, hjust = 0))
MNNPD <- MNNPD + ggtitle("Phylogenetic Distance to the Nearest Native Species (MNNPD)") + theme(plot.title=element_text(size=rel(1.5)))
MNNPD <- MNNPD + annotate("text", label=r2.MNNPD, x=3, y=0.9, size=4) #y=max(as.numeric(as.character(SJ_NN_meltNEW$MinDist.Nearest.native))-20)
MNNPD <- MNNPD + annotate("text", label=p.MNNPD, x=4, y=0.7, size=4) 
MNNPD <- MNNPD + theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
MNNPD
dev.off()

four.islands.phyloObs_melt <- phyloObs_melt[phyloObs_melt$L2 %in% four.islands, ]
pdf("figs/plots/phyloDiv/observed/four.islands.phyloObs.pdf")
MNNPD.four <- ggplot(four.islands.phyloObs_melt, aes(x=reorder(factor(L2),value), y=as.numeric(as.character(MinDist.Nearest.native))), position=position_dodge(width=1))#, col=c("magenta1", "green3"))
MNNPD.four <- MNNPD.four+ geom_boxplot(aes(fill = Species.Status), width = 1)
#MNNPD.four+ geom_smooth(method="lm", se=T, color="black", aes(group=1))
MNNPD.four <- MNNPD.four+ coord_cartesian(ylim=c(.5, 5000))
MNNPD.four <- MNNPD.four+ scale_x_discrete("Island (increasing size)") #, breaks=seq(0, 80, 10)) 
MNNPD.four <- MNNPD.four+ scale_y_log10("Log MNNPD") #+ coord_fixed(ratio=4) 
MNNPD.four <- MNNPD.four+ theme_bw() 
MNNPD.four <- MNNPD.four+ scale_fill_manual(values=c("i"= "magenta1", "n"="green3"), labels=c("i"= "MMNPDi", "n"="MMNPDn")) ##breaks=rev(factor(SJnew$status)),
MNNPD.four <- MNNPD.four+ guides(fill=guide_legend(title=""))
MNNPD.four <- MNNPD.four+ theme(legend.position="top")
MNNPD.four <- MNNPD.four+ theme(axis.text.x = element_text(angle = -45, hjust = 0))
MNNPD.four <- MNNPD.four+ ggtitle("Phylogenetic Distance to the Nearest Native Species (MNNPD)") + theme(plot.title=element_text(size=rel(1.5)))
#MNNPD.four+ annotate("text", label=r2.MNNPD, x=5, y=0.09, size=4) #y=max(as.numeric(as.character(SJ_NN_meltNEW$MinDist.Nearest.native))-20)
#MNNPD.four+ annotate("text", label=p.MNNPD, x=5, y=0.05, size=4) 
MNNPD.four <- MNNPD.four+ theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
MNNPD.four
dev.off()

## Regression of differnence in means to island size
modelMPD <- lm(as.numeric(as.character(MeanDist.NativeCommunity)) ~ value , data=phyloObs_melt)
summary(modelMPD)
anova(modelMPD)
r2.MPD <- paste("R^2 = ", signif(summary(modelMPD)$r.squared, 3)) #Explained variation / Total variation
p.MPD <- paste("p-value = ", signif(anova(modelMPD)[[5]][1], 3), "**", sep="")

pdf("Observed.MPD.islIncSize.pdf", width=20, height=10)
MPD <- ggplot(phyloObs_melt, aes(x=reorder(factor(L2),value), y=as.numeric(as.character(MeanDist.NativeCommunity))), position=position_dodge(width=1))#, col=c("magenta1", "green3"))
MPD <- MPD + geom_boxplot(aes(fill = Species.Status), width = 1)
MPD <- MPD + geom_smooth(method="lm", se=T, color="black", aes(group=1))
MPD <- MPD + coord_cartesian(ylim=c(100, 2000))
MPD <- MPD + scale_x_discrete("Island (increasing size)") #, breaks=seq(0, 80, 10)) 
MPD <- MPD + scale_y_log10("Log MPDN")  
MPD <- MPD + theme_bw() 
MPD <- MPD + scale_fill_manual(values=c("i"= "magenta1", "n"="green3"), labels=c("i"= "MPDNi", "n"="MPDNn")) ##breaks=rev(factor(SJnew$status)),
MPD <- MPD + guides(fill=guide_legend(title=""))
MPD <- MPD + theme(legend.position="top")
MPD <- MPD + theme(axis.text.x = element_text(angle = -45, hjust = 0))
MPD <- MPD + ggtitle("Mean Phylogenetic Distance to Native Community (MPDN)") + theme(plot.title=element_text(size=rel(1.5)))
MPD <- MPD + annotate("text", label=r2.MPD, x=3, y=130, size=4) #y=max(as.numeric(as.character(SJ_NN_meltNEW$MinDist.Nearest.native))-20)
MPD <- MPD + annotate("text", label=p.MPD, x=3.5, y=120, size=4) 
MPD <- MPD + theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
MPD
dev.off()

pdf("four.islands.MPD.islIncSize.pdf")
MPD.four <- ggplot(four.islands.NN_melt, aes(x=reorder(factor(L2),value), y=as.numeric(as.character(MeanDist.NativeCommunity))), position=position_dodge(width=1))#, col=c("magenta1", "green3"))
MPD.four <- MPD.four + geom_boxplot(aes(fill = Species.Status), width = 1)
MPD.four <- MPD.four  + coord_cartesian(ylim=c(100, 2000))
#MPD.four <- MPD.four  + geom_smooth(method="lm", se=T, color="black", aes(group=1))
MPD.four <- MPD.four  + scale_x_discrete("Island (increasing size)") #, breaks=seq(0, 80, 10)) 
MPD.four <- MPD.four  + scale_y_log10("Log MPDN")  
MPD.four <- MPD.four  + theme_bw() 
MPD.four <- MPD.four  + scale_fill_manual(values=c("i"= "magenta1", "n"="green3"), labels=c("i"= "MPDi", "n"="MPDn")) ##breaks=rev(factor(SJnew$status)),
MPD.four <- MPD.four  + guides(fill=guide_legend(title=""))
MPD.four <- MPD.four  + theme(legend.position="top")
MPD.four <- MPD.four  + theme(axis.text.x = element_text(angle = -45, hjust = 0))
MPD.four <- MPD.four  + ggtitle("Mean Phylogenetic Distance to Native Community (MPDN)") + theme(plot.title=element_text(size=rel(1.5)))
#MPD.four <- MPD.four  + annotate("text", label=r2.MPD, x=5, y=100, size=4) #y=max(as.numeric(as.character(SJ_NN_meltNEW$MinDist.Nearest.native))-20)
#MPD.four <- MPD.four  + annotate("text", label=p.MPD, x=5, y=90, size=4) 
MPD.four <- MPD.four  + theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
MPD.four 
dev.off()

###################### Other summary figures

## Proportion of invasive and native species per island, by increasing island size
pdf("figs/misc/Proportion.Inv.islIncSize.pdf", width=20, height=10)
ggplot(phyloObs_melt,aes(x=reorder(factor(L2),value),fill = Species.Status)) + 
  geom_bar(position = "fill") + 
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values=c("i"= "magenta1", "n"="green3"), labels=c("Native", "Introduced")) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  theme(legend.position="top") +
  scale_x_discrete("Island (increasing size)") +
  theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm")) +
  ggtitle("Proportion of Invasive Species Per Island") + theme(plot.title=element_text(size=rel(1.5)))
dev.off()


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
traitSummary(phy=SJfinalTree, traits=SJtrait, col=4) #summarize one trait for all SanJuans 
foo4 <- lapply(trait.list, function(x) traitSummary(SJfinalTree, SJtrait, x)) #summarize all traits
SJtraitSummary <- as.data.frame(do.call(cbind, foo4)) #bring together in a data.frame
SJtraitSummary
# write.csv(SJtraitSummary, file="SJtraitSummary.csv")  ## edited in excel to make Table1_SJtraitSummary.csv


######## Prune trait data to get just species with all trait data for all SJ coommunity (col=1)
SJ_islands <- names(SJcommNew) # names of the islands 
length(SJ_islands) #74
head(SJtraitLog)

################# Summary of each trait, and p-value for difference in mean between native and invasive species

## seed mass
seedmassSummary <- lapply(com.list, function(x) commTraitSummary(phy=SJfinalTree, community=SJcommNew, traits=SJtraitLog, OneTrait=SJtraitLog[2], col=x))
names(seedmassSummary) <- names(SJcommNew)
seedmassSummary <- as.data.frame(do.call(cbind, seedmassSummary), stringsAsFactors=F)
colnames(seedmassSummary) <- colnames(SJcommNew)
seedmassSummary <-  as.data.frame(t(seedmassSummary), stringsAsFactors=F)
head(seedmassSummary)
#write.csv(seedmassSummary, file="seedmassSummary.csv")
length(seedmassSummary$p.value[which(seedmassSummary$p.value >= 0.05)]) # 65 NS
length(seedmassSummary$p.value[which(seedmassSummary$p.value <= 0.05)]) # 9 significant differneces between means 
sigseedmassSummary <- (seedmassSummary[which(seedmassSummary$p.value <= 0.05),]) # 9 significant differneces between means 
9/74 #0.1216216

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigseedmassSummary, by=0)) # 0 small islands sig 
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigseedmassSummary, by=0)) # 6/74 med islands sig  0.08108108
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigseedmassSummary, by=0)) # 3/74 large isl sig  0.04054054


## maximum height
maxheightSummary <- lapply(com.list, function(x) commTraitSummary(phy=SJfinalTree, community=SJcommNew, traits=SJtraitLog, OneTrait=SJtraitLog[3], col=x))
names(maxheightSummary) <- names(SJcommNew)
maxheightSummary <- as.data.frame(do.call(cbind, maxheightSummary), stringsAsFactors=F)
colnames(maxheightSummary) <- colnames(SJcommNew)
maxheightSummary <-  as.data.frame(t(maxheightSummary), stringsAsFactors=F)
head(maxheightSummary)
#write.csv(maxheightSummary, file="maxheightSummary.csv")
length(maxheightSummary$p.value[which(maxheightSummary$p.value >= 0.05)]) # 29 NS
length(maxheightSummary$p.value[which(maxheightSummary$p.value <= 0.05)]) # 45 significant differneces between means 
sigmaxheightSummary <- (maxheightSummary[which(maxheightSummary$p.value <= 0.05),]) # 45 significant differneces between means 
45/74 # 0.6081081

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigmaxheightSummary, by=0)) # 5/74 small islands sig 0.06756757
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigmaxheightSummary, by=0)) # 23/74 med islands sig  0.3108108
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigmaxheightSummary, by=0)) # 17/74 large isl sig  0.2297297


## SLA
SLAtraitSummary <- lapply(com.list, function(x) commTraitSummary(phy=SJfinalTree, community=SJcommNew, traits=SJtraitLog, OneTrait=SJtraitLog[4], col=x))
names(SLAtraitSummary) <- names(SJcommNew)
SLAtraitSummary <- as.data.frame(do.call(cbind, SLAtraitSummary), stringsAsFactors=F)
colnames(SLAtraitSummary) <- colnames(SJcommNew)
SLAtraitSummary <-  as.data.frame(t(SLAtraitSummary), stringsAsFactors=F)
head(SLAtraitSummary)
#write.csv(SLAtraitSummary, file="SLASummary.csv")
length(SLAtraitSummary$p.value[which(SLAtraitSummary$p.value >= 0.05)]) # 31 NS
length(SLAtraitSummary$p.value[which(SLAtraitSummary$p.value <= 0.05)]) # 43 significant differneces between means 
sigSLAtraitSummary <- SLAtraitSummary[which(SLAtraitSummary$p.value <= 0.05),] 
dim(SLAtraitSummary)
head(SLAtraitSummary)
43/74 # 0.5810811

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigSLAtraitSummary, by=0)) # 7/74 small islands sig 0.09459459
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigSLAtraitSummary, by=0)) # 21/74 med islands sig  0.2837838
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigSLAtraitSummary, by=0)) # 15/74 large isl sig  0.2027027

## Leaflet size
leafletSummary <- lapply(com.list, function(x) commTraitSummary(phy=SJfinalTree, community=SJcommNew, traits=SJtraitLog, OneTrait=SJtraitLog[5], col=x))
names(leafletSummary) <- names(SJcommNew)
leafletSummary <- as.data.frame(do.call(cbind, leafletSummary), stringsAsFactors=F)
colnames(leafletSummary) <- colnames(SJcommNew)
leafletSummary <-  as.data.frame(t(leafletSummary), stringsAsFactors=F)
head(leafletSummary)
#write.csv(leafletSummary, file="leafletSummary.csv")
length(leafletSummary$p.value[which(leafletSummary$p.value >= 0.05)]) # 69 NS
length(leafletSummary$p.value[which(leafletSummary$p.value <= 0.05)]) # 5 significant differneces between means 
sigleafletSummary <- (leafletSummary[which(leafletSummary$p.value <= 0.05),]) 
5/74 # 0.06756757

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigleafletSummary, by=0)) # 1/74 small islands sig 0.01351351
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigleafletSummary, by=0)) # 2/74 med islands sig  0.02702703
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigleafletSummary, by=0)) # 2/74 large isl sig  0.02702703

## Leaf N
leafNSummary <- lapply(com.list, function(x) commTraitSummary(phy=SJfinalTree, community=SJcommNew, traits=SJtraitLog, OneTrait=SJtraitLog[6], col=x))
names(leafNSummary) <- names(SJcommNew)
leafNSummary <- as.data.frame(do.call(cbind, leafNSummary), stringsAsFactors=F)
colnames(leafNSummary) <- colnames(SJcommNew)
leafNSummary <-  as.data.frame(t(leafNSummary), stringsAsFactors=F)
head(leafNSummary)
#write.csv(leafNSummary, file="leafNSummary.csv")
length(leafNSummary$p.value[which(leafNSummary$p.value >= 0.05)]) # 43 NS
length(leafNSummary$p.value[which(leafNSummary$p.value <= 0.05)]) # 31 significant differneces between means 
sigleafNSummary<- (leafNSummary[which(leafNSummary$p.value <= 0.05),]) # 31 significant differneces between means 
31/74 # 0.4189189

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigleafNSummary, by=0)) # 2/74 small islands sig 0.02702703
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigleafNSummary, by=0)) # 16/74 med islands sig  0.2162162
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigleafNSummary, by=0)) # 13/74 large isl sig  0.1756757


#### Bring observed trait summaries together all into one dataframe
MyMerge <- function(x, y){
  df <- merge(x, y, by= "row.names", all.x= T, all.y= F)
  rownames(df) <- df$Row.names
  #df$Row.names  <- NULL
  return(df)
}
masterTraitSummary <- Reduce(MyMerge, list(metadataFULL.rownames, SLAtraitSummary, leafletSummary, seedmassSummary, leafNSummary, maxheightSummary))
masterTraitSummary <- as.data.frame(masterTraitSummary)
head(masterTraitSummary)
str(masterTraitSummary)

## Plot richness of species with trait data for each island
pdf("figs/misc/SR.percenttipsTraits.islsizeNEW.pdf", width=20, height=10)
ggplot(masterTraitSummary, aes(x=reorder(factor(Row.names),Area.m2), 
                               y=as.numeric(as.character(Total.species)), shape=Size.cat)) +
  geom_point(aes(x=reorder(factor(Row.names),Area.m2), 
                 y=as.numeric(as.character(Total.tips.with.seedMass)),
                 size=as.numeric(Percent.of.total.tips.with.seedMass)), color="purple") +
  geom_point(data=masterTraitSummary, aes(x=reorder(factor(Row.names),Area.m2), 
                 y=as.numeric(Total.tips.with.maxHeight),
                 size=as.numeric(Percent.of.total.tips.with.maxHeight)), color="red2") + 
  geom_point(aes(x=reorder(factor(Row.names),Area.m2), 
                 y=as.numeric(as.character(Total.tips.with.sla)),
                 size=as.numeric(Percent.of.total.tips.with.sla)), color="gold1") +
  geom_point(aes(x=reorder(factor(Row.names),Area.m2), 
                 y=as.numeric(as.character(Total.tips.with.leafletSize)),
                 size=as.numeric(Percent.of.total.tips.with.leafletSize)), color="green4") +
  geom_point(aes(x=reorder(factor(Row.names),Area.m2), 
                 y=as.numeric(as.character(Total.tips.with.leafN)),
                 size=as.numeric(Percent.of.total.tips.with.leafN)), color="blue") +
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
dev.off()


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
pdf("figs/plots/functionDiv/observed/SeedMass_ObservedValues.pdf", width=20, height=10)
melt.trait.to.meta(list = seedmass.allSJ, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Seed Mass", y.axis.title="Log Seed Mass")
dev.off()

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
######### Calculate observed difference in trait values to nearest native (MNNFD) and to mean value odnative community (MFD) #######################

############# Nearest Neighbor Trait Difference (MNNFD) ##################################################################
## Entire archipelago
seedMass_SJ <- functionDstinct(output=phyloObs[[1]], traits=SJtraitLog, traitname="seedMass")
maxHeight_SJ <- functionDistinct(output=phyloObs[[1]], traits=SJtraitLog, traitname="maxHeight")
SLA_SJ <- functionDistinct(output=phyloObs[[1]], traits=SJtraitLog, traitname="sla") #one community, one trait
leafletSize_SJ <- functionDistinct(output=phyloObs[[1]], traits=SJtraitLog, traitname="leafletSize")
leafN_SJ <- functionDistinct(output=phyloObs[[1]], traits=SJtraitLog, traitname="leafN")
a2 <- plot.functionDistinct.Obs(MNNFDoutput=seedMass_SJ, islandname="", traitname="seedMass", metric = 2)
b2 <- plot.functionDistinct.Obs(MNNFDoutput=maxHeight_SJ, islandname="", traitname="maxHeight", metric = 2)
c2 <- plot.functionDistinct.Obs(MNNFDoutput=SLA_SJ, islandname="", traitname="sla", metric = 2)
d2 <- plot.functionDistinct.Obs(MNNFDoutput=leafletSize_SJ, islandname="", traitname="leafletSize", metric = 2)
e2 <- plot.functionDistinct.Obs(MNNFDoutput=leafN_SJ, islandname="", traitname="leafN", metric = 2)
grid.arrange(a2, b2, c2, d2, e2, ncol=6, main=textGrob(vjust = .8, hjust=.5,"MNNFD \nAll San Juan Islands",gp=gpar(cex=1.5,font=2)))

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
pdf("figs/plots/functionDiv/observed/SeedMass_ObservedMNNFD.pdf", width=20, height=10)
melt.MNNFD.to.meta(list = seedmass.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Seed Mass MNNFD", y.axis.title="Seed Mass")
dev.off()
#maxHeight
pdf("figs/plots/functionDiv/observed/MaxHeight_ObservedMNNFD.pdf", width=20, height=10)
melt.MNNFD.to.meta(list = maxheight.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Maximum Height MNNFD", y.axis.title="Maximum Height")
dev.off()
#SLA
pdf("figs/plots/functionDiv/observed/SLA_ObservedMNNFD.pdf", width=20, height=10)
melt.MNNFD.to.meta(list = sla.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="SLA MNNFD", y.axis.title="SLA")
dev.off()
#leafSize
pdf("figs/plots/functionDiv/observed/Leaflet_ObservedMNNFD.pdf", width=20, height=10)
melt.MNNFD.to.meta(list = leaflet.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Leaflet Size MNNFD", y.axis.title="Leaflet Size")
dev.off()
#leafN
pdf("figs/plots/functionDiv/observed/LeafN_ObservedMNNFD.pdf", width=20, height=10)
melt.MNNFD.to.meta(list = leafN.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Leaf N MNNFD", y.axis.title="% Leaf N")
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
pdf("figs/plots/functionDiv/observed/4Islands_SLA_ObservedMNNFD.pdf")
melt.MNNFD.to.meta(list = four.islands.sla.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="SLA MNNFD", y.axis.title="SLA")
dev.off()
#leafSize
pdf("figs/plots/functionDiv/observed/4Islands_Leaflet_ObservedMNNFD.pdf")
melt.MNNFD.to.meta(list = four.islands.leaflet.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Leaflet Size MNNFD", y.axis.title="Leaflet Size")
dev.off()
#SeedMass
pdf("figs/plots/functionDiv/observed/4Islands_SeedMass_ObservedMNNFD.pdf")
melt.MNNFD.to.meta(list = four.islands.seedmass.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Seed Mass MNNFD", y.axis.title="Seed Mass")
dev.off()
#leafN
pdf("figs/plots/functionDiv/observed/4Islands_LeafN_ObservedMNNFD.pdf")
melt.MNNFD.to.meta(list = four.islands.leafN.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Leaf N MNNFD", y.axis.title="% Leaf N")
dev.off()
#maxHeight
pdf("figs/plots/functionDiv/observed/4Islands_MaxHeight_ObservedMNNFD.pdf")
melt.MNNFD.to.meta(list = four.islands.maxheight.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Maximum Height MNNFD", y.axis.title="Maximum Height")
dev.off()



######################################## Mean Functional Distance (MFD) #############################################
####### All Islands
#SeedMass
pdf("figs/plots/functionDiv/observed/SeedMass_ObservedMFD.pdf", width=20, height=10)
melt.MFD.to.meta(list = seedmass.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Seed Mass Mean Functional Distance to Native Community (MFP)", y.axis.title="Seed Mass")
dev.off()
#maxHeight
pdf("figs/plots/functionDiv/observed/MaxHeight_ObservedMFD.pdf", width=20, height=10)
melt.MFD.to.meta(list = maxheight.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Maximum Height Mean Functional Distance to Native Community (MFP)", y.axis.title="Maximum Height")
dev.off()
#SLA
pdf("figs/plots/functionDiv/observed/SLA_ObservedMFD.pdf", width=20, height=10)
melt.MFD.to.meta(list = sla.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="SLA Mean Functional Distance to Native Community (MFP)", y.axis.title="SLA")
dev.off()
#leafSize
pdf("figs/plots/functionDiv/observed/Leaflet_ObservedMFD.pdf", width=20, height=10)
melt.MFD.to.meta(list = leaflet.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Leaflet Size Mean Functional Distance to Native Community (MFP)", y.axis.title="Leaflet Size")
dev.off()
#leafN
pdf("figs/plots/functionDiv/observed/LeafN_ObservedMFD.pdf", width=20, height=10)
melt.MFD.to.meta(list = leafN.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Leaf N Mean Functional Distance to Native Community (MFP)", y.axis.title="% Leaf N")
dev.off()


####### Four Islands
#SeedMass
pdf("figs/plots/functionDiv/observed/4islands_SeedMass_ObservedMFD.pdf")
melt.MFD.to.meta(list = four.islands.seedmass.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Seed Mass Mean Functional Distance\n to Native Community (MFP)", y.axis.title="Seed Mass")
dev.off()
#maxHeight
pdf("figs/plots/functionDiv/observed/4islands_MaxHeight_ObservedMFD.pdf")
melt.MFD.to.meta(list = four.islands.maxheight.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Maximum Height Mean Functional Distance\n to Native Community (MFP)", y.axis.title="Maximum Height")
dev.off()
#SLA
pdf("figs/plots/functionDiv/observed/4islands_SLA_ObservedMFD.pdf")
melt.MFD.to.meta(list = four.islands.sla.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="SLA Mean Functional Distance\n to Native Community (MFP)", y.axis.title="SLA")
dev.off()
#leafSize
pdf("figs/plots/functionDiv/observed/4islands_Leaflet_ObservedMFD.pdf")
melt.MFD.to.meta(list = four.islands.leaflet.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Leaflet Size Mean Functional Distance\n to Native Community (MFP)", y.axis.title="Leaflet Size")
dev.off()
#leafN
pdf("figs/plots/functionDiv/observed/4islands_LeafN_ObservedMFD.pdf")
melt.MFD.to.meta(list = four.islands.leafN.distObs, metadata = metadata, meta.data.column.name= "Area.m2", plot.title="Leaf N Mean Functional Distance\n to Native Community (MFP)", y.axis.title="% Leaf N")
dev.off()


################### Summarize means for MNNFDi / MNNFDn and MFD inv-nat / MFD nat-nat for each trait on each island

### Seed Mass 
SeedMass.distObsSum <- lapply(com.list, function(x) functionObsSum(seedmass.distObs[[x]])) #apply summary funciton across all communities...summary.MNNPD
names(SeedMass.distObsSum) <- colnames(SJcommNew)
SeedMass.distObsSum <- as.data.frame(do.call(cbind, SeedMass.distObsSum), stringsAsFactors=F)
colnames(SeedMass.distObsSum) <- colnames(SJcommNew)
SeedMass.distObsSum <-  as.data.frame(t(SeedMass.distObsSum), stringsAsFactors=F)
#write.csv(SeedMass.distObsSum, file="SeedMass.distObsSum.MFD.csv")
head(SeedMass.distObsSum)
SeedMass.distObsSum.ranm <- SeedMass.distObsSum[SeedMass.distObsSum$t.MNNFD.p.value != "NA",]

seedtmp <- transform(merge(seedmassSummary, SeedMass.distObsSum, by=0), row.names=Row.names,Row.names=NULL)
seedmass <- (merge(seedtmp, metadata, by=0))
dim(seedtmp)
seedmass[,"Significance"] <- ifelse(seedmass[,"p.value.seedMass"] <= 0.05, 1,
                                    ifelse(seedmass[,"t.MNNFD.p.value"] <= 0.05, 2,
                                           ifelse(seedmass[,"t.MFD.p.value"] <= 0.05, 3, 0)))
head(seedmass)
#### summarize significance data by size categories
seedmass$Size.cat <- as.character(seedmass$Size.cat)
seedmass$Size.cat <- factor(seedmass$Size.cat, levels=c("sm", "med", "lg"))


##### mean MMNPDi / mean MMNPDn 
length(SeedMass.distObsSum.ranm$t.MNNFD.p.value[which(SeedMass.distObsSum.ranm$t.MNNFD.p.value >= 0.05)]) # 72 NS
length(SeedMass.distObsSum.ranm$t.MNNFD.p.value[which(SeedMass.distObsSum.ranm$t.MNNFD.p.value <= 0.05)]) # 1 significant differneces between mean MNNPD 
sigSeedMass.distObsSum <- (SeedMass.distObsSum.ranm[which(SeedMass.distObsSum.ranm$t.MNNFD.p.value <= 0.05), ]) # 1 significant differneces between mean MNNPD 
1/73 #0.01369863
#### mean MPDin / mean MPDnn
length(SeedMass.distObsSum.ranm$t.MFD.p.value[which(SeedMass.distObsSum.ranm$t.MFD.p.value >= 0.05)]) # 62 NS
length(SeedMass.distObsSum.ranm$t.MFD.p.value[which(SeedMass.distObsSum.ranm$t.MFD.p.value <= 0.05)]) # 11 significant differneces between mean MPD
sigSeedMasssummaryMPD <- (SeedMass.distObsSum.ranm[which(SeedMass.distObsSum.ranm$t.MFD.p.value <= 0.05),]) # 11 significant differneces between mean MPD
11/73 #0.1506849

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigSeedMass.distObsSum, by=0)) # 0 small islands sig MPD
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigSeedMass.distObsSum, by=0)) # 1/73 med islands sig MPD 0.01369863
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigSeedMass.distObsSum, by=0)) # 0 large isl sig MPD

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigSeedMasssummaryMPD, by=0)) #2/73 small islands sig MPD 0.02739726
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigSeedMasssummaryMPD, by=0)) #4/73 med islands sig MPD 0.05479452
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigSeedMasssummaryMPD, by=0)) # 5/73 large isl sig MPD 0.06849315


### Maximum height
Heightsummary.distObsSum <- lapply(com.list, function(x) functionObsSum(maxheight.distObs[[x]])) #apply summary funciton across all communities...summary.MNNPD
names(Heightsummary.distObsSum) <- colnames(SJcommNew)
Heightsummary.distObsSum <- as.data.frame(do.call(cbind, Heightsummary.distObsSum), stringsAsFactors=F)
colnames(Heightsummary.distObsSum) <- colnames(SJcommNew)
Heightsummary.distObsSum <-  as.data.frame(t(Heightsummary.distObsSum), stringsAsFactors=F)
#write.csv(Heightsummary.distObsSum, file="Heightsummary.distObsSum.MFD.csv")
head(Heightsummary.distObsSum)
Heightsummary.distObsSum.ranm <- Heightsummary.distObsSum[Heightsummary.distObsSum$t.MNNFD.p.value != "NA",]

heighttmp <- transform(merge(maxheightSummary, Heightsummary.distObsSum, by=0), row.names=Row.names,Row.names=NULL)
height <- (merge(heighttmp, metadata, by=0))
head(height)
height[,"Significance"] <- ifelse(height[,"p.value.maxHeight"] <= 0.05, 1,
                                  ifelse(height[,"t.MNNFD.p.value"] <= 0.05, 2,
                                         ifelse(height[,"t.MFD.p.value"] <= 0.05, 3, 0)))

#### summarize significance data by size categories
height$Size.cat <- as.character(height$Size.cat)
height$Size.cat <- factor(height$Size.cat, levels=c("sm", "med", "lg"))


##### mean MMNPDi / mean MMNPDn 
length(Heightsummary.distObsSum.ranm$t.MNNFD.p.value[which(Heightsummary.distObsSum.ranm$t.MNNFD.p.value >= 0.05)]) # 51 NS
length(Heightsummary.distObsSum.ranm$t.MNNFD.p.value[which(Heightsummary.distObsSum.ranm$t.MNNFD.p.value <= 0.05)]) # 20 significant differneces between mean MNNPD 
sigHeightsummary.distObsSum <- (Heightsummary.distObsSum.ranm[which(Heightsummary.distObsSum.ranm$t.MNNFD.p.value <= 0.05),]) # 20 significant differneces between mean MNNPD 
20/71 #0.2816901
#### mean MPDin / mean MPDnn
length(Heightsummary.distObsSum.ranm$t.MFD.p.value[which(Heightsummary.distObsSum.ranm$t.MFD.p.value >= 0.05)]) # 54 NS
length(Heightsummary.distObsSum.ranm$t.MFD.p.value[which(Heightsummary.distObsSum.ranm$t.MFD.p.value <= 0.05)]) # 17 significant differneces between mean MPD
sigHeightsummaryMPD <- (Heightsummary.distObsSum.ranm[which(Heightsummary.distObsSum.ranm$t.MFD.p.value <= 0.05),]) # 17 significant differneces between mean MPD
17/71 #0.2394366

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigHeightsummary.distObsSum, by=0)) # 4/71 small islands sig MPD 0.05633803
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigHeightsummary.distObsSum, by=0)) # 12/71 med islands sig MPD 0.1690141
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigHeightsummary.distObsSum, by=0)) # 4/71 large isl sig MPD  0.05633803

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigHeightsummaryMPD, by=0)) # 3/71 small islands sig MPD 0.04225352
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigHeightsummaryMPD, by=0)) # 8/71 med islands sig MPD 0.1126761
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigHeightsummaryMPD, by=0)) # 6/71 large isl sig MPD 0.08450704

### SLA
SLAsummary.distObsSum  <- lapply(com.list, function(x) functionObsSum(sla.distObs[[x]])) #apply summary funciton across all communities...summary.MNNPD
names(SLAsummary.distObsSum) <- colnames(SJcommNew)
SLAsummary.distObsSum <- as.data.frame(do.call(cbind, SLAsummary.distObsSum), stringsAsFactors=F)
colnames(SLAsummary.distObsSum) <- colnames(SJcommNew)
SLAsummary.distObsSum <-  as.data.frame(t(SLAsummary.distObsSum), stringsAsFactors=F)
#write.csv(SLAsummary.distObsSum, file="SLAsummary.distObsSum.MFD.csv")
head(SLAsummary.distObsSum)
SLAsummary.distObsSum.ranm <- SLAsummary.distObsSum[SLAsummary.distObsSum$meanMNNFDnatives != "NA",]

SLAtmp <- transform(merge(SLAtraitSummary, SLAsummary.distObsSum, by=0), row.names=Row.names, Row.names=NULL)
SLA <- (merge(SLAtmp, metadata, by=0))

SLA[,"Significance"] <- ifelse(SLA[,"p.value.sla"] <= 0.05, 1,
                              ifelse(SLA[,"t.MNNFD.p.value"] <= 0.05, 2,
                               ifelse(SLA[,"t.MFD.p.value"] <= 0.05, 3, 0)))

#### summarize significance data by size categories
SLA$Size.cat <- as.character(SLA$Size.cat)
SLA$Size.cat <- factor(SLA$Size.cat, levels=c("sm", "med", "lg"))

##### mean MMNPDi / mean MMNPDn 
length(SLAsummary.distObsSum.ranm$t.MNNFD.p.value[which(SLAsummary.distObsSum.ranm$t.MNNFD.p.value >= 0.05)]) # 51 NS
length(SLAsummary.distObsSum.ranm$t.MNNFD.p.value[which(SLAsummary.distObsSum.ranm$t.MNNFD.p.value <= 0.05)]) # 10 
sigSLAsummary.distObsSum <- (SLAsummary.distObsSum.ranm[which(SLAsummary.distObsSum.ranm$t.MNNFD.p.value <= 0.05),]) # 10 
10/61 #0.1639344
#### mean MPDin / mean MPDnn
length(SLAsummary.distObsSum.ranm$t.MFD.p.value[which(SLAsummary.distObsSum.ranm$t.MFD.p.value >= 0.05)]) # 57 NS
length(SLAsummary.distObsSum.ranm$t.MFD.p.value[which(SLAsummary.distObsSum.ranm$t.MFD.p.value <= 0.05)]) # 4 significant differneces between mean MPD
sigSLAsummaryMPD <- SLAsummary.distObsSum.ranm[which(SLAsummary.distObsSum.ranm$t.MFD.p.value <= 0.05),] # 4 significant differneces between mean MPD
4/61 #0.06557377

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigSLAsummary.distObsSum, by=0)) # 1/61 small islands sig MPD 0.01639344
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigSLAsummary.distObsSum, by=0)) # 6/61 med islands sig MPD 0.09836066
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigSLAsummary.distObsSum, by=0)) # 3/61 large isl sig MPD 0.04918033

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigSLAsummaryMPD, by=0)) # 2/61 small islands sig MPD 0.03278689
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigSLAsummaryMPD, by=0)) # 1/61 med islands sig MPD 0.01639344
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigSLAsummaryMPD, by=0)) # 1/61 large isl sig MPD 0.01639344


### Leaflet 
Leafletsummary.distObsSum <- lapply(com.list, function(x) functionObsSum(leaflet.distObs[[x]])) #apply summary funciton across all communities...summary.MNNPD
names(Leafletsummary.distObsSum) <- colnames(SJcommNew)
Leafletsummary.distObsSum <- as.data.frame(do.call(cbind, Leafletsummary.distObsSum), stringsAsFactors=F)
colnames(Leafletsummary.distObsSum) <- colnames(SJcommNew)
Leafletsummary.distObsSum <-  as.data.frame(t(Leafletsummary.distObsSum), stringsAsFactors=F)
#write.csv(Leafletsummary.distObsSum, file="Leafletsummary.distObsSum.MFD.csv")
head(Leafletsummary.distObsSum)
Leafletsummary.distObsSum.ranm <- Leafletsummary.distObsSum[Leafletsummary.distObsSum$t.MNNFD.p.value != "NA",]

leaflettmp <-transform(merge(leafletSummary, Leafletsummary.distObsSum, by=0), row.names=Row.names,Row.names=NULL)
leaflet <- (merge(leaflettmp, metadata,by=0))
head(leaflet)
leaflet[,"Significance"] <- ifelse(leaflet[,"p.value.leafletSize"] <= 0.05, 1,
                               ifelse(leaflet[,"t.MNNFD.p.value"] <= 0.05, 2,
                                      ifelse(leaflet[,"t.MFD.p.value"] <= 0.05, 3, 0)))

#### summarize significance data by size categories
leaflet$Size.cat <- as.character(leaflet$Size.cat)
leaflet$Size.cat <- factor(leaflet$Size.cat, levels=c("sm", "med", "lg"))

##### mean MMNPDi / mean MMNPDn 
length(Leafletsummary.distObsSum.ranm$t.MNNFD.p.value[which(Leafletsummary.distObsSum.ranm$t.MNNFD.p.value >= 0.05)]) # 54 NS
length(Leafletsummary.distObsSum.ranm$t.MNNFD.p.value[which(Leafletsummary.distObsSum.ranm$t.MNNFD.p.value <= 0.05)]) # 3 significant differneces between mean MNNPD 
Leafletsummary.distObsSumsig <- (Leafletsummary.distObsSum[which(Leafletsummary.distObsSum.ranm$t.MNNFD.p.value <= 0.05),]) # 3 significant differneces between mean MNNPD 
3/57 #0.05263158
#### mean MPDin / mean MPDnn
length(Leafletsummary.distObsSum.ranm$t.MFD.p.value[which(Leafletsummary.distObsSum.ranm$t.MFD.p.value >= 0.05)]) # 53 NS
length(Leafletsummary.distObsSum.ranm$t.MFD.p.value[which(Leafletsummary.distObsSum.ranm$t.MFD.p.value <= 0.05)]) # 4 significant differneces between mean MPD
sigLeafletsummaryMPD <- Leafletsummary.distObsSum.ranm[which(Leafletsummary.distObsSum.ranm$t.MFD.p.value <= 0.05),] # 4 significant differneces between mean MPD
4/57 #0.07017544

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], Leafletsummary.distObsSumsig, by=0)) # 1/57 small islands sig MPD 0.01754386
nrow(merge(metadata[which(metadata$Size.cat == "med"),], Leafletsummary.distObsSumsig, by=0)) # 1 med islands sig MPD 0.01754386
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], Leafletsummary.distObsSumsig, by=0)) # 1 large isl sig MPD  0.01754386

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigLeafletsummaryMPD, by=0)) # 1 small islands sig MPD 0.01754386
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigLeafletsummaryMPD, by=0)) # 2/57 med islands sig MPD  0.03508772
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigLeafletsummaryMPD, by=0)) # 1 large isl sig MPD 0.01754386


### Leaf N
LeafNsummary.distObsSum <- lapply(com.list, function(x) functionObsSum(leafN.distObs[[x]])) #apply summary funciton across all communities...summary.MNNPD
names(LeafNsummary.distObsSum) <- colnames(SJcommNew)
LeafNsummary.distObsSum <- as.data.frame(do.call(cbind, LeafNsummary.distObsSum), stringsAsFactors=F)
colnames(LeafNsummary.distObsSum) <- colnames(SJcommNew)
LeafNsummary.distObsSum <-  as.data.frame(t(LeafNsummary.distObsSum), stringsAsFactors=F)
#write.csv(LeafNsummary.distObsSum, file="LeafNsummary.distObsSum.MFD.csv")
head(LeafNsummary.distObsSum)
LeafNsummary.distObsSum.ranm <- LeafNsummary.distObsSum[LeafNsummary.distObsSum$t.MNNFD.p.value != "NA",]

leafNtmp <- transform(merge(leafNSummary, LeafNsummary.distObsSum, by=0), row.names=Row.names,Row.names=NULL)
leafN <- (merge(leafNtmp, metadata, by=0))
head(leafN)
leafN[,"Significance"] <- ifelse(leafN[,"p.value.leafN"] <= 0.05, 1,
                                    ifelse(leafN[,"t.MNNFD.p.value"] <= 0.05, 2,
                                           ifelse(leafN[,"t.MFD.p.value"] <= 0.05, 3, 0)))

#### summarize significance data by size categories
leafN$Size.cat <- as.character(leafN$Size.cat)
leafN$Size.cat <- factor(leafN$Size.cat, levels=c("sm", "med", "lg"))

##### mean MMNPDi / mean MMNPDn 
length(LeafNsummary.distObsSum.ranm$t.MNNFD.p.value[which(LeafNsummary.distObsSum.ranm$t.MNNFD.p.value >= 0.05)]) # 44 NS
length(LeafNsummary.distObsSum.ranm$t.MNNFD.p.value[which(LeafNsummary.distObsSum.ranm$t.MNNFD.p.value <= 0.05)]) # 3 significant differneces between mean MNNPD 
sigLeafNsummary.distObsSum <- (LeafNsummary.distObsSum.ranm[which(LeafNsummary.distObsSum.ranm$t.MNNFD.p.value <= 0.05),]) # 3 significant differneces between mean MNNPD 
3/47 #0.06382979
#### mean MPDin / mean MPDnn
length(LeafNsummary.distObsSum.ranm$t.MFD.p.value[which(LeafNsummary.distObsSum.ranm$t.MFD.p.value >= 0.05)]) # 45 NS
length(LeafNsummary.distObsSum.ranm$t.MFD.p.value[which(LeafNsummary.distObsSum.ranm$t.MFD.p.value <= 0.05)]) # 2 significant differneces between mean MPD
LeafNsummaryMPD <- (LeafNsummary.distObsSum.ranm[which(LeafNsummary.distObsSum.ranm$t.MFD.p.value <= 0.05),]) # 2 significant differneces between mean MPD
2/47 #0.04255319

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigLeafNsummary.distObsSum, by=0)) # 0 small islands sig MPD 
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigLeafNsummary.distObsSum, by=0)) # 3/47 med islands sig MPD 0.06382979
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigLeafNsummary.distObsSum, by=0)) # 0 large isl sig MPD  

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], LeafNsummaryMPD, by=0)) # 0 small islands sig MPD 
nrow(merge(metadata[which(metadata$Size.cat == "med"),], LeafNsummaryMPD, by=0)) # 0 med islands sig MPD 
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], LeafNsummaryMPD, by=0)) # 2/47 large isl sig MPD 0.04255319


##################### TOTAL OBSERVED SUMMARY BAR PLOTS #########################
sigDNNS = cbind("Island"=phyloObsSum[,1], "Size.cat"=as.character(phyloObsSum[,"Size.cat"]), "Metric"=rep(x = "MNNPD", times = nrow(phyloObsSum)), "Significance"=ifelse(phyloObsSum[,"t.MNNPD.p.value"] <= 0.05, 1,0))
sigMPD = cbind("Island"=phyloObsSum[,1], "Size.cat"=as.character(phyloObsSum[,"Size.cat"]), "Metric"=rep(x = "MPD", times = nrow(phyloObsSum)), "Significance"=ifelse(phyloObsSum[,"t.MPD.p.value"] <= 0.05, 1,0))

sigSeed = cbind("Island"=seedmass[,1], "Size.cat"=as.character(seedmass[,"Size.cat"]), "Metric"=rep(x = "Seed Mass", times = nrow(seedmass)), "Significance"=ifelse(seedmass[,"p.value.seedMass"] <= 0.05, 1,0))
sigMNNFDSeed = cbind("Island"=seedmass[,1], "Size.cat"=as.character(seedmass[,"Size.cat"]), "Metric"=rep(x = "MNNFD Seed Mass", times = nrow(seedmass)), "Significance"=ifelse(seedmass[,"t.MNNFD.p.value"] <= 0.05, 1,0))
sigMFDSeed = cbind("Island"=seedmass[,1], "Size.cat"=as.character(seedmass[,"Size.cat"]), "Metric"=rep(x = "MFD Seed Mass", times = nrow(seedmass)), "Significance"=ifelse(seedmass[,"t.MFD.p.value"] <= 0.05, 1,0))

sigHeight = cbind("Island"=height[,1], "Size.cat"=as.character(height[,"Size.cat"]), "Metric"=rep(x = "Height", times = nrow(height)), "Significance"=ifelse(height[,"p.value.maxHeight"] <= 0.05, 1,0))
sigMNNFDheight = cbind("Island"=height[,1], "Size.cat"=as.character(height[,"Size.cat"]), "Metric"=rep(x = "MNNFD Height", times = nrow(height)), "Significance"=ifelse(height[,"t.MNNFD.p.value"] <= 0.05, 1,0))
sigMFDheight = cbind("Island"=height[,1], "Size.cat"=as.character(height[,"Size.cat"]), "Metric"=rep(x = "MFD Height", times = nrow(height)), "Significance"=ifelse(height[,"t.MFD.p.value"] <= 0.05, 1,0))

sigSLA = cbind("Island"=SLA[,1], "Size.cat"=as.character(SLA[,"Size.cat"]), "Metric"=rep(x = "SLA", times = nrow(SLA)), "Significance"=ifelse(SLA[,"p.value.sla"] <= 0.05, 1,0))
sigMNNFDSLA = cbind("Island"=SLA[,1], "Size.cat"=as.character(SLA[,"Size.cat"]), "Metric"=rep(x = "MNNFD SLA", times = nrow(SLA)), "Significance"=ifelse(SLA[,"t.MNNFD.p.value"] <= 0.05, 1,0))
sigMFDSLA = cbind("Island"=SLA[,1], "Size.cat"=as.character(SLA[,"Size.cat"]), "Metric"=rep(x = "MFD SLA", times = nrow(SLA)), "Significance"=ifelse(SLA[,"t.MFD.p.value"] <= 0.05, 1,0))

sigleaflet = cbind("Island"=leaflet[,1], "Size.cat"=as.character(leaflet[,"Size.cat"]), "Metric"=rep(x = "Leaf Size", times = nrow(leaflet)), "Significance"=ifelse(leaflet[,"p.value.leafletSize"] <= 0.05, 1,0))
sigMNNFDleaflet = cbind("Island"=leaflet[,1], "Size.cat"=as.character(leaflet[,"Size.cat"]), "Metric"=rep(x = "MNNFD Leaf Size", times = nrow(leaflet)), "Significance"=ifelse(leaflet[,"t.MNNFD.p.value"] <= 0.05, 1,0))
sigMFDleaflet = cbind("Island"=leaflet[,1], "Size.cat"=as.character(leaflet[,"Size.cat"]), "Metric"=rep(x = "MFD Leaf Size", times = nrow(leaflet)), "Significance"=ifelse(leaflet[,"t.MFD.p.value"] <= 0.05, 1,0))

sigleafN = cbind("IleafNnd"=leafN[,1], "Size.cat"=as.character(leafN[,"Size.cat"]), "Metric"=rep(x = "Leaf N", times = nrow(leafN)), "Significance"=ifelse(leafN[,"p.value.leafN"] <= 0.05, 1,0))
sigMNNFDleafN = cbind("IleafNnd"=leafN[,1], "Size.cat"=as.character(leafN[,"Size.cat"]), "Metric"=rep(x = "MNNFD Leaf N", times = nrow(leafN)), "Significance"=ifelse(leafN[,"t.MNNFD.p.value"] <= 0.05, 1,0))
sigMFDleafN = cbind("IleafNnd"=leafN[,1], "Size.cat"=as.character(leafN[,"Size.cat"]), "Metric"=rep(x = "MFD Leaf N", times = nrow(leafN)), "Significance"=ifelse(leafN[,"t.MFD.p.value"] <= 0.05, 1,0))

obs <- as.data.frame(rbind(sigSeed, sigMNNFDSeed, sigMFDSeed, 
                           sigHeight, sigMNNFDheight, sigMFDheight,
                           sigSLA, sigMNNFDSLA, sigMFDSLA,
                           sigleaflet, sigMNNFDleaflet, sigMFDleaflet,
                          sigleafN, sigMNNFDleafN, sigMFDleafN))

neworder <- c('Seed Mass', 'MNNFD Seed Mass', 'MFD Seed Mass', 
              'Height', 'MNNFD Height', 'MFD Height',
              'SLA', 'MNNFD SLA', 'MFD SLA',
              'Leaf Size', 'MNNFD Leaf Size', 'MFD Leaf Size',
              'Leaf N', 'MNNFD Leaf N', 'MFD Leaf N')

obs2 <- arrange(transform(obs, Metric=factor(Metric,levels=neworder)),Metric)

pdf("figs/plots/functionDiv/observed/obs.functionDistinct.SummaryBar.pdf")
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
dev.off()


