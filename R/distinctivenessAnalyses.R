################### Final Analysis of San Juans Dataset

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
length(phyloObs) #74
names(phyloObs) <- colnames(SJcommNew)
names(phyloObs[1]) #"All_SanJuanIslands"
head(phyloObs[[1]])
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

# indicate islands with significant difference in observed means between status groups for each island
phyloObsSum[,"Significance"] <- ifelse(phyloObsSum[,"t.DNNS.p.value"] <= 0.05, 1,
                                            ifelse(phyloObsSum[,"t.MDNS.p.value"] <= 0.05, 2, 0)) 

#### summarize significance data by size categories
phyloObsSum$Size.cat <- as.character(phyloObsSum$Size.cat)
phyloObsSum$Size.cat <- factor(phyloObsSum$Size.cat, levels=c("sm", "med", "lg"))


sigDNNS = cbind("Island"=phyloObsSum[,1], "Size.cat"=as.character(phyloObsSum[,"Size.cat"]), "Metric"=rep(x = "DNNS", times = nrow(phyloObsSum)), "Significance"=ifelse(phyloObsSum[,"t.DNNS.p.value"] <= 0.05, 1,0))
sigMDNS = cbind("Island"=phyloObsSum[,1], "Size.cat"=as.character(phyloObsSum[,"Size.cat"]), "Metric"=rep(x = "MDNS", times = nrow(phyloObsSum)), "Significance"=ifelse(phyloObsSum[,"t.MDNS.p.value"] <= 0.05, 1,0))
obs.phylo <- as.data.frame(rbind(sigDNNS, sigMDNS))

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
length(phyloObsSum$t.DNNS.p.value[which(phyloObsSum$t.DNNS.p.value >= 0.05)]) # 34 islands NS
length(phyloObsSum$t.DNNS.p.value[which(phyloObsSum$t.DNNS.p.value <= 0.05)]) 
# 40 isalnds with significant differnece in mean betweenmean MMNPDi / mean MMNPDn 
phyloObsSum.DNNS.sig <- phyloObsSum[which(phyloObsSum$t.DNNS.p.value <= 0.05),]
phyloObsSum[which(phyloObsSum$t.DNNS.p.value <= 0.005),] # 3
40/74 # 0.5405405 

nrow(merge(metadataFULL[which(metadataFULL$Size.cat == "sm"),], phyloObsSum.DNNS.sig, by=1)) # 7/74 small islands sig DNNS 0.09459459
nrow(merge(metadataFULL[which(metadataFULL$Size.cat == "med"),], phyloObsSum.DNNS.sig, by=1)) # 20/74 medium isl sig DNNS 0.2702703
nrow(merge(metadataFULL[which(metadataFULL$Size.cat == "lg"),], phyloObsSum.DNNS.sig, by=1)) # 13/74 large isl sig DNNS 0.1756757

#### significance in observed difference bewteen mean MDNSin / mean MDNSnn
length(phyloObsSum$t.MDNS.p.value[which(phyloObsSum$t.MDNS.p.value >= 0.05)]) # 53 NS
length(phyloObsSum$t.MDNS.p.value[which(phyloObsSum$t.MDNS.p.value <= 0.05)]) # 21 significant differneces between mean MDNS
phyloObsSum.MDNS.sig <- phyloObsSum[which(phyloObsSum$t.MDNS.p.value <= 0.05),]
phyloObsSum[which(phyloObsSum$t.MDNS.p.value <= 0.005),] # 5
21/74 # 0.2837838

nrow(merge(metadataFULL[which(metadataFULL$Size.cat == "sm"),], phyloObsSum.MDNS.sig, by=1)) # 2/74 small islands sig MDNS 0.02702703
nrow(merge(metadataFULL[which(metadataFULL$Size.cat == "med"),], phyloObsSum.MDNS.sig, by=1)) # 4/74 med islands sig MDNS 0.05405405
nrow(merge(metadataFULL[which(metadataFULL$Size.cat == "lg"),], phyloObsSum.MDNS.sig, by=1)) # 15/74 large isl sig MDNS 0.2027027



########################  Plot observed phylogenetic distinctiveness for each island, ordered by increasing island size

### append metadata for island area to column in list of DNNS/MDNS for each island
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
r2.DNNS <- paste("R^2 = ", signif(summary(model)$r.squared, 3), sep="") #Explained variation / Total variation
p.DNNS <- paste("p-value = ", signif(anova(model)[[5]][1], 3), "***", sep="")

head(phyloObs_melt)
pdf("figs/plots/phyloDiv/observed/Observed.DNNS.islIncSize.pdf", width=20, height=10)
DNNS <- ggplot(phyloObs_melt, aes(x=reorder(factor(L2),value), y=as.numeric(as.character(MinDist.Nearest.native))), position=position_dodge(width=1))#, col=c("magenta1", "green3"))
DNNS <- DNNS + geom_boxplot(aes(fill = Species.Status), width = 1)
DNNS <- DNNS + geom_smooth(method="lm", se=T, color="black", aes(group=1))
DNNS <- DNNS + coord_cartesian(ylim=c(.5, 5000)) 
DNNS <- DNNS + scale_x_discrete("Island (increasing size)") #, breaks=seq(0, 80, 10)) 
DNNS <- DNNS + scale_y_log10("Log DNNS") #+ coord_fixed(ratio=4) 
DNNS <- DNNS + theme_bw() 
DNNS <- DNNS + scale_fill_manual(values=c("i"= "magenta1", "n"="green3"), labels=c("i"= "MMNPDi", "n"="MMNPDn")) ##breaks=rev(factor(SJnew$status)),
DNNS <- DNNS + guides(fill=guide_legend(title=""))
DNNS <- DNNS + theme(legend.position="top")
DNNS <- DNNS + theme(axis.text.x = element_text(angle = -45, hjust = 0))
DNNS <- DNNS + ggtitle("Phylogenetic Distance to the Nearest Native Species (DNNS)") + theme(plot.title=element_text(size=rel(1.5)))
DNNS <- DNNS + annotate("text", label=r2.DNNS, x=3, y=0.9, size=4) #y=max(as.numeric(as.character(SJ_NN_meltNEW$MinDist.Nearest.native))-20)
DNNS <- DNNS + annotate("text", label=p.DNNS, x=4, y=0.7, size=4) 
DNNS <- DNNS + theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
DNNS
dev.off()

four.islands.phyloObs_melt <- phyloObs_melt[phyloObs_melt$L2 %in% four.islands, ]
pdf("figs/plots/phyloDiv/observed/four.islands.phyloObs.pdf")
DNNS.four <- ggplot(four.islands.phyloObs_melt, aes(x=reorder(factor(L2),value), y=as.numeric(as.character(MinDist.Nearest.native))), position=position_dodge(width=1))#, col=c("magenta1", "green3"))
DNNS.four <- DNNS.four+ geom_boxplot(aes(fill = Species.Status), width = 1)
#DNNS.four+ geom_smooth(method="lm", se=T, color="black", aes(group=1))
DNNS.four <- DNNS.four+ coord_cartesian(ylim=c(.5, 5000))
DNNS.four <- DNNS.four+ scale_x_discrete("Island (increasing size)") #, breaks=seq(0, 80, 10)) 
DNNS.four <- DNNS.four+ scale_y_log10("Log DNNS") #+ coord_fixed(ratio=4) 
DNNS.four <- DNNS.four+ theme_bw() 
DNNS.four <- DNNS.four+ scale_fill_manual(values=c("i"= "magenta1", "n"="green3"), labels=c("i"= "MMNPDi", "n"="MMNPDn")) ##breaks=rev(factor(SJnew$status)),
DNNS.four <- DNNS.four+ guides(fill=guide_legend(title=""))
DNNS.four <- DNNS.four+ theme(legend.position="top")
DNNS.four <- DNNS.four+ theme(axis.text.x = element_text(angle = -45, hjust = 0))
DNNS.four <- DNNS.four+ ggtitle("Phylogenetic Distance to the Nearest Native Species (DNNS)") + theme(plot.title=element_text(size=rel(1.5)))
#DNNS.four+ annotate("text", label=r2.DNNS, x=5, y=0.09, size=4) #y=max(as.numeric(as.character(SJ_NN_meltNEW$MinDist.Nearest.native))-20)
#DNNS.four+ annotate("text", label=p.DNNS, x=5, y=0.05, size=4) 
DNNS.four <- DNNS.four+ theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
DNNS.four
dev.off()

## Regression of differnence in means to island size
modelMDNS <- lm(as.numeric(as.character(MeanDist.NativeCommunity)) ~ value , data=phyloObs_melt)
summary(modelMDNS)
anova(modelMDNS)
r2.MDNS <- paste("R^2 = ", signif(summary(modelMDNS)$r.squared, 3)) #Explained variation / Total variation
p.MDNS <- paste("p-value = ", signif(anova(modelMDNS)[[5]][1], 3), "**", sep="")

pdf("figs/plots/phyloDiv/observed/Observed.MDNS.islIncSize.pdf", width=20, height=10)
MDNS <- ggplot(phyloObs_melt, aes(x=reorder(factor(L2),value), y=as.numeric(as.character(MeanDist.NativeCommunity))), position=position_dodge(width=1))#, col=c("magenta1", "green3"))
MDNS <- MDNS + geom_boxplot(aes(fill = Species.Status), width = 1)
MDNS <- MDNS + geom_smooth(method="lm", se=T, color="black", aes(group=1))
MDNS <- MDNS + coord_cartesian(ylim=c(100, 2000))
MDNS <- MDNS + scale_x_discrete("Island (increasing size)") #, breaks=seq(0, 80, 10)) 
MDNS <- MDNS + scale_y_log10("Log MDNS")  
MDNS <- MDNS + theme_bw() 
MDNS <- MDNS + scale_fill_manual(values=c("i"= "magenta1", "n"="green3"), labels=c("i"= "MDNSNi", "n"="MDNSNn")) ##breaks=rev(factor(SJnew$status)),
MDNS <- MDNS + guides(fill=guide_legend(title=""))
MDNS <- MDNS + theme(legend.position="top")
MDNS <- MDNS + theme(axis.text.x = element_text(angle = -45, hjust = 0))
MDNS <- MDNS + ggtitle("Mean Phylogenetic Distance to Native Community (MDNS)") + theme(plot.title=element_text(size=rel(1.5)))
MDNS <- MDNS + annotate("text", label=r2.MDNS, x=3, y=130, size=4) #y=max(as.numeric(as.character(SJ_NN_meltNEW$MinDist.Nearest.native))-20)
MDNS <- MDNS + annotate("text", label=p.MDNS, x=3.5, y=120, size=4) 
MDNS <- MDNS + theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
MDNS
dev.off()

pdf("figs/plots/phyloDiv/observed/four.islands.MDNS.islIncSize.pdf")
MDNS.four <- ggplot(four.islands.phyloObs_melt, aes(x=reorder(factor(L2),value), y=as.numeric(as.character(MeanDist.NativeCommunity))), position=position_dodge(width=1))#, col=c("magenta1", "green3"))
MDNS.four <- MDNS.four + geom_boxplot(aes(fill = Species.Status), width = 1)
MDNS.four <- MDNS.four  + coord_cartesian(ylim=c(100, 2000))
#MDNS.four <- MDNS.four  + geom_smooth(method="lm", se=T, color="black", aes(group=1))
MDNS.four <- MDNS.four  + scale_x_discrete("Island (increasing size)") #, breaks=seq(0, 80, 10)) 
MDNS.four <- MDNS.four  + scale_y_log10("Log MDNSN")  
MDNS.four <- MDNS.four  + theme_bw() 
MDNS.four <- MDNS.four  + scale_fill_manual(values=c("i"= "magenta1", "n"="green3"), labels=c("i"= "MDNSi", "n"="MDNSn")) ##breaks=rev(factor(SJnew$status)),
MDNS.four <- MDNS.four  + guides(fill=guide_legend(title=""))
MDNS.four <- MDNS.four  + theme(legend.position="top")
MDNS.four <- MDNS.four  + theme(axis.text.x = element_text(angle = -45, hjust = 0))
MDNS.four <- MDNS.four  + ggtitle("Mean Phylogenetic Distance to Native Community (MDNS)") + theme(plot.title=element_text(size=rel(1.5)))
#MDNS.four <- MDNS.four  + annotate("text", label=r2.MDNS, x=5, y=100, size=4) #y=max(as.numeric(as.character(SJ_NN_meltNEW$MinDist.Nearest.native))-20)
#MDNS.four <- MDNS.four  + annotate("text", label=p.MDNS, x=5, y=90, size=4) 
MDNS.four <- MDNS.four  + theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
MDNS.four 
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
#write.csv(SJtraitSummary, file="figs/plots/SJtraitSummary.csv")  ## edited in excel to make Table1_SJtraitSummary.csv


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
######### Calculate observed difference in trait values to nearest native (NNFD) and to mean value odnative community (MFD) #######################

############# Nearest Neighbor Trait Difference (NNFD) ##################################################################
## Entire archipelago
seedMass_SJ <- functionDistinct(output=phyloObs[[1]], traits=SJtraitLog, traitname="seedMass")
maxHeight_SJ <- functionDistinct(output=phyloObs[[1]], traits=SJtraitLog, traitname="maxHeight")
SLA_SJ <- functionDistinct(output=phyloObs[[1]], traits=SJtraitLog, traitname="sla") #one community, one trait
leafletSize_SJ <- functionDistinct(output=phyloObs[[1]], traits=SJtraitLog, traitname="leafletSize")
leafN_SJ <- functionDistinct(output=phyloObs[[1]], traits=SJtraitLog, traitname="leafN")
a2 <- plot.functionDistinct.Obs(NNFDoutput=seedMass_SJ, islandname="", traitname="seedMass", metric = 2)
b2 <- plot.functionDistinct.Obs(NNFDoutput=maxHeight_SJ, islandname="", traitname="maxHeight", metric = 2)
c2 <- plot.functionDistinct.Obs(NNFDoutput=SLA_SJ, islandname="", traitname="sla", metric = 2)
d2 <- plot.functionDistinct.Obs(NNFDoutput=leafletSize_SJ, islandname="", traitname="leafletSize", metric = 2)
e2 <- plot.functionDistinct.Obs(NNFDoutput=leafN_SJ, islandname="", traitname="leafN", metric = 2)
grid.arrange(a2, b2, c2, d2, e2, ncol=6, main=textGrob(vjust = .8, hjust=.5,"NNFD \nAll San Juan Islands",gp=gpar(cex=1.5,font=2)))

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
SeedMass.distObsSum <- lapply(com.list, function(x) functionObsSum(seedmass.distObs[[x]])) #apply summary funciton across all communities...summary.DNNS
names(SeedMass.distObsSum) <- colnames(SJcommNew)
SeedMass.distObsSum <- as.data.frame(do.call(cbind, SeedMass.distObsSum), stringsAsFactors=F)
colnames(SeedMass.distObsSum) <- colnames(SJcommNew)
SeedMass.distObsSum <-  as.data.frame(t(SeedMass.distObsSum), stringsAsFactors=F)
#write.csv(SeedMass.distObsSum, file="SeedMass.distObsSum.MFD.csv")
head(SeedMass.distObsSum)
SeedMass.distObsSum.ranm <- SeedMass.distObsSum[SeedMass.distObsSum$t.NNFD.p.value != "NA",]

seedtmp <- transform(merge(seedmassSummary, SeedMass.distObsSum, by=0), row.names=Row.names,Row.names=NULL)
seedmass <- (merge(seedtmp, metadata, by=0))
dim(seedtmp)
seedmass[,"Significance"] <- ifelse(seedmass[,"p.value.seedMass"] <= 0.05, 1,
                                    ifelse(seedmass[,"t.NNFD.p.value"] <= 0.05, 2,
                                           ifelse(seedmass[,"t.MFD.p.value"] <= 0.05, 3, 0)))
head(seedmass)
#### summarize significance data by size categories
seedmass$Size.cat <- as.character(seedmass$Size.cat)
seedmass$Size.cat <- factor(seedmass$Size.cat, levels=c("sm", "med", "lg"))


##### mean MMNPDi / mean MMNPDn 
length(SeedMass.distObsSum.ranm$t.NNFD.p.value[which(SeedMass.distObsSum.ranm$t.NNFD.p.value >= 0.05)]) # 72 NS
length(SeedMass.distObsSum.ranm$t.NNFD.p.value[which(SeedMass.distObsSum.ranm$t.NNFD.p.value <= 0.05)]) # 1 significant differneces between mean DNNS 
sigSeedMass.distObsSum <- (SeedMass.distObsSum.ranm[which(SeedMass.distObsSum.ranm$t.NNFD.p.value <= 0.05), ]) # 1 significant differneces between mean DNNS 
1/73 #0.01369863
#### mean MDNSin / mean MDNSnn
length(SeedMass.distObsSum.ranm$t.MFD.p.value[which(SeedMass.distObsSum.ranm$t.MFD.p.value >= 0.05)]) # 62 NS
length(SeedMass.distObsSum.ranm$t.MFD.p.value[which(SeedMass.distObsSum.ranm$t.MFD.p.value <= 0.05)]) # 11 significant differneces between mean MDNS
sigSeedMasssummaryMDNS <- (SeedMass.distObsSum.ranm[which(SeedMass.distObsSum.ranm$t.MFD.p.value <= 0.05),]) # 11 significant differneces between mean MDNS
11/73 #0.1506849

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigSeedMass.distObsSum, by=0)) # 0 small islands sig MDNS
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigSeedMass.distObsSum, by=0)) # 1/73 med islands sig MDNS 0.01369863
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigSeedMass.distObsSum, by=0)) # 0 large isl sig MDNS

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigSeedMasssummaryMDNS, by=0)) #2/73 small islands sig MDNS 0.02739726
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigSeedMasssummaryMDNS, by=0)) #4/73 med islands sig MDNS 0.05479452
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigSeedMasssummaryMDNS, by=0)) # 5/73 large isl sig MDNS 0.06849315


### Maximum height
Heightsummary.distObsSum <- lapply(com.list, function(x) functionObsSum(maxheight.distObs[[x]])) #apply summary funciton across all communities...summary.DNNS
names(Heightsummary.distObsSum) <- colnames(SJcommNew)
Heightsummary.distObsSum <- as.data.frame(do.call(cbind, Heightsummary.distObsSum), stringsAsFactors=F)
colnames(Heightsummary.distObsSum) <- colnames(SJcommNew)
Heightsummary.distObsSum <-  as.data.frame(t(Heightsummary.distObsSum), stringsAsFactors=F)
#write.csv(Heightsummary.distObsSum, file="Heightsummary.distObsSum.MFD.csv")
head(Heightsummary.distObsSum)
Heightsummary.distObsSum.ranm <- Heightsummary.distObsSum[Heightsummary.distObsSum$t.NNFD.p.value != "NA",]

heighttmp <- transform(merge(maxheightSummary, Heightsummary.distObsSum, by=0), row.names=Row.names,Row.names=NULL)
height <- (merge(heighttmp, metadata, by=0))
head(height)
height[,"Significance"] <- ifelse(height[,"p.value.maxHeight"] <= 0.05, 1,
                                  ifelse(height[,"t.NNFD.p.value"] <= 0.05, 2,
                                         ifelse(height[,"t.MFD.p.value"] <= 0.05, 3, 0)))

#### summarize significance data by size categories
height$Size.cat <- as.character(height$Size.cat)
height$Size.cat <- factor(height$Size.cat, levels=c("sm", "med", "lg"))


##### mean MMNPDi / mean MMNPDn 
length(Heightsummary.distObsSum.ranm$t.NNFD.p.value[which(Heightsummary.distObsSum.ranm$t.NNFD.p.value >= 0.05)]) # 51 NS
length(Heightsummary.distObsSum.ranm$t.NNFD.p.value[which(Heightsummary.distObsSum.ranm$t.NNFD.p.value <= 0.05)]) # 20 significant differneces between mean DNNS 
sigHeightsummary.distObsSum <- (Heightsummary.distObsSum.ranm[which(Heightsummary.distObsSum.ranm$t.NNFD.p.value <= 0.05),]) # 20 significant differneces between mean DNNS 
20/71 #0.2816901
#### mean MDNSin / mean MDNSnn
length(Heightsummary.distObsSum.ranm$t.MFD.p.value[which(Heightsummary.distObsSum.ranm$t.MFD.p.value >= 0.05)]) # 54 NS
length(Heightsummary.distObsSum.ranm$t.MFD.p.value[which(Heightsummary.distObsSum.ranm$t.MFD.p.value <= 0.05)]) # 17 significant differneces between mean MDNS
sigHeightsummaryMDNS <- (Heightsummary.distObsSum.ranm[which(Heightsummary.distObsSum.ranm$t.MFD.p.value <= 0.05),]) # 17 significant differneces between mean MDNS
17/71 #0.2394366

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigHeightsummary.distObsSum, by=0)) # 4/71 small islands sig MDNS 0.05633803
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigHeightsummary.distObsSum, by=0)) # 12/71 med islands sig MDNS 0.1690141
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigHeightsummary.distObsSum, by=0)) # 4/71 large isl sig MDNS  0.05633803

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigHeightsummaryMDNS, by=0)) # 3/71 small islands sig MDNS 0.04225352
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigHeightsummaryMDNS, by=0)) # 8/71 med islands sig MDNS 0.1126761
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigHeightsummaryMDNS, by=0)) # 6/71 large isl sig MDNS 0.08450704

### SLA
SLAsummary.distObsSum  <- lapply(com.list, function(x) functionObsSum(sla.distObs[[x]])) #apply summary funciton across all communities...summary.DNNS
names(SLAsummary.distObsSum) <- colnames(SJcommNew)
SLAsummary.distObsSum <- as.data.frame(do.call(cbind, SLAsummary.distObsSum), stringsAsFactors=F)
colnames(SLAsummary.distObsSum) <- colnames(SJcommNew)
SLAsummary.distObsSum <-  as.data.frame(t(SLAsummary.distObsSum), stringsAsFactors=F)
#write.csv(SLAsummary.distObsSum, file="SLAsummary.distObsSum.MFD.csv")
head(SLAsummary.distObsSum)
SLAsummary.distObsSum.ranm <- SLAsummary.distObsSum[SLAsummary.distObsSum$meanNNFDnatives != "NA",]

SLAtmp <- transform(merge(SLAtraitSummary, SLAsummary.distObsSum, by=0), row.names=Row.names, Row.names=NULL)
SLA <- (merge(SLAtmp, metadata, by=0))

SLA[,"Significance"] <- ifelse(SLA[,"p.value.sla"] <= 0.05, 1,
                              ifelse(SLA[,"t.NNFD.p.value"] <= 0.05, 2,
                               ifelse(SLA[,"t.MFD.p.value"] <= 0.05, 3, 0)))

#### summarize significance data by size categories
SLA$Size.cat <- as.character(SLA$Size.cat)
SLA$Size.cat <- factor(SLA$Size.cat, levels=c("sm", "med", "lg"))

##### mean MMNPDi / mean MMNPDn 
length(SLAsummary.distObsSum.ranm$t.NNFD.p.value[which(SLAsummary.distObsSum.ranm$t.NNFD.p.value >= 0.05)]) # 51 NS
length(SLAsummary.distObsSum.ranm$t.NNFD.p.value[which(SLAsummary.distObsSum.ranm$t.NNFD.p.value <= 0.05)]) # 10 
sigSLAsummary.distObsSum <- (SLAsummary.distObsSum.ranm[which(SLAsummary.distObsSum.ranm$t.NNFD.p.value <= 0.05),]) # 10 
10/61 #0.1639344
#### mean MDNSin / mean MDNSnn
length(SLAsummary.distObsSum.ranm$t.MFD.p.value[which(SLAsummary.distObsSum.ranm$t.MFD.p.value >= 0.05)]) # 57 NS
length(SLAsummary.distObsSum.ranm$t.MFD.p.value[which(SLAsummary.distObsSum.ranm$t.MFD.p.value <= 0.05)]) # 4 significant differneces between mean MDNS
sigSLAsummaryMDNS <- SLAsummary.distObsSum.ranm[which(SLAsummary.distObsSum.ranm$t.MFD.p.value <= 0.05),] # 4 significant differneces between mean MDNS
4/61 #0.06557377

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigSLAsummary.distObsSum, by=0)) # 1/61 small islands sig MDNS 0.01639344
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigSLAsummary.distObsSum, by=0)) # 6/61 med islands sig MDNS 0.09836066
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigSLAsummary.distObsSum, by=0)) # 3/61 large isl sig MDNS 0.04918033

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigSLAsummaryMDNS, by=0)) # 2/61 small islands sig MDNS 0.03278689
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigSLAsummaryMDNS, by=0)) # 1/61 med islands sig MDNS 0.01639344
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigSLAsummaryMDNS, by=0)) # 1/61 large isl sig MDNS 0.01639344


### Leaflet 
Leafletsummary.distObsSum <- lapply(com.list, function(x) functionObsSum(leaflet.distObs[[x]])) #apply summary funciton across all communities...summary.DNNS
names(Leafletsummary.distObsSum) <- colnames(SJcommNew)
Leafletsummary.distObsSum <- as.data.frame(do.call(cbind, Leafletsummary.distObsSum), stringsAsFactors=F)
colnames(Leafletsummary.distObsSum) <- colnames(SJcommNew)
Leafletsummary.distObsSum <-  as.data.frame(t(Leafletsummary.distObsSum), stringsAsFactors=F)
#write.csv(Leafletsummary.distObsSum, file="Leafletsummary.distObsSum.MFD.csv")
head(Leafletsummary.distObsSum)
Leafletsummary.distObsSum.ranm <- Leafletsummary.distObsSum[Leafletsummary.distObsSum$t.NNFD.p.value != "NA",]

leaflettmp <-transform(merge(leafletSummary, Leafletsummary.distObsSum, by=0), row.names=Row.names,Row.names=NULL)
leaflet <- (merge(leaflettmp, metadata,by=0))
head(leaflet)
leaflet[,"Significance"] <- ifelse(leaflet[,"p.value.leafletSize"] <= 0.05, 1,
                               ifelse(leaflet[,"t.NNFD.p.value"] <= 0.05, 2,
                                      ifelse(leaflet[,"t.MFD.p.value"] <= 0.05, 3, 0)))

#### summarize significance data by size categories
leaflet$Size.cat <- as.character(leaflet$Size.cat)
leaflet$Size.cat <- factor(leaflet$Size.cat, levels=c("sm", "med", "lg"))

##### mean MMNPDi / mean MMNPDn 
length(Leafletsummary.distObsSum.ranm$t.NNFD.p.value[which(Leafletsummary.distObsSum.ranm$t.NNFD.p.value >= 0.05)]) # 54 NS
length(Leafletsummary.distObsSum.ranm$t.NNFD.p.value[which(Leafletsummary.distObsSum.ranm$t.NNFD.p.value <= 0.05)]) # 3 significant differneces between mean DNNS 
Leafletsummary.distObsSumsig <- (Leafletsummary.distObsSum[which(Leafletsummary.distObsSum.ranm$t.NNFD.p.value <= 0.05),]) # 3 significant differneces between mean DNNS 
3/57 #0.05263158
#### mean MDNSin / mean MDNSnn
length(Leafletsummary.distObsSum.ranm$t.MFD.p.value[which(Leafletsummary.distObsSum.ranm$t.MFD.p.value >= 0.05)]) # 53 NS
length(Leafletsummary.distObsSum.ranm$t.MFD.p.value[which(Leafletsummary.distObsSum.ranm$t.MFD.p.value <= 0.05)]) # 4 significant differneces between mean MDNS
sigLeafletsummaryMDNS <- Leafletsummary.distObsSum.ranm[which(Leafletsummary.distObsSum.ranm$t.MFD.p.value <= 0.05),] # 4 significant differneces between mean MDNS
4/57 #0.07017544

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], Leafletsummary.distObsSumsig, by=0)) # 1/57 small islands sig MDNS 0.01754386
nrow(merge(metadata[which(metadata$Size.cat == "med"),], Leafletsummary.distObsSumsig, by=0)) # 1 med islands sig MDNS 0.01754386
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], Leafletsummary.distObsSumsig, by=0)) # 1 large isl sig MDNS  0.01754386

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigLeafletsummaryMDNS, by=0)) # 1 small islands sig MDNS 0.01754386
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigLeafletsummaryMDNS, by=0)) # 2/57 med islands sig MDNS  0.03508772
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigLeafletsummaryMDNS, by=0)) # 1 large isl sig MDNS 0.01754386


### Leaf N
LeafNsummary.distObsSum <- lapply(com.list, function(x) functionObsSum(leafN.distObs[[x]])) #apply summary funciton across all communities...summary.DNNS
names(LeafNsummary.distObsSum) <- colnames(SJcommNew)
LeafNsummary.distObsSum <- as.data.frame(do.call(cbind, LeafNsummary.distObsSum), stringsAsFactors=F)
colnames(LeafNsummary.distObsSum) <- colnames(SJcommNew)
LeafNsummary.distObsSum <-  as.data.frame(t(LeafNsummary.distObsSum), stringsAsFactors=F)
#write.csv(LeafNsummary.distObsSum, file="LeafNsummary.distObsSum.MFD.csv")
head(LeafNsummary.distObsSum)
LeafNsummary.distObsSum.ranm <- LeafNsummary.distObsSum[LeafNsummary.distObsSum$t.NNFD.p.value != "NA",]

leafNtmp <- transform(merge(leafNSummary, LeafNsummary.distObsSum, by=0), row.names=Row.names,Row.names=NULL)
leafN <- (merge(leafNtmp, metadata, by=0))
head(leafN)
leafN[,"Significance"] <- ifelse(leafN[,"p.value.leafN"] <= 0.05, 1,
                                    ifelse(leafN[,"t.NNFD.p.value"] <= 0.05, 2,
                                           ifelse(leafN[,"t.MFD.p.value"] <= 0.05, 3, 0)))

#### summarize significance data by size categories
leafN$Size.cat <- as.character(leafN$Size.cat)
leafN$Size.cat <- factor(leafN$Size.cat, levels=c("sm", "med", "lg"))

##### mean MMNPDi / mean MMNPDn 
length(LeafNsummary.distObsSum.ranm$t.NNFD.p.value[which(LeafNsummary.distObsSum.ranm$t.NNFD.p.value >= 0.05)]) # 44 NS
length(LeafNsummary.distObsSum.ranm$t.NNFD.p.value[which(LeafNsummary.distObsSum.ranm$t.NNFD.p.value <= 0.05)]) # 3 significant differneces between mean DNNS 
sigLeafNsummary.distObsSum <- (LeafNsummary.distObsSum.ranm[which(LeafNsummary.distObsSum.ranm$t.NNFD.p.value <= 0.05),]) # 3 significant differneces between mean DNNS 
3/47 #0.06382979
#### mean MDNSin / mean MDNSnn
length(LeafNsummary.distObsSum.ranm$t.MFD.p.value[which(LeafNsummary.distObsSum.ranm$t.MFD.p.value >= 0.05)]) # 45 NS
length(LeafNsummary.distObsSum.ranm$t.MFD.p.value[which(LeafNsummary.distObsSum.ranm$t.MFD.p.value <= 0.05)]) # 2 significant differneces between mean MDNS
LeafNsummaryMDNS <- (LeafNsummary.distObsSum.ranm[which(LeafNsummary.distObsSum.ranm$t.MFD.p.value <= 0.05),]) # 2 significant differneces between mean MDNS
2/47 #0.04255319

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], sigLeafNsummary.distObsSum, by=0)) # 0 small islands sig MDNS 
nrow(merge(metadata[which(metadata$Size.cat == "med"),], sigLeafNsummary.distObsSum, by=0)) # 3/47 med islands sig MDNS 0.06382979
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], sigLeafNsummary.distObsSum, by=0)) # 0 large isl sig MDNS  

nrow(merge(metadata[which(metadata$Size.cat == "sm"),], LeafNsummaryMDNS, by=0)) # 0 small islands sig MDNS 
nrow(merge(metadata[which(metadata$Size.cat == "med"),], LeafNsummaryMDNS, by=0)) # 0 med islands sig MDNS 
nrow(merge(metadata[which(metadata$Size.cat == "lg"),], LeafNsummaryMDNS, by=0)) # 2/47 large isl sig MDNS 0.04255319


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


