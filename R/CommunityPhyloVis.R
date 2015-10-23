
###### Pretty Trees  #######
###### 20 Nov 2013 ########
###### Hannah E. Marx #######

source("R/treeFunctions.R")
source("R/trait.plot.colorTips.R")

## visualization of invasive and native species, significant nodes, fossil calibrated nodes
## uses confruifier to get nodes that are dated, and modifies color.plot.phylo (picante)


SJfinalTree <- read.tree("data/SJtreePL.bootstrap.tre") ## FINAL TREE (after prepPipeline.R)
trait <- read.csv("data/SJtraits.csv", as.is=T, row.names=1) 
treedat <- treedata(SJfinalTree, trait)
tra <- as.data.frame(treedat$data)
head(tra)
dim(tra) # 366 8

#### Plot Nodes
# Match bootstrap node labels from ML tree to treePL scaled tree
# Modified from: http://treethinkers.blogspot.com/2008/10/labeling-trees-posterior-probability.html

plot(ladderize(SJfinalTree, right=F), cex=.25)
SJfinalTree <- drop.tip(SJfinalTree, "Arctostaphylos_media")
SJfinalTree
# Read in tree with nodelabels == bootstrap values
bootTree <- read.tree(file="output/6_Trees/Concatenated/RAxML_bipartitions.align.concat0530.1000.unweight")
bootTree # 367 tips
bootTree <- root(bootTree, outgroup="Selaginella_wallacei", resolve.root=T)
plot(ladderize(bootTree, right=F), cex=.25)
treedatBoot <- treedata(bootTree, trait)
SJfinalBoot <- drop.tip(bootTree, "Arctostaphylos_media") # to keep branch lenghts; treedata drops them for some reason....
SJfinalBoot #366
plotTreePLBoot(treePL=SJfinalTree, bootTree=SJfinalBoot) # writes SJtreePL.bootstrap.FINAL.tre

#### Need to re-run Congruifier to get cal object (calibrated nodes)
genetree=read.nexus("output/6_Trees/Concatenated/RAxML_bipartitions.align.concat0530.1000.unweight.nex") #(saved as rooted)
genetree #367
genetree.drop= drop.tip(phy=genetree, tip="Arctostaphylos_media") #not in data
genetree.drop #366 
tax=read.csv(file="output/7_Scaling/Congruify/fleshed_genera.csv", as.is=TRUE, row=1) ##"linkage table", From Jon (NESCent working group on plant rates and traits)
tips=sapply(genetree.drop$tip.label, function(x){
  unlist(strsplit(x,"_",fixed=TRUE))[1]
})
ll=match(tips, rownames(tax))
SJ_tax=tax[ll,]
rownames(SJ_tax)=names(tips)
SJ_tax=as.matrix(SJ_tax)
SJ_tax[is.na(SJ_tax)]=""

atol=read.tree("output/7_Scaling/Congruify/out_dates.tre") #dataed reference tree, Soltis et al. 2011
ftax=tax[match(atol$tip.label, rownames(tax)),]
ftax[,2]="Spermatophyta"
fatol=subset(atol, ftax, "family")

phy=genetree.drop
swaptips=paste(1:length(tips),tips,sep="_")
phy$tip.label=swaptips
tax=SJ_tax
rownames(tax)=swaptips
res=congruify.phylo(fatol, phy, tax, tol=0, scale="PATHd8") # need to use PATHd8 to get res$phy
res
cal=res$calibrations
vec=mrcaID(phy, cal) ## get vector for calibrated nodes
sum(vec)==nrow(cal) ## check on whether the function is working appropriately

treePLboots <- read.tree("figs/trees/treePL.bootstrap.FINAL.tre")  # Read in the treePL scaled tree with bootstraps
treePLboots # 366, nodelabes = boots, rooted
# Vector for bootstrap support
p2 <- character(length(treePLboots$node.label)) # create a vector of colors that correspond to different support ranges 
p2[] <- "#0000ff00" ## Transparent color: "#RRGGBBAA" and the AA portion is the opacity/trasparency.
p2[treePLboots$node.label >=  95] <- "black"
p2[treePLboots$node.label ==  100] <- "black"
p2[treePLboots$node.label < 95 & treePLboots$node.label >= 75] <- "slategray4"
#p2[treePLboots$node.label < 95 & treePLboots$node.label >= 75] <- "slategray"
#p2[treePLboots$node.label < 75] <- "red1"
p2[treePLboots$node.label == ""] <- "black" # node = 353 = 100
p2[[1]] <- "#0000ff00" # "NA" => blank
p2[[2]] <- "#0000ff00"  # treePLboots$node.label[2]; 96 => blank
p2[[3]] <- "black" # 97 => 96
p2[[348]] <- "black" # 100 => 97
p2 

# See bootstraps nodes on treePL tree
plot.phylo(treePLboots, show.node.label=F, type="fan", cex=.2)
nodelabels(as.numeric(1:length(treePLboots$node.label)), cex=.5)

tra2 <- cbind(rownames(tra), tra)
names(tra2) <- c("Sp.names", "Status", "seedMass","maxHeight", "sla","leafletSize","leafN")  
names(tra2)

df <- (tra[,2:ncol(tra)])
df <- as.data.frame(tra[,2:ncol(tra)], stringsAsFactors=F)
str(df)
df[,1:ncol(df)]<- sapply(df[,1:ncol(df)],as.numeric) 
df[!is.na(df)] <- 1
df[is.na(df)] <- 0
#df <- cbind(tra[,1], df)
names(df) <- c( "seedMass","maxHeight", "sla","leafletSize","leafN")  
df

df2 <- as.data.frame(tra[2:ncol(tra)], stringsAsFactors=F)
df2[,1:ncol(df2)]<- sapply(df2[,1:ncol(df2)],as.numeric) 
head(df2)
df2[!is.na(df2)] <- 1
df2[is.na(df2)] <- 0
df2 <- cbind(as.numeric(tra[,1]), df2)
names(df2) <- c( "Status", "seedMass","maxHeight", "sla","leafletSize","leafN")  
df2

##Plot phylo with calibrated nodes IDed and nodes lables with taxonomy, dated with PATHd8
#pdf(file="figs/trees/SJtreeNodesTraitsStatus.pdf", height=10, width=10) 
trait.plot.colorTip(tree = treePLboots, dat = df, cols = list(seedMass = c("white", "purple"),
                                                        maxHeight = c("white", "red2"),
                                                        sla = c("white", "gold1"),
                                                        leafletSize = c("white", "green4"),
                                                        leafN = c("white", "blue")), 
                    datTr = tra2, trait = "Status", taxa.names = "Sp.names", col.names = c("magenta1", "green4"), 
                    cex.lab = 0.3, font.lab = 1)
#dev.off()

##Plot phylo with calibrated nodes IDed and nodes lables with taxonomy, dated with PATHd8
#pdf(file="figs/trees/SJtreeNodesTraitsStatusNoTips2.pdf", height=10, width=10) 
trait.plot.colorTip(tree = treePLboots, dat = df2, cols = list(Status = c("magenta1", "green4"),
                                                         seedMass = c("white", "yellow"),
                                                         maxHeight = c("white", "blue"),
                                                         sla = c("white", "orange"),
                                                         leafletSize = c("white", "lightskyblue3"),
                                                         leafN = c("white", "purple")), 
                    datTr = tra2, trait = "Status", taxa.names = "Sp.names", col.names = c("white", "white"), 
                    cex.lab = 0.3, font.lab = 1)
#dev.off()

