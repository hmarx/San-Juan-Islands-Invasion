###### Pretty Trees  #######
###### 20 Nov 2013 ########
###### Hannah E. Marx #######

## visualization of invasive and native species, significant nodes, fossil calibrated nodes
## uses confruifier to get nodes that are dated, and modifies color.plot.phylo (picante)

setwd("~/Dropbox/Work/TankLab/Projects/SanJuans/Manuscript/Drafts/Figs.v5/")

## FINAL TREE:  SJtreePL.bootstrap.tre
SJfinalTree <- read.tree("~/Dropbox/Work/TankLab/Projects/SanJuans/FINAL/7_Scaling/SJtreePL.bootstrap.Drop.tre")
trait <- read.csv("~/Dropbox/Work/TankLab/Projects/SanJuans/Manuscript/Drafts/Figs.v4/AppendixS2_SJtraits.newHeight.Final.numeric.csv", as.is=T, row.names=1) 
treedat <- treedata(SJfinalTree, trait)
tra <- as.data.frame(treedat$data)
head(tra)
dim(tra) # 366 8

### Plot Nodes
##### Match bootstrap node labels from ML tree to treePL scaled tree ####
#Modified from: http://treethinkers.blogspot.com/2008/10/labeling-trees-posterior-probability.html
#library(ape)

plot(ladderize(SJfinalTree, right=F), cex=.25)
bootTree <- read.tree(file="~/Dropbox/Work/TankLab/Projects/SanJuans/FINAL/6_Trees/Unweight/Concat/RAxML_bipartitions.align.concat0530.1000.unweight")
bootTree #367
#plot(bootTree, cex=.2, label.offset=.2)
#nodelabels(cex=.2,adj=c(1.2, -0.3), frame="n")
#bootTree$edge
#bootTree$edge[bootTree$tip.label=="Selaginella_wallacei"]
bootTree <- root(bootTree, outgroup="Selaginella_wallacei", resolve.root=T)
plot(ladderize(bootTree, right=F), cex=.25)
treedatBoot <- treedata(bootTree, trait)
SJfinalBoot <- drop.tip(bootTree, "Arctostaphylos_media") # to keep branch lenghts; treedata drops them for some reason....
SJfinalBoot #366

##The getAllSubTrees function below is a necessary subfunction that atomizes a tree into each individual subclade and was provided compliments of Luke Harmon.
getAllSubtrees<-function(phy, minSize=2) {
  res<-list()
  count=1
  ntip<-length(phy$tip.label)
  for(i in 1:phy$Nnode) {
    l<-tips(phy, ntip+i)
    bt<-match(phy$tip.label, l)
    if(sum(is.na(bt))==0) {
      st<-phy} 
    else st<-drop.tip(phy, phy$tip.label[is.na(bt)])
    if(length(st$tip.label)>=minSize) {
      res[[count]]<-st
      count<-count+1
    }
  }
  res
}

plotTreePLBoot <- function(treePL,bootTree) {
  getAllSubtrees(treePL)->treePLSub
  getAllSubtrees(bootTree)->bootSub
  #bootList<-matrix("<50",Nnode(treePL),1)
  bootList<-matrix("", Nnode(treePL),1)
  
  #The commands below compare all the subclades in the Bayes tree to all the subclades in the bootstrap tree, and vice versa, and identifies all those clades that are identical.
  for(i in 1:Nnode(treePL)) {
    for(j in 1:Nnode(bootTree)) {
      match(treePLSub[[i]]$tip.label[order(treePLSub[[i]]$tip.label)], bootSub[[j]]$tip.label[order(bootSub[[j]]$tip.label)])->shared
      match(bootSub[[j]]$tip.label[order(bootSub[[j]]$tip.label)], treePLSub[[i]]$tip.label[order(treePLSub[[i]]$tip.label)])->shared2
      if(sum(is.na(c(shared,shared2)))==0) {
        bootTree$node.label[j]->bootList[i]
      }}}
  treePLBS <- treePL
  treePLBS$node.label <- bootList
  plot(ladderize(treePLBS, right=F), cex=.10, lwd=0.1) #Plots your Bayesian consensus tree
  nodelabels(treePLBS$node.label, adj=c(1.2, -0.3), frame="n", cex=.2, font=2) #Adds bootstrap values.
  tree <- as.phylo(treePLBS)
  write.tree(tree, file="treePL.bootstrap.FINAL.tre") #SAVE as tree
  
}

plotTreePLBoot(treePL=SJfinalTree, bootTree=SJfinalBoot) 
# SJtreePL.bootstrap.FINAL.tre


#### Need to re-run Congruifier to get cal object (calibrated nodes)
genetree=read.nexus("~/Dropbox/Work/TankLab/Projects/SanJuans/FINAL/6_Trees/Unweight/Concat/RAxML_bipartitions.align.concat0530.1000.unweight.nex") #(saved as rooted)
genetree #367
genetree.drop= drop.tip(phy=genetree, tip="Arctostaphylos_media") #not in data
genetree.drop #366 
tax=read.csv(file="~/Dropbox/Work/TankLab/Programs/Congruify/congruify_hannah/fleshed_genera.csv", as.is=TRUE, row=1) ##"linkage table", From Jon (NESCent working group on plant rates and traits)
tips=sapply(genetree.drop$tip.label, function(x){
  unlist(strsplit(x,"_",fixed=TRUE))[1]
})
ll=match(tips, rownames(tax))
SJ_tax=tax[ll,]
rownames(SJ_tax)=names(tips)
SJ_tax=as.matrix(SJ_tax)
SJ_tax[is.na(SJ_tax)]=""

atol=read.tree("~/Dropbox/Work/TankLab/Programs/Congruify/congruify_hannah/out_dates.tre") #dataed reference tree, Soltis et al. 2011
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

#### Vector of congruified nodes
mrcaID=function(phy, cal){
  cal=as.matrix(cal)
  res=sapply(1:nrow(cal), function(idx){ ## loop over rows
    tips=cal[idx, c("taxonA", "taxonB")] ## fetch spanning taxa
    return(geiger:::.mrca(tips, phy)) ## MRCA of spanning taxa (node ID)
  })
  N=Ntip(phy)
  n=Nnode(phy)
  nn=integer(N+n) ## create empty vector of same length as branches in tree
  nn[res]=1 ## identify nodes that appear within calibrations
  nn=nn[-c(1:N)] ## exclude tip branches
  return(nn) ## return vector (ordered from first to last internal node in the tree)
}
vec=mrcaID(phy, cal) ## get vector for calibrated nodes
sum(vec)==nrow(cal) ## check on whether the function is working appropriately
# plot.phylo(phy, type="fan", show.tip=FALSE, edge.width=0.1)
## plot box at node only if calibrated
# nodelabels(text=NULL, cex=ifelse(vec==0, NA, 2), frame="n", bg="lightskyblue", col="lightgray", pch=22)


treePLboots <- read.tree("~/Dropbox/Work/TankLab/Projects/SanJuans/FINAL/7_Scaling/SJtreePL.bootstrap.FINAL.tre")  # Read in the treePL scaled tree with bootstraps
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

## Vector for Traits:
statusLabelsla <- character(length(treedat$phy$tip.label)) #We're going to make a new matrix to store the colors will use to label our tip taxa
names(statusLabelsla) <- rownames(tra)
length(tra$Status)
i <- as.data.frame(tra[tra$Status=="i",])
dim(i) #140
n <- as.data.frame(tra[tra$Status=="n",])
dim(n) #226

tra$seedMass
statusLabelseed <- character(length(treedat$phy$tip.label)) #We're going to make a new matrix to store the colors will use to label our tip taxa
names(statusLabelseed) <- rownames(tra)
statusLabelseed[] <- "white"
statusLabelseed[tra$seedMass!="NA"] <- "purple"

tra$maxHeight
statusLabelheight <- character(length(treedat$phy$tip.label)) #We're going to make a new matrix to store the colors will use to label our tip taxa
names(statusLabelheight) <- rownames(tra)
statusLabelheight[] <- "white"
statusLabelheight[tra$maxHeight!="NA"] <- "red2"

statusLabelsla[] <- "white"
#statusLabelsla[i$sla!="NA"] <- "red"
#statusLabelsla[n$sla!="NA"] <- "green"
statusLabelsla[tra$sla!="NA"] <- "gold1"

statusLabelleaf <- character(length(treedat$phy$tip.label)) #We're going to make a new matrix to store the colors will use to label our tip taxa
names(statusLabelleaf) <- rownames(tra)
statusLabelleaf[] <- "white"
statusLabelleaf[tra$leafletSize!="NA"] <- "green4"

tra$leafN
statusLabeln <- character(length(treedat$phy$tip.label)) #We're going to make a new matrix to store the colors will use to label our tip taxa
names(statusLabeln) <- rownames(tra)
statusLabeln[] <- "white"
statusLabeln[tra$leafN!="NA"] <- "blue"



###### Color taxa Native / Invasive: Code modified from poster.R:
#data
head(tra)
dim(tra)
names(tra)

tra2 <- cbind(rownames(tra), tra)
names(tra2) <- c("Sp.names", "Status","Habit","woodiness","sla","leafletSize","seedMass","leafN","maxHeight")  
names(tra2)

trsettings<-plot(treePLboots, typ="fan", cex=.3)

##Plot phylo with calibrated nodes IDed and nodes lables with taxonomy, dated with PATHd8
pdf(file="SJtreePLbootsStatusDropTraitsFinal.pdf", height=10, width=10) 
#jpeg(file="nodes3.jpeg")
color.plot.phylo3(par=par(mar = c(0, 0, 0, 1)), phylo=treePLboots, df=tra2, trait="Status",
                  taxa.names="Sp.names", col.names = c("magenta1", "green4"), cut.labs =c("Introduced", "Native"), 
                  leg.cex=.65) ##"turquoise4", "tomato"
dev.off()


color.plot.phylo3 <- function (par, phylo, df, trait, taxa.names, num.breaks = ifelse(is.factor(df[, trait]), length(levels(df[, trait])), 12), col.names = rainbow(ifelse(length(num.breaks) > 1, length(num.breaks) - 1, num.breaks)), cut.labs = NULL, leg.title = NULL, main = "", leg.cex = 1, tip.labs = NULL, ...) 
{
  init.par <- par 
  stopifnot(trait %in% names(df), taxa.names %in% names(df), 
            class(df) == "data.frame", class(phylo) == "phylo")
  len.tips <- length(phylo$tip.label)
  len.taxa <- length(df[, taxa.names])
  if (len.tips != len.taxa | sum(phylo$tip.label %in% df[,taxa.names]) != len.taxa) {
    stop("ERROR. Missing taxa in tree or data frame; # tips: ", 
         len.tips, "# taxa: ", len.taxa, "# tips in df: ", 
         sum(phylo$tip.label %in% df[, taxa.names]))
  }
  order <- match(phylo$tip.label, df[, taxa.names])
  ordered.trait <- df[trait][order, ]
  if (is.factor(ordered.trait)) {
    levs <- levels(ordered.trait)
    tip.color <- rep("black", times = len.taxa)
    tip.color <- col.names[match(ordered.trait, levs)]
  }
  else {
    tip.color = as.character(cut(ordered.trait, breaks = num.breaks, 
                                 labels = col.names))
    levs <- levels(cut(ordered.trait, breaks = num.breaks))
  }
  if (!is.null(tip.labs)) {
    phylo$tip.label <- df[tip.labs][order, ]
  }
  
  plot.phylo(phylo, label.offset=2, type="fan", cex=.28, plot=F)  ##increase = tighter to tips #0.15
  tiplabels(pch=16, cex=.3, col=statusLabeln[phylo$tip.label])
  par(new=TRUE)
  plot.phylo(phylo, label.offset=2, type="fan", cex=.26, plot=F) #.12
  tiplabels(pch=16, cex=.3, col=statusLabelsla[phylo$tip.label])
  par(new=TRUE)
  plot.phylo(phylo, label.offset=2, type="fan", cex=.24,  plot=F) #.09
  tiplabels(pch=16, cex=.3, col=statusLabelheight[phylo$tip.label])
  par(new=TRUE)
  plot.phylo(phylo, label.offset=2, type="fan", cex=.22,  plot=F) #.06
  tiplabels(pch=16, cex=.3, col=statusLabelleaf[phylo$tip.label])
  par(new=TRUE)
  plot.phylo(phylo, label.offset=2, type="fan", cex=.20, plot=F) #.03
  tiplabels(pch=16, cex=.3, col=statusLabelseed[phylo$tip.label])
  par(new=TRUE)
  plot(phylo, type="fan", label.offset=35, cex=.3, tip.color=tip.color, main=main,...)
  #plot.phylo(phylo, type="fan", label.offset=35, cex=.3)
  
  #plot.phylo(ladderize(phylo), type="fan", cex=phy.cex, label.offset=.02, tip.color=tip.color, main=main,...)
  nodelabels(text=NULL, cex=ifelse(vec==0, NA, 1), frame="n", col="black", pch=21, bg="white") # Plot congruified nodes
  nodelabels(pch = 21, cex = 0.3, col=p2, bg = p2) # Plot bootstrap support
  title(line = 0)
  if (is.null(cut.labs)) 
    cut.labs <- levs
  legend("topleft", cut.labs, fill = col.names, inset = 0.025, 
         title = leg.title, cex = leg.cex, bty = "n")
  co <- c("black", "slategray4", "black", "blue", "green4", "gold1", "red2", "purple")
  legend("bottomright", cut.labs, legend = c(expression(BS >= 95, 75 <= BS * " < 95", Dated), 
          "Leaf N", "Leaf Size", "SLA", "Max Height", "Seed Mass"), pch = c(21, 21, 1, 21,21,21,21,21), pt.bg = co, 
         inset = 0.0025, cex=leg.cex, bty = "n", col=c("white", "white", "black", "white","white","white","white","white")) #pt.cex=1.5, bty = "n"
  #co2 <- c("orange", "purple", "gold1", "green4", "red2")
  #legend("bottomleft", cut.labs, legend = c("Leaf N", "SLA", "Max Height", "Leaf Size", "Seed Mass"), pch = c(21, 21, 1), pt.bg = co2, 
         #inset = 0.025, cex=leg.cex) #pt.cex=1.5, bty = "n"
  on.exit(par(init.par))
}



#### Without Traits
#pdf("SJtreePLbootsStatusDropNEW.pdf")
color.plot.phylo3NoTrait(treePLboots, tra, "status", "Row.names", col.names = c("magenta1", "green4"), label.offset=.02, phy.cex=.3, cut.labs =c("Introduced", "Native"), leg.cex=.8) ##"turquoise4", "tomato"
#dev.off()

color.plot.phylo3NoTrait <- function (phylo, df, trait, taxa.names, label.offset, phy.cex=.3, num.breaks = ifelse(is.factor(df[, trait]), length(levels(df[, trait])), 12), col.names = rainbow(ifelse(length(num.breaks) > 1, length(num.breaks) - 1, num.breaks)), cut.labs = NULL, leg.title = NULL, main = trait, leg.cex = 1, tip.labs = NULL, ...) 
{
  init.par <- par(mar = c(0, 0, 1, 0))
  stopifnot(trait %in% names(df), taxa.names %in% names(df), 
            class(df) == "data.frame", class(phylo) == "phylo")
  len.tips <- length(phylo$tip.label)
  len.taxa <- length(df[, taxa.names])
  if (len.tips != len.taxa | sum(phylo$tip.label %in% df[,taxa.names]) != len.taxa) {
    stop("ERROR. Missing taxa in tree or data frame; # tips: ", 
         len.tips, "# taxa: ", len.taxa, "# tips in df: ", 
         sum(phylo$tip.label %in% df[, taxa.names]))
  }
  order <- match(phylo$tip.label, df[, taxa.names])
  ordered.trait <- df[trait][order, ]
  if (is.factor(ordered.trait)) {
    levs <- levels(ordered.trait)
    tip.color <- rep("black", times = len.taxa)
    tip.color <- col.names[match(ordered.trait, levs)]
  }
  else {
    tip.color = as.character(cut(ordered.trait, breaks = num.breaks, 
                                 labels = col.names))
    levs <- levels(cut(ordered.trait, breaks = num.breaks))
  }
  if (!is.null(tip.labs)) {
    phylo$tip.label <- df[tip.labs][order, ]
  }
  
  
  plot.phylo(ladderize(phylo), type="fan", cex=phy.cex, label.offset=label.offset, tip.color=tip.color, main=main,...) #label.offset=.02
  nodelabels(text=NULL, cex=ifelse(vec==0, NA, 1), frame="n", col="black", pch=21, bg="white") # Plot congruified nodes
  nodelabels(pch = 21, cex = 0.3, col=p2, bg = p2) # Plot bootstrap support
  title(line = 0)
  if (is.null(cut.labs)) 
    cut.labs <- levs
  legend("topleft", cut.labs, fill = col.names, inset = 0.025, 
         title = leg.title, cex = leg.cex, bty = "n")
  co <- c("black", "slategray4", "black")
  legend("bottomright", cut.labs, legend = c(expression(BS >= 95, 75 <= BS * " < 95", Dated)), pch = c(21, 21, 1, 21,21,21,21,21), pt.bg = co, 
         inset = 0.0025, cex=leg.cex, bty = "n", col=c("white", "white", "black", "white")) #pt.cex=1.5, bty = "n"
  #co2 <- c("orange", "purple", "gold1", "green4", "red2")
  #legend("bottomleft", cut.labs, legend = c("Leaf N", "SLA", "Max Height", "Leaf Size", "Seed Mass"), pch = c(21, 21, 1), pt.bg = co2, 
  #inset = 0.025, cex=leg.cex) #pt.cex=1.5, bty = "n"
  on.exit(par(init.par))
}










