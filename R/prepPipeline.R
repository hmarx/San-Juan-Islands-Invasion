#########################################################################################################
###########################  6 Nov 2014     #############################################################
###########################  Hannah E. Marx #############################################################

### Pre-processing for Phylogeny of the flora of the San Jaun Islands (data/SJtreePL.bootstrap.tre) 
### For analyses of Phylogenetic and Functional Distinctiveness, just load final community phylogeny 

######################### 2_SpeciesList:  ##############################################################

# Use ParsePHLAWD.R to collapse infraspecific taxa, and keep longest sequences for each species

######################### 4_Weight:  ####################################################################
######### Change Zorro weights to integers for RAxML...did not imporve phylogeny estimate, so excluded weights

######################### 6_Trees:  #####################################################################
### Use Congruifier (geiger) to ASSESS Trees 
### Modified from Jon Eastman

## Gene Trees: all *.nex files could be put into a folder, and then the code following for the concatenated could be looped over each

## Concatenated Tree:
genetree=read.nexus("output/6_Trees/Concatenated/RAxML_bipartitions.align.concat0530.1000.unweight.nex") # saved as rooted in FigTree

## USE CURATED FILE to resolve taxonomy 
##"linkage table", From Jon (NESCent working group on plant rates and traits)
tax=read.csv(file="output/7_Scaling/Congruify/fleshed_genera.csv", as.is=TRUE, row=1) 
tips=sapply(genetree$tip.label, function(x){
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
  
phy=genetree
swaptips=paste(1:length(tips),tips,sep="_")
phy$tip.label=swaptips
tax=SJ_tax
rownames(tax)=swaptips
res=congruify.phylo(fatol, phy, tax, tol=0, scale="PATHd8") # need to use PATHd8 to get res$phy
res
## DETERMINE WHICH NODES are congruified (and add asterisks to nodelabels)
out=res$phy
congruif=out$node.label%in%res$calibrations$MRCA
out$node.label=NULL
out=nodelabel.phylo(out, tax, strict=FALSE)
out$node.label=ifelse(congruif, paste(out$node.label, "*", sep=""), out$node.label)
### CHANGE tip labels back
out$tip.label=genetree$tip.label[match(swaptips, res$phy$tip.label)] 
### ADD Family to tip labels
out$tip.label=paste(genetree$tip.label[match(swaptips, res$phy$tip.label)], tax[,"family"], sep="=") 
  
pdf(file = "figs/trees/SJconcat0530_congruification.pdf", width=10, height=10) 
tree <- ladderize(out, right=F)
plot(tree, type="phylogram", cex=0.15, ) #x.lim=c(-10, 500))
nodelabels(out$node.label, frame="n", col="red", cex=0.25)
dev.off()

#see how the best node for lineage differs from the clade definition:
#missing from the clade in your tree OR unexpected within clade (but is found here in the tree)
out$FUN("Solanales")
out$FUN("Poaceae")

# * = Congruified Nodes
# Noded labels with ""  implies some inconsistency between the tips expected as defined in fleshed_genera.csv 
# and the subtended tips in the tree at the best matching node


######################### 7_Scaling:  ###################################################################
##### WRITES treePL for TreeScaling: Re-run congruifier, changing scale = "NA", and writing out .treePL files

setwd("/Users/hannahmarx/Documents/Idaho/Tank/Projects/SanJuans/FINAL/7_Scaling/")

SJ_tree=read.nexus("output/6_Trees/Concatenated/RAxML_bipartitions.align.concat0530.1000.unweight.nex") #saved rooted

tax=read.csv(file="output/7_Scaling/Congruify/fleshed_genera.csv", as.is=TRUE, row=1) 
tips=sapply(genetree$tip.label, function(x){
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
fatol=subset.phylo(atol, ftax, "family")

phy=SJ_tree
swaptips=paste(1:length(tips),tips,sep="_")
phy$tip.label=swaptips
tax=SJ_tax
rownames(tax)=swaptips
res1=congruify.phylo(fatol, phy, tax, tol=0, scale="NA") 

## WRITE OUT TREEPL FILE -- you'll need to do more than just run the exported file, but this gives you a start
# i.e. after prime, need to input scaling parameters 
nsites=5745 #SHOULD be number of nucleotides in alignment 
#write.treePL(res1$target, res1$calibrations, nsites, base="SJ0530.1000.unweight.treePL", opts=list(prime=TRUE))

### Read in treePL output to change tip labels back
treePL <- read.tree("SJ0530.1000.unweight.treePL.dated.tre")
treePL$tip.label=SJ_tree$tip.label[match(swaptips, res$phy$tip.label)] ### CHANGE tip labels back

#Write .tre file
########### USE THIS TREE FOR treePL!!!
#write.tree(treePL, file="SJ0530.1000.unweight.treePL.dated.rename.tre")


###### Check the treePL scaled ML tree using congruifier
## WRITE OUT FILE for matching nodes ######
res=congruify.phylo(fatol, phy, tax, tol=0, scale="PATHd8") ##SO IS PATHd8 THE ONLY OPTION FOR SCALE AT THIS POINT?
tt=cbind(original=SJ_tree$tip.label[match(swaptips, res$phy$tip.label)], faked=res$phy$tip.label)
#write.csv(tt, file="SJ_tips.csv", row.names=FALSE)

## DETERMINE WHICH NODES are congruified (and add asterisks to nodelabels)
out=res$phy
congruif=out$node.label%in%res$calibrations$MRCA
out$node.label=NULL
out=nodelabel.phylo(out, tax, strict=FALSE)
out$node.label=ifelse(congruif, paste(out$node.label, "*", sep=""), out$node.label)

out$tip.label=SJ_tree$tip.label[match(swaptips, res$phy$tip.label)] ### CHANGE tip labels back

## Congruified with PATHd8
#pdf("SJ_congruificationPATHd8.pdf") #width=10, height=10)
plot(out, type="phylogram", cex=0.25)# x.lim=c(-10, 500))
nodelabels(out$node.label, frame="n", col="red", cex=0.5)
axisPhylo(cex.axis=0.75)
#dev.off()
#write.tree(out, file="SJcongruifyPATHd8.tre") ## STORE tree with nodelabels

## Congruified with NA and dataed with treePL
#pdf("SJtreePL1000.pdf") #width=10, height=10)
plot(treePL, type="phylogram", cex=0.20)# x.lim=c(-10, 500))
nodelabels(out$node.label, frame="n", col="red", cex=0.5)
axisPhylo(cex.axis=0.75)
#dev.off()



