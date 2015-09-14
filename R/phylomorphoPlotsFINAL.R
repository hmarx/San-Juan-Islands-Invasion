##### Final Trait visulaization San Juan Isalnds Native : Invasives

setwd("~/Dropbox/Work/TankLab/Projects/SanJuans/Manuscript/Drafts/Figs.v1/")
#install.packages("~/Downloads/phytools_0.3-93.tar", type="source")
require(mvtnorm)
require(snowfall)
require(nlme)
require(qpcR)
require(reshape)
require(calibrate)
library(ggbiplot)
require(adephylo)
library(devtools)
library(ggbiplot)
require(phytools)
require(phylobase)
library(picante)

######################################## Read in Final Datasets #########
SJfinalTree <- read.tree("~/Dropbox/Work/TankLab/Projects/SanJuans/FINAL/7_Scaling/SJtreePL.bootstrap.tre")
SJfinalTree #367,  scaled with treePL, bootstrap support added to nodes

SJcomm <- read.csv("AppendixS1_SpList.csv", as.is=T, row.names=1) # community data (just pres/abs)
head(SJcomm)
class(SJcomm$Aleck_Rock)
dim(SJcomm) #415 80

SJtrait <- read.csv("AppendixS2_Traits.csv", as.is=T, row.names=1) # trait data (including "Status" column)
head(SJtrait)

dim(SJtrait) #415 7 
class(SJtrait)
#SJtrait <- apply(SJtrait, 2, list)
length(which(SJtrait["Status"] =="i")) #150 invasive species
length(which(SJtrait["Status"] =="n")) #265 native species

head(SJtrait)
statusCommunity(comm=SJcomm, col=4, lookup=SJtrait)
SJtrait[1, "Status"]
allislands <- c(1:ncol(SJcomm))
SJfoo <- lapply(allislands, function(x) statusCommunity(SJcomm, x, SJtrait))  #convert presnece (1) to n/i across community 
SJcommNew <- as.data.frame(do.call(cbind, SJfoo))
colnames(SJcommNew) <- colnames(SJcomm)
head(SJcommNew)
tail(SJcommNew)
SJcommNew$Aleck_Rock
class(SJcommNew$Aleck_Rock)
#SJtrait <- as.list(SJtrait)
com.list <- 1:ncol(SJcommNew) # make a list of the columns in the community data you want the fucntion to itterate across

##############################  Prune trait data to get just species with all trait data
head(SJtrait)
all.traits <- 1:length(SJtrait)
PCA_all <- lapply(all.traits, function(x) pruneTrait(phy=SJfinalTree, community=SJcommNew, traits=SJtrait, OneTrait=SJtrait[x], col=1))
head(PCA_all)
dim(PCA_all[[1]]) #366   3


############################### Get complete data for just Hight, SLA, seed mass
head(PCA_all) #output from pruneTrait; all traits for one community (here, all San Juans)
prunedData <- PCA_all #all traits
finalTree <- SJfinalTree

ThreetraitsComplete <- getCompleteDataset(prunedData=prunedData, finalTree=SJfinalTree,log=T, logCol=c(2,3), phylo4d=T)
head(ThreetraitsComplete) ### Log SLA, Leaf N
summary(ThreetraitsComplete) #81 taxa 
tipData <- tdata(ThreetraitsComplete, "tip")
treData <- extractTree(ThreetraitsComplete)
treData <- as(treData, "phylo")
plot(treData)
class(treData)

#ThreetraitsComplete2 <- getCompleteDataset(prunedData=prunedData, finalTree=SJfinalTree, log=F, phylo4d=F) 


############################ Phylogenetic PCA ############################ 
phylo.pca_SJ <-phyl.pca(tree=treData, Y=tipData, mode="cov")
summary(phylo.pca_SJ)
#Importance of components:
##  PC1       PC2        PC3        PC4        PC5
#Standard deviation     0.3220225 0.2667566 0.13001248 0.09009749 0.06083219
#Proportion of Variance 0.5093772 0.3495406 0.08303043 0.03987426 0.01817750
#Cumulative Proportion  0.5093772 0.8589178 0.94194825 0.98182250 1.00000000
pdf("5traits_biplotPCA.pdf")
biplot(phylo.pca_SJ)
dev.off()

phylo.pca_SJ2 <- ppca(ThreetraitsComplete,scale=FALSE,scannf=FALSE,nfposi=1,nfnega=1, method="Abouheif")
pdf("5traits_ppca.pdf")
plot(phylo.pca_SJ2, cex.label=0.05, tip.color=data_Complete.status, cex.label=0.1)
dev.off()

############################ Regular PCA ############################ 

# apply PCA - scale. = TRUE is highly advisable, but default is FALSE. 
pca_SJ <- prcomp(ThreetraitsComplete2[,2:6], center = TRUE, scale. = TRUE)
data_Complete.status <- ThreetraitsComplete2[,1]
names(data_Complete.status) <- rownames(ThreetraitsComplete2)

# plot method: plot varances and PC
plot(pca_SJ, type = "l")

# summary method
summary(pca_SJ)

#Importance of components:
#  PC1    PC2    PC3     PC4    PC5
#Standard deviation     1.4749 1.1544 0.8786 0.69218 0.4909 
#Proportion of Variance 0.4351 0.2665 0.1544 0.09582 0.0482
#Cumulative Proportion  0.4351 0.7016 0.8560 0.95180 1.0000
# First 2 PC account for 70% of the variance of the data


# Project the data on the first two PCs. Other PCs can be chosen through the argument choices of the function. 
#It colors each point according to the flowers??? species and draws a Normal contour line with 
#ellipse.prob probability (default to {68\%}) for each group.
pdf("5traits_regPCA.pdf")
g <- ggbiplot(pca_SJ, obs.scale = 1, var.scale = 1, 
              groups = data_Complete.status, ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
dev.off()

# plot each variables coefficients inside a unit circle to get insight on a possible interpretation for PCs 
pdf("5traits_regPCAunits.pdf")
theta <- seq(0,2*pi,length.out = 100)
circle <- data.frame(x = cos(theta), y = sin(theta))
p <- ggplot(circle,aes(x,y)) + geom_path()
loadings <- data.frame(pca_SJ$rotation, 
                       .names = row.names(pca_SJ$rotation))
p + geom_text(data=loadings, 
              mapping=aes(x = PC1, y = PC2, label = .names, colour = .names)) +
  coord_fixed(ratio=1) +
  labs(x = "PC1", y = "PC2")
dev.off()


############################ Test for Phylogenetic Signal in one trait ############################ 
names(SJtrait)
head(PCA_all) #list of traits including only species with data for each trait

PCA_all[[2]][[3]] <- sapply(PCA_all[[2]][[3]], log) 
PCA_all[[3]][[3]] <- sapply(PCA_all[[3]][[3]], log) 

nrow(PCA_all[[1]])  # Make a vector of trait meaurements
PCA_all[[1]]
statusCommunityREV <- function(comm, col){
  newCom <- data.frame()
  for (i in 1:nrow(comm)){
    if (comm[i, col]=="n"){
      y <- 0
    }
    if (comm[i, col]=="i"){
      y <- 1
    }
    newCom <- c(newCom, y)
  }
  names(newCom) <- rownames(comm)
  return(newCom)
}
n <- statusCommunityREV(comm=PCA_all[[1]], col=3)
PCA_all[[1]][[3]] <- n

commphyloSignal <- function(PCA_all, SJfinalTree){
  T <- list()
  for (i in 1:length(PCA_all)){
    prune <- treedata(phy=SJfinalTree, data=PCA_all[[i]])
    phy <- prune$phy
    trait.pr <- as.data.frame(prune$data)
    
    trait <- as.numeric(trait.pr[,3])  # Make a vector of trait meaurements
    names(trait) <- row.names(trait.pr) # Name the vector elements.
    Stat <- as.character(trait.pr[,2])
    names(Stat) <- row.names(trait.pr)
    #Map the continuous traits on the phylogeny 
    #dev.off()
    #contMap(phy,trait,res=30, fsize=0.4, lwd=2)
    
    #### Test for phylogenetic signal in traits
    trait.lambda <- phylosig(tree=phy, x=trait, method="lambda")
    
    trait.K <- phylosig(phy, trait, method="K", nsim=1000)
    
    tmp <- rbind(length(trait), trait.lambda$lambda, trait.lambda$logL, trait.K)
    colnames(tmp) = colnames(trait.pr[3]) 
    rownames(tmp) = c("n", "lambda", "lambda.logL", "K")
    T[[i]] <- tmp
    
  }
  return(T)
}
cps <- commphyloSignal(PCA_all=PCA_all, SJfinalTree=SJfinalTree)
cps.df <- do.call(cbind, cps)
cps.df

####NOT LOG SLA..cm.2.g.     leafN.... 
#Status  SLA..cm.2.g.     leafN.... LogSeedMass..mg. LogLeafletSize..cm.2. LogMaxHeight..m.
#n            366.00000000  2.000000e+02  110.00000000      322.0000000          190.00000000      244.0000000
#lambda         0.71320513  7.790356e-01    0.74855686        0.9888964            0.85698184        0.9718427
#lambda.logL -233.10083152 -1.183372e+03 -148.57808240     -579.6651033         -370.46689628     -317.3632170
#K              0.04022484  9.597868e-02    0.05223155        0.2208516            0.08215263        0.2128810

#### LOG SLA..cm.2.g.     leafN.... 
#Status SLA..cm.2.g.    leafN.... LogSeedMass..mg. LogLeafletSize..cm.2. LogMaxHeight..m.
#n            366.00000000  200.0000000 110.00000000      322.0000000          190.00000000      244.0000000
#lambda         0.71320513    0.8909395   0.82768145        0.9888964            0.85698184        0.9718427
#lambda.logL -233.10083152  -99.3322941 -43.92090351     -579.6651033         -370.46689628     -317.3632170
#K              0.04022484    0.1786478   0.07088418        0.2208516            0.08215263        0.2128810



