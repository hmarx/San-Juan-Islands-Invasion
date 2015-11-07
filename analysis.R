##### San Juans Final Datasets ######
########## 21 Oct. 2014 #############
########## Hannah E. Marx ###########

## loads and clean datasets

source("R/cleanData.R")

## The following packages will be required 
library(geiger)
library(ggplot2) 
library(gridExtra)
library(reshape2)
library(mvtnorm)
library(snowfall)
library(nlme)
library(qpcR)
library(reshape)
library(calibrate)
library(adephylo)
library(devtools)
library(phytools)
library(phylobase)
library(picante)
library(scales)
library(plyr)
library(maps)
library(maptools)
library(ggplot2)
library(grid)
library(maps)
library(mapdata)
library(sp)
library(rgdal)
library(rgeos)
library(ggmap)
library(mapplots)
library(cati)


################# PHYLOGENY
#### Final ML tree, scaled with treePL, bootstrap support added to nodes
SJfinalTree <- read.tree("data/SJtreePL.bootstrap.tre")
SJfinalTree #367
#plot(ladderize(SJfinalTree), type="fan", cex=.1)

################# FUNCTIONAL TRAITS
# from SJtraitsFINAL.csv == edited NESCent (edited to average values for synonymous taxa); 
# rows == species names (matching SJfinalTree$tip.label)
# columns == trait data (including "Status") 
# Status   seedMass  maxHeight sla   leafletSize    leafN  

SJtrait <- read.csv("data/SJtraits.csv", as.is=T, row.names=1) 
head(SJtrait)
dim(SJtrait) #415 8 
summary(SJtrait)
length(which(SJtrait["Status"] =="i")) #150 invasive species
length(which(SJtrait["Status"] =="n")) #265 native species
150/415 #36.14 %invasive


## plot the trait data to visualize
#pairs(SJtrait[,2:length(SJtrait)], main="Raw Trait Data")  # check distribution of data
#pairs(log(SJtrait[,2:length(SJtrait)]), main="Log Trait Data") 

SJtraitLog <- cbind(SJtrait[1],log10(SJtrait[,2:6])) #very skewed -> log
pairs(SJtraitLog[,2:length(SJtrait)], main="Log Trait Data") 
colnames(SJtraitLog) <- c("Status", "seedMass", "maxHeight", "sla", "leafletSize", "leafN")
summary(SJtraitLog)

################# COMMUNITY
# community data (just pres/abs)
# All San Juans + 80 uninhabited islands
SJcomm <- read.csv("data/SJspList.csv", as.is=T, row.names=1) 
head(SJcomm)
dim(SJcomm) #415 81

## convert presnece (1) to n/i across community 
allislands <- c(1:ncol(SJcomm)) # make a list of the columns in the community data you want the fucntion to itterate across
SJfoo <- lapply(allislands, function(x) statusCommunity(comm=SJcomm, col=x, lookup=SJtrait))  #convert presnece (1) to n/i across community 
SJcommNew <- as.data.frame(do.call(cbind, SJfoo))
colnames(SJcommNew) <- colnames(SJcomm)
dim(SJcommNew) #415  81
head(SJcommNew)

#### Remove islands with few specices that don't make sense for nearest native comparisons  
#East_Sucia_5_Island #no invasive species
#Little_Oak__2 #only one species
#Shag_Reef #no invasive species 
#Smaller_Island_near_Charles #no invasive species 
#Smallest_unnamed_island_by_Long_Island #only one native species
#Swirl_Rock_East #only one species
#Swirl_Rock_West #no invasive species 
#Unnamed_west_of_Castle_Island...only one invasive species

head(SJcommNew[, c("East_Sucia_5_Island", "Little_Oak__2", "Shag_Reef", "Smaller_Island_near_Charles", 
                   "Smallest_unnamed_island_by_Long_Island", "Swirl_Rock_East", "Swirl_Rock_West")])
remove.islands <- c("East_Sucia_5_Island", "Little_Oak__2", "Shag_Reef", "Smaller_Island_near_Charles", 
                    "Smallest_unnamed_island_by_Long_Island", "Swirl_Rock_East", "Swirl_Rock_West") 
SJcommNew  <- SJcommNew[, -which(names(SJcommNew) %in% remove.islands)]
length(SJcommNew) #74
head(SJcommNew) #74

### A lsit of just four islands for figures in manuscript
com.list.4 <- c(21, 16, 44, 1) 
four.islands <- c("East_Sucia_3_Island", "Deadman_Island", "Matia_Island", "All_SanJuanIslands")

################# META DATA
## Define size categories by quantile of island sizes
metadataFULL.rownames <- read.csv("data/SJcommSummary.csv", as.is=T, row.names=1) 
dim(metadataFULL.rownames) #81  12
metadataFULL.rownames <- metadataFULL.rownames[-which(rownames(metadataFULL.rownames) %in% remove.islands),]
head(metadataFULL.rownames)

metadataFULL <- read.csv("data/SJcommSummary.csv", as.is=T)
head(metadataFULL)
dim(metadataFULL) #81  13
metadataFULL <- metadataFULL[-which(metadataFULL$Row.names %in% remove.islands),]
dim(metadataFULL)

metadata <- read.csv("data/SJcommSummary.csv", as.is=T, row.names=1) 
remove.islands <- c("East_Sucia_5_Island", "Little_Oak__2", "Shag_Reef", "Smaller_Island_near_Charles", 
                    "Smallest_unnamed_island_by_Long_Island", "Swirl_Rock_East", "Swirl_Rock_West") 
metadata <- metadata[-which(rownames(metadata) %in% remove.islands),]
dim(metadata) # 74 12




