
###### 6 Nov 2014  ##########
###### Hannah E. Marx #######

################## Data Cleaning ##################
## Functions that clean, prep, and summarize data

#### A function that converts your community data to incorporate status (native/invasive)
statusCommunity <- function(comm, col, lookup){
  newCom <- data.frame()
  for (i in 1:nrow(comm)){
    if (comm[i, col]==0){
      y <- 0
    }
    if (comm[i, col]==1){
      y <- lookup[rownames(comm[i,]), "Status"]
    }
    newCom <- c(newCom, y)
  }
  names(newCom) <- rownames(comm)
  return(newCom)
}

#### A function to summarize all trait data on one island 
traitSummary  <- function(phy, traits, col){
  # Prune data to get species with trait measure
  dataPruned <- (traits)[which(!traits[col] == "NA"),] ## prune data to trait
  dataPruned.mean <- mean(dataPruned[,col])
  dataPruned.median <- median(dataPruned[,col])
  dataPruned.min <- min(dataPruned[,col])
  dataPruned.max <- max(dataPruned[,col])
  
  n.dataPruned <- dataPruned[dataPruned$Status == "n",] # get all native
  i.dataPruned <- dataPruned[dataPruned$Status == "i",] # get all invasive
  
  # Numbers of species with trait 
  n.sp <- length(rownames(n.dataPruned))
  i.sp <- length(rownames(i.dataPruned))
  tot.sp <- length(rownames(dataPruned))
  
  n.sp.mean <- mean(n.dataPruned[,col])
  n.sp.median <- median(n.dataPruned[,col])
  n.sp.min <- min(n.dataPruned[,col])
  n.sp.max <- max(n.dataPruned[,col])
  
  i.sp.mean <- mean(i.dataPruned[,col])
  i.sp.median <- median(i.dataPruned[,col])
  i.sp.min <- min(i.dataPruned[,col])
  i.sp.max <- max(i.dataPruned[,col])
  
  # Which species in phylogeny have data for this trait
  tmp <- treedata(phy, dataPruned) 
  new.phy <- tmp$phy ## prune the tree 
  com.data <- as.data.frame(tmp$data)  ## subset the community data, to make sure in the right order
  tot.tips <- nrow(com.data)
  tot.tips.mean <- mean(as.numeric(as.character(com.data[,col])))
  tot.tips.median <- median(as.numeric(as.character(com.data[,col])))
  tot.tips.min <- min(as.numeric(as.character(com.data[,col])))
  tot.tips.max <- max(as.numeric(as.character(com.data[,col])))
  
  n.tips <- (com.data)[com.data$Status == "n", ] # get all the tips which are native
  i.tips <- (com.data)[com.data$Status == "i", ] # get all the tips which are invasive
  length.n.tips <- nrow(n.tips)
  length.i.tips <- nrow(i.tips)
  
  n.tips.mean <- mean(as.numeric(as.character(n.tips[,col])))
  n.tips.median <- median(as.numeric(as.character(n.tips[,col])))
  n.tips.min <- min(as.numeric(as.character(n.tips[,col])))
  n.tips.max <- max(as.numeric(as.character(n.tips[,col])))
  
  i.tips.mean <- mean(as.numeric(as.character(i.tips[,col])))
  i.tips.median <- median(as.numeric(as.character(i.tips[,col])))
  i.tips.min <- min(as.numeric(as.character(i.tips[,col])))
  i.tips.max <- max(as.numeric(as.character(i.tips[,col])))
  
  
  summ <- rbind(tot.sp, dataPruned.mean, dataPruned.median, dataPruned.min, dataPruned.max,
                n.sp, n.sp.mean, n.sp.median, n.sp.min, n.sp.max,
                i.sp, i.sp.mean, i.sp.median, i.sp.min, i.sp.max,
                tot.tips, tot.tips.mean, tot.tips.median, tot.tips.min, tot.tips.max,
                length.n.tips, n.tips.mean, n.tips.median, n.tips.min, n.tips.max,
                length.i.tips, i.tips.mean, i.tips.median, i.tips.min, i.tips.max)
  rownames(summ) <- c("Total species with trait", "Mean (tot)", "Median (tot)","Min (tot)","Max (tot)",
                      "Total natives with Trait", "Mean (nat)", "Median (nat)", "Min (nat)",  "Max (nat)",
                      "Total invasives with Trait", "Mean (inv)", "Median (inv)",  "Min (inv)",  "Max (inv)",  
                      "Total tips with trait", "Mean (tot tips)", "Median (tot tips)", "Min (tot tips)",  "Max (tot tips)",
                      "Total native tips with trait", "Mean (nat tips)", "Median (nat tips)", "Min (nat tips)",  "Max (nat tips)", 
                      "Totla invasive tips with trait", "Mean (inv tips)", "Median (inv tips)",  "Min (inv tips)",  "Max (inv tips)")
  
  colnames(summ) <- names(traits[col])
  return(summ)
  
}


