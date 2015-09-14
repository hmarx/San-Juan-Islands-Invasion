################################ Test for phylogenetic signal of traits in community ################################ 
SJfinalTree
head(SJtraitLog)

traits.inv <- SJtraitLog[SJtraitLog$Status == "i",]
traits.nat <- SJtraitLog[SJtraitLog$Status == "n",]

phylo.signal <- multiPhylosignal(SJtraitLog[c(2:ncol(SJtraitLog))], SJfinalTree)

# From picante_Intro:
# Phylogenetic signal is a quantitative measure of the degree to which phylogeny predicts the ecological similarity 
# of species. The K statistic is a measure of phylogenetic signal that compares the observed signal in a trait to 
# the signal under a Brownian motion model of trait evolution on a phylogeny (Blomberg et al. 2003). 
# K values of 1 correspond to a Brownian motion process, which implies some degree of phylogenetic signal or conservatism. 
# K values closer to zero correspond to a random or convergent pattern of evolution, 
# while K values greater than 1 indicate strong phylogenetic signal and conservatism of traits. 
# The statistical significance of phylo- genetic signal can be evaluated by comparing observed patterns of the 
# variance of independent contrasts of the trait to a null model of shuffling taxa labels across the tips of the phylogeny.

# The higher the K statistic, the more phylogenetic signal in a trait. 
# PIC.variance.P is the quantile of the observed phylogenetically independent contrast variance versus the null 
# distribution, which can be used as a 1-tailed P-value to test for greater phylogenetic signal than expected. 
# Traits with PIC.variance.P < 0.05 have non-random phylogenetic signal.

### From Swensen 2014:
## K = observed mean sq error / expected under browian motion model of evolution across given phylogeny


phylo.signal.inv <- multiPhylosignal(traits.inv[c(2:ncol(traits.inv))], SJfinalTree)
phylo.signal.nat <- multiPhylosignal(traits.nat[c(2:ncol(traits.nat))], SJfinalTree)

######### Pagel's lamda
### From Swenson 2014:
#search for the lambda value that transforms the original phylogeny such that the observed distribution
#of traits on the tips of the phylogeny is mirrored by that expected under
#Brownian Motion on the transformed phylogeny 
# A very low lambda increases the distance between sister taxa and therefore their expected trait difference under
#Brownian Motion. Thus, a lambda-based test of phylogenetic signal that results in
#low lambda values indicates very little phylogenetic signal in the trait data given the
#original phylogeny and a high lambda value indicates relatively more phylogenetic
#signal in the trait data given the original tree.

#Note that the test compares your reported lambda to the null hypothesis of a lambda equal to zero or
#a “star phylogeny.” It is therefore a true null hypothesis where relatedness explains
#none of the trait similarity between species.

seed.mass.lambda <- phylosig(x = SJtraitLog[SJfinalTree$tip.label, 2], tree = SJfinalTree, method="lambda", test=T)
maxHeight.mass.lambda <- phylosig(x = SJtraitLog[SJfinalTree$tip.label, 3], tree = SJfinalTree, method="lambda", test=T)
sla.mass.lambda <- phylosig(x = SJtraitLog[SJfinalTree$tip.label, 4], tree = SJfinalTree, method="lambda", test=T)
leafletSize.mass.lambda <- phylosig(x = SJtraitLog[SJfinalTree$tip.label, 5], tree = SJfinalTree, method="lambda", test=T)
leafN.mass.lambda <- phylosig(x = SJtraitLog[SJfinalTree$tip.label, 6], tree = SJfinalTree, method="lambda", test=T)


seed.complete <- treedata(phy = SJfinalTree, data = SJtraitLog[complete.cases(SJtraitLog$seedMass),])
height.complete <- treedata(phy = SJfinalTree, data = SJtraitLog[complete.cases(SJtraitLog$maxHeight),])
sla.complete <- treedata(phy = SJfinalTree, data = SJtraitLog[complete.cases(SJtraitLog$sla),])
leafsize.complete <- treedata(phy = SJfinalTree, data = SJtraitLog[complete.cases(SJtraitLog$leafletSize),])
leafN.complete <- treedata(phy = SJfinalTree, data = SJtraitLog[complete.cases(SJtraitLog$leafN),])

x.seed = phylo4d(seed.complete$phy, as.data.frame(seed.complete$data[,2]))
x.height = phylo4d(height.complete$phy, as.data.frame(height.complete$data[,3]))
x.sla = phylo4d(sla.complete$phy, as.data.frame(sla.complete$data[,4]))
x.leafsize = phylo4d(leafsize.complete$phy, as.data.frame(leafsize.complete$data[,5]))
x.leafN = phylo4d(leafN.complete$phy, as.data.frame(leafN.complete$data[,6]))

## Abouheif's tests for each trait
abo.test.seed <- abouheif.moran(x.seed)
abo.test.height <- abouheif.moran(x.height)
abo.test.sla <- abouheif.moran(x.sla)
abo.test.leafsize <- abouheif.moran(x.leafsize)
abo.test.leafN <- abouheif.moran(x.leafN)

plot(abo.test.seed)




