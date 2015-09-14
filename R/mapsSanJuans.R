################################ Map of San Juan Islands with Sampled Highlighted ################################ 
##### Following the online tutorial: http://mazamascience.com/WorkingWithData/?p=1277


library(maps)
library(mapdata)
library(sp)
library(rgdal)
library(rgeos)
library(ggmap)
library(mapplots)

source("R/addScaleBarMap.R")

########## Load San Juan County Island Shoreline Boundary Layer  ########## 
## Data Source: https://sanjuanco.com/gis/gislib.aspx#BOUNDARY
localDir <- "output/11_Maps/"

# layerName is the name of the unzipped shapefile without file type extensions 
layerName <- "NOAA_Shorelines"  
# Read in the data https://sanjuanco.com/gis/gislib.aspx
data_projected <- readOGR(dsn="output/11_Maps/NOAA_Shorelines/", layer=layerName) 

layerName <- "PSbasemap"  
## http://encdirect.noaa.gov
data_projected <- readOGR(dsn="output/11_Maps/SourceShapefiles/PSbasemap/", layer=layerName) 

layerName <- "Approach_Land_Area"  
## Layer of Washington State
### http://viewer.nationalmap.gov/viewer/
### https://www.sciencebase.gov/catalog/item/533c33b0e4b0f4f326e38352
data_projected <- readOGR(dsn="output/11_Maps/islands/", layer=layerName)

layerName <- "ne_10m_coastline"  
## http://www.naturalearthdata.com
data_projected <- readOGR(dsn="output/11_Maps/ExtraShapefiles/10m_physical/", layer=layerName) 

layerName <- "CUSPLine"  
## http://www.ngs.noaa.gov/NSDE/
data_projected <- readOGR(dsn="output/11_Maps/NSDE75317/", layer=layerName) 

layerName <- "composite_shoreline_final"  
## http://shoreline.noaa.gov/data/datasheets/composite.html
data_projected <- readOGR(dsn="output/11_Maps/noaa_composite/", layer=layerName) 


############### Fianl Edited Shapefile
### Used layerName <- "PSbasemap"  
## http://encdirect.noaa.gov
## Edited in QGIS to remove unnecsary area:
##### merged with shapefile from layerName <- "NOAA_Shorelines"  (not all islands are in this county subset...GISislnamesEdit.csv)
##### add attribute column "names" 
##### Labeled names using google earth as reference (SanJuanIslandsAddedManualQGIS.xlsx)

layerName <- "SanJuanIslandsMap"  
## http://shoreline.noaa.gov/data/datasheets/composite.html
data_projected <- readOGR(dsn="output/11_Maps/SanJuanIslandsMap/", layer=layerName) 

# What is this thing and what's in it?
class(data_projected)
slotNames(data_projected)
# It's an S4 "SpatialPolygonsDataFrame" object with the following slots:
# [1] "data"        "polygons"    "plotOrder"   "bbox"        "proj4string"

# What does the data look like with the default plotting command? 
plot(data_projected)

# Could use names(data_projected@data) or just:
names(data_projected)
#"FERRY"      "Name"       "SHAPE_Leng" "SHAPE_Area" "Acres"      "Sq_Miles"  
data_projected$Name

# Identify the attributes to keep and associate new names with them
#attributes <- c("FERRY", "Name", "SHAPE_Leng", "SHAPE_Area", "Acres", "Sq_Miles")

# user friendly names 
#newNames <- c("ferry", "name", "length", "area", "acres", "sqmiles")

# Subset the full dataset extracting only the desired attributes
#data_projected_subset <- data_projected[,attributes]

# Assign the new attribute names
#names(data_projected_subset) <- newNames


###### Transform shapefiel to dataframe
# add to data a new column termed "id" composed of the rownames of data
data_projected@data$id <- rownames(data_projected@data)

# create a data.frame from our spatial object
sanjuanPoints <- fortify(data_projected, region = "id")
head(sanjuanPoints)
# merge the "fortified" data with the data from our spatial object
sanjuansDF <- merge(sanjuanPoints, data_projected@data, by = "id")
head(sanjuansDF)
# NOTE : If we so choose, we could have loaded the plyr library to use the
#      : join() function. For those familiar with SQL, this may be a more
#      : intuitive way to understand the merging of two data.frames. An
#      : equivalent SQL statement might look something like this:
#      : SELECT *
#      : FROM dataProjected@data
#      : INNER JOIN watershedPoints
#      : ON dataProjected@data$id = watershedPoints$id

# library(plyr)
# watershedDF <- join(watershedPoints, dataProjected@data, by = "id")

p <- ggplot()
p <- p + geom_polygon(data = sanjuansDF, aes(x=long, y=lat, group=group)) 
p


# reproject the data onto a "longlat" projection
subsetTransform <- spTransform(data_projected, CRS("+proj=longlat"))

# determine the bounding box of the spatial object
b <- bbox(subsetTransform)

map("worldHires", "usa", xlim = c(-123.23778, -122.54390), ylim= c(48.41291, 48.78919))

# get and plot a map
washingtonState <- ggmap(get_map(location = c(-123.11, 48.41291, -122.54390, 48.78919), zoom=10, color="bw"))

washingtonState <- ggmap(get_map(location = c(-123.11, 48.41291, -122.54390, 48.78919), zoom=10, color="bw"))



a <- c(-126, 45, -116, 49)
washingtonStateOut <- ggmap(get_map(location = a,source="google", maptype = "hybrid", color="bw", zoom=6))


subsetTransformFortified <- fortify(subsetTransform, region = "id")
subsetTransformFortified <- merge(subsetTransformFortified,
                                  subsetTransform@data, by.x = "id")


washingtonState + geom_polygon(data = subsetTransformFortified,
                                  aes(x = long, y = lat, group = group, color="red"),  fill="NA") +
  theme(legend.position = "none", title = element_blank())



######################### Just sampled islands ######################### 
head(metadata)
namestmp <- gsub(rownames(metadata), pattern = "___", replacement = " - ")
namestmp <- gsub(namestmp, pattern = "__", replacement = "_#")
namestmp <- gsub(namestmp, pattern = "_", replacement = " ")

origsamp.tmp <- read.csv("output/11_Maps/SanJuansIslandOrigData.csv")
origsamp <- origsamp.tmp$Island

sampledMask <- which(SANJUANS$Name %in% origsamp)
sampledMaskName <- SANJUANS$Name[sampledMask]

# Create a SpatialPolygonsDataFrame subset
sampledIsls <- SANJUANS[sampledMask,]

# add to data a new column termed "id" composed of the rownames of data
sampledIsls@data$id <- rownames(sampledIsls@data)

# create a data.frame from our spatial object
sanjuanPointsSampled <- fortify(sampledIsls, region = "id")
head(sanjuanPointsSampled)
# merge the "fortified" data with the data from our spatial object
sanjuansDFsampled <- merge(sanjuanPointsSampled, sampledIsls@data, by = "id")
head(sanjuansDFsampled)

pdf(file="output/11_Maps/map.sampledIslands.pdf")
washingtonState + geom_polygon(data = sanjuansDFsampled,
                               aes(x = long, y = lat, group = group, color="red"), fill="NA") +
  theme(legend.position = "none", title = element_blank()) 
dev.off()

length(unique(sanjuansDFsampled$Name))


pdf(file="output/11_Maps/map.washingtonSampled.pdf")
washingtonStateOut + geom_polygon(data = sanjuansDFsampled,
                                  aes(x = long, y = lat, group = group, color="red"), fill="NA") +
  scale_x_continuous(limits = c(a[1],a[3])) +
  scale_y_continuous(limits = c(a[2],a[4])) +
  theme(legend.position = "none", title = element_blank())
dev.off()




############################


wa <- subset(usa, region %in% "washington")
wa.map <- ggplot() + geom_polygon(data = wa, aes(x = long, y = lat, group = group)) + coord_map()
#P + scaleBar(lon = -124, lat = 45, distanceLon = 50, distanceLat = 10, distanceLegend = 20, dist.unit = "km", orientation = FALSE)
wa.map <- wa.map + scaleBar(lon = -124, lat = 45, distanceLon = 50, distanceLat = 10, distanceLegend = 20, dist.unit = "km", 
             arrow.length = 50, arrow.distance = 30, arrow.North.size = 4)
wa.map <- wa.map + theme(panel.grid.minor = element_line(colour = NA), panel.grid.minor = element_line(colour = NA),
          panel.background = element_rect(fill = NA, colour = NA), axis.text.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(), axis.title = element_blank(),
          rect = element_blank(),
          plot.margin = unit(0 * c(-1.5, -1.5, -1.5, -1.5), "lines"))
pdf(file="output/11_Maps/wa.general.pdf")
wa.map
dev.off()

usa.map <- map_data("state")
usa.map <- ggplot() + geom_polygon(data = usa.map, aes(x = long, y = lat, group = group), fill="NA", color="grey") + coord_map()
usa.map <- usa.map + geom_polygon(data = wa, aes(x = long, y = lat, group = group), fill="grey") 
usa.map <- usa.map + scaleBar(lon = -130, lat = 26, distanceLon = 500, distanceLat = 100, distanceLegend = 200, dist.unit = "km")
usa.map <- usa.map + theme(panel.grid.minor = element_line(colour = NA), panel.grid.minor = element_line(colour = NA),
                           panel.background = element_rect(fill = NA, colour = NA), axis.text.x = element_blank(),
                           axis.text.y = element_blank(), axis.ticks.x = element_blank(),
                           axis.ticks.y = element_blank(), axis.title = element_blank(),
                           rect = element_blank(),
                           plot.margin = unit(0 * c(-1.5, -1.5, -1.5, -1.5), "lines"))
pdf(file="output/11_Maps/usa.general.pdf")
usa.map
dev.off()

sj.map <- washingtonState + geom_polygon(data = sanjuansDFsampled,
                                                 aes(x = long, y = lat, group = group, color="red"), fill="NA") 
#P + scaleBar(lon = -124, lat = 45, distanceLon = 50, distanceLat = 10, distanceLegend = 20, dist.unit = "km", orientation = FALSE)
sj.map <- sj.map + scaleBar(lon = -123.255, lat = 48.35, distanceLon = 5, distanceLat = 1, distanceLegend = 2, dist.unit = "km", 
                            arrow.length = 5, arrow.distance = 3, arrow.North.size = 4)
sj.map <- sj.map + theme(legend.position = "none", title = element_blank()) 
pdf(file="output/11_Maps/map.sampledIslands.scale.pdf")
sj.map
dev.off()



sj.map2 <- ggplot() + geom_polygon(data = subsetTransformFortified, aes(x = long, y = lat, group = group), colour="gray", fill="NA") + coord_map()
sj.map2 <- sj.map2 +  geom_polygon(data = sanjuansDFsampled, aes(x = long, y = lat, group = group), fill="black") 
sj.map2 <- sj.map2 + theme(legend.position = "none", title = element_blank())

sj.map2 <- sj.map2 + scaleBar(lon = -123.255, lat = 48.35, distanceLon = 5, distanceLat = 1, distanceLegend = 2, dist.unit = "km", 
                            arrow.length = 5, arrow.distance = 3, arrow.North.size = 4)
sj.map2 <- sj.map2 + theme(legend.position = "none", title = element_blank(), panel.grid.minor = element_line(colour = NA), panel.grid.minor = element_line(colour = NA),
                           panel.background = element_rect(fill = NA, colour = "grey")) 
pdf(file="output/11_Maps/map.sampledIslands.scale.2.pdf")
sj.map2
dev.off()








# Create a dataframe name (potentially different from layerName)
data_name <- "SANJUANS"

# Reproject the data onto a "longlat" projection and assign it to the new name
#assign(data_name,spTransform(data_projected_subset, CRS("+proj=longlat")))
assign(data_name,spTransform(data_projected, CRS("+proj=longlat")))


# The WRIA dataset is now projected in latitude longitude coordinates as a
# SpatialPolygonsDataFrame.  We save the converted data as .RData for faster
# loading in the future.
save(list=c(data_name),file=paste(localDir,"SANJUANSproj.RData",sep="/"))

# Upon inspecting the metadata you can see that the first 19 islands in the San Juans
SANJUANS$Name[1:19]

# For fun, save a subset including only only these 19 areas
PSWRIANumbers <- c(1:19)
WRIAPugetSound <- SANJUANS[!is.na(SANJUANS$Name),]
# Sanity check, plot WRIAPugetSound to make sure it looks like the subset we want
plot(WRIAPugetSound)

# Save Puget Sound data as an RData file which we can quickly "load()"
data_name <- "SANJUANS_SAMPLED"
save(list=c(data_name),file=paste(localDir,"SANJUANS_SAMPLED.RData",sep="/"))
SANJUANS$Shape_Leng
head(SANJUANS)

usa <- map_data("state")
head(usa)
wa <- subset(usa, region %in% "washington")



head(metadata)
namestmp <- gsub(rownames(metadata), pattern = "___", replacement = " - ")
namestmp <- gsub(namestmp, pattern = "__", replacement = "_#")
namestmp <- gsub(namestmp, pattern = "_", replacement = " ")

sampledMask <- which(SANJUANS$Name %in% namestmp)
sampledMaskName <- SANJUANS$Name[sampledMask]

# Create a SpatialPolygonsDataFrame subset
sampledIsls <- SANJUANS[sampledMask,]

# plot biggest area in blue with a title
plot(sampledIsls, col="blue")


p <- ggplot()
p <- p + geom_polygon(data=wa, aes(x=long, y=lat, group=group), colour="grey10", fill="lightgreen")
p <- p + geom_polygon(data=sampledIsls, colour="white", fill="grey10")
p <- p + geom_point(data = envAlps, aes(x=X_WGS84, y=Y_WGS84, color='corall'))
p <- p + geom_point(aes(x=2.3508, y=48.8567), color="red")
p <- p + geom_point(aes(x=5.7222, y=45.2002), color="red")
p <- p + annotate("text", x=5.7222, y=45.5, color="red", label="Grenoble")
p <- p + annotate("text", x=2.3508, y=49.5, color="red", label="Paris")
p <- p + theme(legend.position = "none", axis.title = element_blank())
p


ggmap(ecrins) + 
  geom_point(data = env.dispersion, aes(x=X_WGS84, y=Y_WGS84, color=mntd.obs.z, size=ntaxa.x)) +  
  scale_colour_gradientn(colours=rainbow(5)) +
  scale_size_continuous(range = c(3, 10)) +
  ggtitle("SES mntd")


#####

# Load the data that was converted to SpatialPolygonsDataFrame
# NOTE: This can be skipped (but does not have to be) if the spatial
# NOTE: objects are still in your workspace.

# Load and show the names of the attributes in WRIA
file <- paste(localDir,"WAWRIAs.RData",sep="/")
load(file) 
names(WRIA)

file <- paste(localDir,"WAWRIAs.RData",sep="/")
load(file)
names(WRIAPugetSound)

# Sweet.  We can see that this is the WRIA dataset we saved earlier

# NOTE: For more advanced users, slotNames(WRIA) will list the structures 
# in WRIA. Using the @ command allows you to grab a particular slot from the
# spatial object.  If you really want the HUGE AND GORY details of what's
# in this object, you can examine the full structure with str(WRIA).

# Here is how you would extract the contents of slots for easier use.
WriaData <- SANJUANS@data
WriaBbox <- SANJUANS@bbox

# We have a pretty good idea of what kind of data we are working with 
# and what it looks like. Now its time for the data to answer some
# questions and tell us a story.

# What is the biggest water resource area in Washington? 
maxArea <- max(SANJUANS$area)
# Create a 'mask' identifying the biggest area so we can find out its name
# NOTE:  Eaach 'mask' we create is a vector of TRUE/FALSE values that we
# NOTE:  will use to subset the dataframe.
biggestAreaMask <- which(SANJUANS$area == maxArea)
biggestAreaName <- SANJUANS$name[biggestAreaMask]
biggestAreaName

# Create a SpatialPolygonsDataFrame subset
biggestArea <- SANJUANS[biggestAreaMask,]

# plot biggest area in blue with a title
plot(biggestArea, col="blue")
title(biggestAreaName) 

# NOTE: Many more plot arguments can be explored by investigating 
# NOTE: the "SpatialPolygons" "plot-method" in the sp package

# I have heard of a water resource management area in Washington State
# called Pend Oreille.  Where is it located in this dataframe?
which(WriaData$name == "Swirl")

# Now we have isolated the watershed with the largest area as well as the
# fabled Pend Oreille.  Lets figure out how to highlight these regions when
# plotting all  regions. I have also heard that Lake Chelan is Beautiful.
# Lets isolate it as well.

# Each of the following makes a spatialPolygonsDataFrame subset, selecting 
# a specific region based on some selected attribute in WRIA.

WRIA_puget <- SANJUANS
WRIA_swirl <- SANJUANS[SANJUANS$name == "Swirl",]
WRIA_orcas <- SANJUANS[SANJUANS$name == "Orcas",]

# Check out what they look like plotted individually 
plot(WRIA_puget)
plot(WRIA_swirl)
plot(WRIA_orcas)

# For fun we will make 8 different watersheds 8 different colors!
watersheds <- c(1:8)
watershed_colors <- c("burlywood1","forestgreen","burlywood3","darkolivegreen3",
                      "cadetblue4","sienna3","cadetblue3","darkkhaki")

watershed_color_variation <- WRIA[WRIA$number %in% watersheds,]

# Plot some of the created spatial objects together
plot(WRIA_puget)
plot(WRIA_swirl,add=TRUE,col="yellow")
#plot(watershed_color_variation, add=TRUE, col=watershed_colors)
plot(WRIA_orcas,add=TRUE,col="red4")

# NOTE:  gCentroid is from the 'rgeos' package
library(rgeos)
# Find the center of each region and label lat and lon of centers
centroids <- gCentroid(SANJUANS, byid=TRUE)
centroidLons <- coordinates(centroids)[,1]
centroidLats <- coordinates(centroids)[,2]

# Find regions with center East of -121 longitude (East of Cascade mountains) 
eastSideMask <- centroidLons > -121

# Create spatialPolygonsDataFrame for these regions
WRIA_nonPacific <- SANJUANS[eastSideMask,]

# Find watersheds with area 75th percentile or greater
basinAreaThirdQuartile <- quantile(WRIA$area,0.75)
largestBasinsMask <- WRIA$area >= basinAreaThirdQuartile

WRIA_largest <- WRIA[largestBasinsMask,]

# To get legend and labels to fit on the figure we can change the size of the
# area plotting. bbox(WRIA) shows the bounding lat and long of WRIA.
bBox <- bbox(SANJUANS)
xlim <- c(bBox[1,1], bBox[1,2] )
ylim <- c(bBox[2,1], bBox[2,2] )

# Plot some of the different spatial objects and show lat-long axis, 
# label watersheds
pdf(file = "output/11_Maps/SanJuanLabels.pdf")
plot(SANJUANS,axes=TRUE,xlim=xlim,ylim=ylim)
#plot(WRIA_nonPacific,add=TRUE,col="red") 
#plot(WRIA_orcas,add=TRUE,col="navy")
#plot(WRIA_largest,add=TRUE,density=20,col="green1")
text(centroidLons, centroidLats, labels=c(1:length(SANJUANS$name)), col="red", cex=.1)
dev.off()

pdf(file = "output/11_Maps/SanJuanNames.pdf")
plot(SANJUANS,axes=TRUE,xlim=xlim,ylim=ylim)
#plot(WRIA_nonPacific,add=TRUE,col="red") 
#plot(WRIA_orcas,add=TRUE,col="navy")
#plot(WRIA_largest,add=TRUE,density=20,col="green1")
text(centroidLons, centroidLats, labels=SANJUANS$name, col="red", cex=.1)
dev.off()

SANJUANS$name[[336]]
write.csv(SANJUANS$name, file="output/11_Maps/GISislnames.csv")

title(main="San Juan Islands")
#labels <- c("Puget Sound Watersheds", "Washington's biggest watersheds",
            "drain to Pacific via Columbia river") 
#label_cols <- c("navy","green1","red")
#legend("bottomleft", NULL, labels, fill=label_cols, density, bty="n", cex=.8)



