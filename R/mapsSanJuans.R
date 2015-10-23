
##### Map of San Juan Islands with Sampled Islands Highlighted ##### 
##### Following the online tutorial by Steven Brey, 14 April 2013: http://mazamascience.com/WorkingWithData/?p=1277
##### April 16, 2014 by Bethany Yollin http://mazamascience.com/WorkingWithData/?p=1494

source("R/addScaleBarMap.R")

localDir <- "output/11_Maps/"

############### Source of Orginal San Juan County Island Shoreline Boundary Layer  ########## 
## Data Source: https://sanjuanco.com/gis/gislib.aspx#BOUNDARY

# layerName is the name of the unzipped shapefile without file type extensions 
layerName <- "PSbasemap"  
## http://encdirect.noaa.gov
data_projected_orig <- readOGR(dsn="output/11_Maps/SourceShapefiles/PSbasemap/", layer=layerName) 

layerName <- "NOAA_Shorelines"  
# Read in the data https://sanjuanco.com/gis/gislib.aspx
data_projected_orig2 <- readOGR(dsn="output/11_Maps/SourceShapefiles/NOAA_Shorelines/", layer=layerName) 


############### Final Edited Shapefile ############### 
### Used layerName <- "PSbasemap"  
## http://encdirect.noaa.gov
## Edited in QGIS to remove unnecsary area: SanJuansMapComb.qgs
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


###### Transform shapefile to dataframe ###### 
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

# reproject the data onto a "longlat" projection
subsetTransform <- spTransform(data_projected, CRS("+proj=longlat"))

# determine the bounding box of the spatial object
b <- bbox(subsetTransform)

# get and plot a map
washingtonState <- ggmap(get_map(location = c(-123.11, 48.41291, -122.54390, 48.78919), zoom=10, color="bw", source="google", maptype = "satellite"))

subsetTransformFortified <- fortify(subsetTransform, region = "id")
subsetTransformFortified <- merge(subsetTransformFortified,
                                  subsetTransform@data, by.x = "id")

washingtonState + geom_polygon(data = subsetTransformFortified,
                                  aes(x = long, y = lat, group = group, color="red"),  fill="NA") 
                + theme(legend.position = "none", title = element_blank())



######################### Just sampled islands ######################### 
# Create a dataframe name (potentially different from layerName)
data_name <- "SANJUANS"

# Reproject the data onto a "longlat" projection and assign it to the new name
#assign(data_name,spTransform(data_projected_subset, CRS("+proj=longlat")))
assign(data_name,spTransform(data_projected, CRS("+proj=longlat")))

# The WRIA dataset is now projected in latitude longitude coordinates as a
# SpatialPolygonsDataFrame.  We save the converted data as .RData for faster
# loading in the future.
# save(list=c(data_name),file=paste(localDir,"SANJUANSproj.RData",sep="/"))

# Upon inspecting the metadata you can see that the first 19 islands in the San Juans
SANJUANS$Name[1:19]
WRIAPugetSound <- SANJUANS[!is.na(SANJUANS$Name),]
# Sanity check, plot WRIAPugetSound to make sure it looks like the subset we want
plot(WRIAPugetSound)

# Save Puget Sound data as an RData file which we can quickly "load()"
data_name <- "SANJUANS_SAMPLED"
save(list=c(data_name),file=paste(localDir,"SANJUANS_SAMPLED.RData",sep="/"))
SANJUANS$Shape_Leng
head(SANJUANS)

# Load metadata for San Juan Islands
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


length(unique(sanjuansDFsampled$Name)) #79


##################### PLOT ##################### 
usa <- map_data("state")
head(usa)

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
#pdf(file="output/11_Maps/wa.general.pdf")
wa.map
#dev.off()

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
#pdf(file="output/11_Maps/usa.general.pdf")
usa.map
#dev.off()

sj.map <- washingtonState + geom_polygon(data = sanjuansDFsampled,
                                                 aes(x = long, y = lat, group = group, color="red"), fill="NA") 
#P + scaleBar(lon = -124, lat = 45, distanceLon = 50, distanceLat = 10, distanceLegend = 20, dist.unit = "km", orientation = FALSE)
sj.map <- sj.map + scaleBar(lon = -123.255, lat = 48.35, distanceLon = 5, distanceLat = 1, distanceLegend = 2, dist.unit = "km", 
                            arrow.length = 5, arrow.distance = 3, arrow.North.size = 4)
sj.map <- sj.map + theme(legend.position = "none", title = element_blank()) 
#pdf(file="output/11_Maps/map.sampledIslands.scale.pdf")
sj.map
#dev.off()

sj.map2 <- ggplot() + geom_polygon(data = subsetTransformFortified, aes(x = long, y = lat, group = group), colour="gray", fill="NA") + coord_map()
sj.map2 <- sj.map2 +  geom_polygon(data = sanjuansDFsampled, aes(x = long, y = lat, group = group), fill="red") 
sj.map2 <- sj.map2 + theme(legend.position = "none", title = element_blank())

sj.map2 <- sj.map2 + scaleBar(lon = -123.255, lat = 48.35, distanceLon = 5, distanceLat = 1, distanceLegend = 2, dist.unit = "km", 
                            arrow.length = 5, arrow.distance = 3, arrow.North.size = 4)
sj.map2 <- sj.map2 + theme(legend.position = "none", title = element_blank(), panel.grid.minor = element_line(colour = NA), panel.grid.minor = element_line(colour = NA),
                           panel.background = element_rect(fill = NA, colour = "grey")) 
#pdf(file="output/11_Maps/map.sampledIslands.scale.3.pdf")
sj.map2
#dev.off()

