rm(list=setdiff(ls(),c("ZoneBlanche","ZoneSurv")))
library(sp)
library(rgeos)
library(rgdal)
library(raster)
library(tmap)
library(tmaptools)
# Projection Lambert 93
CRS= "+proj=lcc +lat_1=45.898918964419 +lat_2=47.696014502038 +lat_0=46.8 +lon_0=2.337229104484 +x_0=600000 +y_0=2200000 +ellps=clrk80 +units=m +no_defs"




forets08<-readOGR(dsn="./Data/Shapefiles/BDFORET/08/1_DONNEES/BDF_D08/FORMATION_VEGETALE_08.SHP",
                  layer="FORMATION_VEGETALE_08")
forets02<-readOGR(dsn="./Data/Shapefiles/BDFORET/02/1_DONNEES/BDF_D02/FORMATION_VEGETALE_02.SHP",
                  layer="FORMATION_VEGETALE_02")
forets51<-readOGR(dsn="./Data/Shapefiles/BDFORET/51/1_DONNEES/BDF_D51/FORMATION_VEGETALE_51.SHP",
                  layer="FORMATION_VEGETALE_51")

forets54<-readOGR(dsn="./Data/Shapefiles/BDFORET/54/1_DONNEES/BDF_D54/FORMATION_VEGETALE_54.SHP",
                  layer="FORMATION_VEGETALE_54")
forets55<-readOGR(dsn="./Data/Shapefiles/BDFORET/55/1_DONNEES/BDF_D55/FORMATION_VEGETALE_55.SHP",
                  layer="FORMATION_VEGETALE_55")

forets57<-readOGR(dsn="./Data/Shapefiles/BDFORET/57/1_DONNEES/BDF_D57/FORMATION_VEGETALE_57.SHP",
                  layer="FORMATION_VEGETALE_57")

forets59<-readOGR(dsn="./Data/Shapefiles/BDFORET/59/1_DONNEES/BDF_D59/FORMATION_VEGETALE_59.SHP",
                  layer="FORMATION_VEGETALE_59")

forets62<-readOGR(dsn="./Data/Shapefiles/BDFORET/62/1_DONNEES/BDF_D62/FORMATION_VEGETALE_62.SHP",
                  layer="FORMATION_VEGETALE_62")

forets67<-readOGR(dsn="./Data/Shapefiles/BDFORET/67/1_DONNEES/BDF_D67/FORMATION_VEGETALE_67.SHP",
                  layer="FORMATION_VEGETALE_67")

forets88<-readOGR(dsn="./Data/Shapefiles/BDFORET/88/1_DONNEES/BDF_D88/FORMATION_VEGETALE_88.SHP",
                  layer="FORMATION_VEGETALE_88")

forets<-bind(forets02,forets08,forets51,
             forets54,forets55,forets57,
             forets59,forets62,forets67,
             forets88)

forets = spTransform(forets, CRS)
forets = aggregate(forets, "LIBELLE2")

rm(list = setdiff(ls(), c("ZoneBlanche", "ZoneSurv", "forets", "CRS")))
tmap_mode('view')

# Select the areas to study -----------------------------------------------
#Here we select the 10 neighboring departments of the WhiteZone

data <- readOGR(dsn = "./Data/Shapefiles/COMMUNE/COMMUNE.shp", layer = "COMMUNE")
colnames(data@data)
data <- data[which(data@data$INSEE_DEP %in% c("08","02","51","54","55","57","59","62","57","88")),]
raster::aggregate(data,"INSEE_DEP") -> a
a <- spTransform(a, CRS)
zoneLarge <- gUnaryUnion(a)


# Load surveillance, observation and white zones --------------------------
SurvZone <- readOGR(dsn = "./Data/Shapefiles/SurveillanceZone/SurveillanceZone.shp", layer = "SurveillanceZone")
SurvZone = spTransform(SurvZone, CRS)
WhiteZone<-readOGR(dsn=  "./Data/Shapefiles/SurveillanceZone/WhiteZone.shp", layer = "WhiteZone")
WhiteZone = spTransform(WhiteZone, CRS)
ObsZone <- readOGR(dsn = "./Data/Shapefiles/SurveillanceZone/ObservationZone.shp", layer = "ObservationZone")
ObsZone = spTransform(ObsZone, CRS)



# Load Roads data and extract the ones that are delimitating our area --------
Routes2 <- readOGR(dsn = "./Data/Shapefiles/ROUTE120/1_DONNEES_SHP/DONNEES/RESEAU_ROUTIER/TRONCON_ROUTE.SHP", layer = "TRONCON_ROUTE")
tm_shape(zoneLarge) + tm_borders()
Routes2[Routes2@data$NUM_ROUTE %in% c("A4","A34","A304","N51","N43","D986"), ] -> L
L = spTransform(L, CRS)

# Define the new working polygons - some effort to close it...  -----------
as(zoneLarge, 'SpatialLines') -> s
s = spTransform(s, CRS)
d <- function(x, y) sqrt((x[1]-y[1])^2 + (x[2]-y[2])^2)
dd = NULL
for (i in (1:length(s@lines[[1]]@Lines))){
  dd[i] = min(apply(s@lines[[1]]@Lines[[i]]@coords,1,d,L@lines[[1]]@Lines[[1]]@coords[1,]))
}
segment = which.min(dd)
val0 = min(apply(s@lines[[1]]@Lines[[segment]]@coords,1,d,L@lines[[1]]@Lines[[1]]@coords[1,]))
val1 = min(apply(s@lines[[1]]@Lines[[segment]]@coords,1,d,tail(L@lines[[1]]@Lines[[1]]@coords,1)))
index = which.min(apply(s@lines[[1]]@Lines[[segment]]@coords,1,d,L@lines[[1]]@Lines[[1]]@coords[1,]))
road.to.Boundary <- rbind(s@lines[[1]]@Lines[[segment]]@coords[index,],L@lines[[1]]@Lines[[1]]@coords[1,])
colnames(road.to.Boundary) <- c("Lon", "Lat")
road.to.Boundary <- data.frame(road.to.Boundary)
coordinates(road.to.Boundary) <- ~Lon+Lat
road.to.Boundary <- SpatialLines(list(Lines(list(Line(road.to.Boundary)), ID = "road.to.Boundary")))
proj4string(road.to.Boundary) <- proj4string(s)
plot(road.to.Boundary, col = 'green')
Roads <- gUnion(L, road.to.Boundary)

# Once the boundaries of the polygon are closed define the spatialP --------
lpi <- gIntersection(zoneLarge, Roads)               # intersect your line with 
blpi <- gBuffer(lpi, width = 10)  # create a very thin polygon 
dpi <- gDifference(zoneLarge, blpi)                # split using gDifference

# Plot the different polygons to find the one we're looking for
plot(SpatialPolygons(list(Polygons(list(dpi@polygons[[1]]@Polygons[[6]]), "4"))), add = TRUE, col = "lightgreen")
zoneLarge=SpatialPolygons(list(Polygons(list(dpi@polygons[[1]]@Polygons[[7]]), "4")))
# Our polygon is the seventh!

proj4string(zoneLarge) <- proj4string(s)
gBuffer(zoneLarge, width = 1000) -> tmp
Roads = gIntersection(tmp, Roads)


tm_shape(zoneLarge) + tm_borders() + tm_shape(Roads) + tm_lines("blue") + tm_shape(ObsZone) + 
  tm_borders("red") + tm_shape(WhiteZone) + tm_borders("red")


# Rasterize the zone ------------------------------------------------------

b = bbox(a)
r <- raster(xmn= b[1,1], ymn= b[2,1], xmx = b[1,2], ymx = b[2,2], resolution = 3000, crs = CRS)
r <- setValues(r, 1:ncell(r))
grid <- as(r,'SpatialPolygonsDataFrame')
names(grid@data) <- 'ID'

# Limit to the Zone including National Parks of Ardennes and North Vosges --------
intersect(zoneLarge,grid) -> z
Grid2 = grid[grid$ID %in% z$ID, ]
Grid2$ID = 1:length(Grid2)

library(doParallel)
library(plyr)
library(foreach)
rgeos::gSimplify(forets, tol = 0.00001) -> a

nodes <- detectCores()
cl <- makeCluster(round(3*nodes/4))
registerDoParallel(cl)
Forest_cover = foreach (i = 1:length(Grid2), .combine = "rbind")%dopar%{
  require(raster)
  # gSimplify(Gtmp,tol=0.00001)->Gtmp
  tmp = intersect(Grid2[i, ], a)
  if (!is.null(tmp)) {
    to_res = c(i, sum(area(tmp))/area(Grid2[i, ])*100)
  }else{
    to_res = c(i, 0)
  }
}
stopCluster(cl)
gc() 

Forest_cover <- Forest_cover[order(Forest_cover[,1]), ]
Grid2$Forest_Cov = Forest_cover[ ,2]
Threshold = 10
Grid2$habitat = NA
Grid2$habitat[Grid2$Forest_Cov>Threshold] <- 1
Grid2$habitat[Grid2$Forest_Cov<Threshold] <- 2
Grid2$habitat[Grid2$Forest_Cov<Threshold/10] <- 3
colors = c("green","orange","red")
Grid2$color = colors[Grid2$habitat]
observationZone <- gUnaryUnion(ObsZone)
habitats <- Grid2
rm(list = setdiff(ls(), c("Roads","WhiteZone","observationZone","habitats","zoneLarge")))
save.image(file = "Data_9kmGrid.RData")


