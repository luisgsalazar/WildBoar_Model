####                                  Creating inputs for Departments 64 and 65

library(sp)
library(rgeos)
library(rgdal)
library(raster)
library(tmap)
library(tmaptools)

# Projection Lambert 93
CRS = "+proj=lcc +lat_1=45.898918964419 +lat_2=47.696014502038 +lat_0=46.8 +lon_0=2.337229104484 +x_0=600000 +y_0=2200000 +ellps=clrk80 +units=m +no_defs"

#Upload data from IGN
forest64 <- readOGR(dsn = "BDFORET/64/1_DONNEES/BDF_D064/FORMATION_VEGETALE.shp", layer = "FORMATION_VEGETALE")
#forest65 <- readOGR(dsn = "BDFORET/65/1_DONNEES/BDF_D065/FORMATION_VEGETALE.shp", layer = "FORMATION_VEGETALE")

#forests <- bind(forest64, forest65)
forest = spTransform(forest64, CRS)

tmap_mode('view')
tm_shape(forest64) + tm_fill("ESSENCE")

#Take only the department 64

cs = rep(3000, 2)
bb <- bbox(forest)
bb[ ,1] = bb[ ,1] - 10000
bb[ ,2] = bb[ ,2] + 10000
cc <- bb[ ,1] + (cs/2)  # cell offset
cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
grd <- GridTopology(cellcentre.offset = cc, cellsize = cs, cells.dim = cd)
sp_grd <- SpatialGridDataFrame(grd,
                               data = data.frame(id = 1:prod(cd)),
                               proj4string = CRS(proj4string(forest)))
#as(sp_grd, "SpatialPolygonsDataFrame") -> GRID

#Rasterize the zone
r <- raster(xmn = bb[1, 1], 
            ymn = bb[2, 1], 
            xmx = bb[1, 2], 
            ymx = bb[2, 2], 
            resolution = 3000, 
            crs = CRS)
r <- setValues(r, 1:ncell(r))
Gridr <- as(r,'SpatialPolygonsDataFrame')
names(Gridr@data) <- 'ID'

#Limit the zone to the dep 64
dep64 <- readOGR(dsn = "COMMUNE/Departement_64/geoflar-communes-2016.shp", layer = "geoflar-communes-2016")
dep64 <- spTransform(dep64, CRS)
z <- intersect(dep64, Gridr)
Grid2 <- Gridr[Gridr$ID %in% z$ID, ]
Grid2$ID = 1:length(Grid2)

#Calculate the forest cover from forest shapefile

library(doParallel)
library(plyr)
library(foreach)

rgeos::gSimplify(forest, tol = 0.00001) -> a
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


Forest_cover <- Forest_cover[order(Forest_cover[ ,1]), ]
Grid2$Forest_Cover = Forest_cover[ ,2]
Threshold = 10
Grid2$habitat = NA
Grid2$habitat[Grid2$Forest_Cover > Threshold]    <- 1
Grid2$habitat[Grid2$Forest_Cover < Threshold]    <- 2
Grid2$habitat[Grid2$Forest_Cover < Threshold/10] <- 3
colors = c("green","orange","red")
Grid2$color = colors[Grid2$habitat]
Depart64 <- gUnaryUnion(dep64)
habitats <- Grid2
rm(list = setdiff(ls(), c("habitats", 'Forest_cover', "z", "Depart64")))
save.image(file = "Department64Grid.RData")


# Compute neighbors of cells ----------------------------------------------

tmap_mode('view')
tm_shape(habitats) + tm_polygons('color')

library(doParallel)
library(plyr)
library(foreach)

nodes <- detectCores()
cl    <- makeCluster(nodes)
registerDoParallel(cl)
Neighbours <- foreach(i=1:length(habitats),.combine = 'rbind') %dopar% {
  library(rgeos)
  neighbour = rep(NA, 9) #8 neighbors + 1 to get the track on the ID (first element of row)
  indexes = which(gTouches(habitats[i, ], habitats, byid = T))
  neighbour[1:(length(indexes) + 1)] = c(i, indexes)
  neighbour
}
stopCluster(cl)
gc()

#Add the Neighbours data set to the grid
habitats@data[ ,5:12] = Neighbours[ ,2:9]

tmap_mode('view')
tm_shape(habitats) + tm_polygons("Habitat", fill)


#Get the coordinates of the polygons in the WB Grid
Habitatscoord <- coordinates(habitats)
habitats@data <- cbind(Habitatscoord, habitats@data)

colnames(habitats@data)
#colnames(habitats@data) <- c("ID", "Lon", "Lat", "Colors", "Forest_Cover", "Nb1", "Nb2", 'Nb3', 'Nb4', 'Nb5',
                             #'Nb6', 'Nb7', 'Nb8', "Habitat")

#habitats@data <- habitats@data[ ,c(3,1,2,6,4,7:14,5)] # get order as in the original script
#c("ID", "Lon", "Lat", "Forest_cover", 'Habitat', "Nb1", "Nb2", 'Nb3', 'Nb4', 'Nb5',
#               'Nb6', 'Nb7', 'Nb8',  'Colors')



write.csv(x = habitats@data, file = 'Input_dep64.csv', sep = ',', col.names = T, row.names = T )

#Remove everything except object on the list above
#rm(list=setdiff(ls(),c('habitats')))

#colnames(habitats@data) <- c("ID", "Lon", "Lat", 'Colors', "Forest_Cover", "Nb1", "Nb2", 'Nb3', 'Nb4', 'Nb5',
#                             'Nb6', 'Nb7', 'Nb8', "Habitat")

#habitats@data <- habitats@data[ ,c(1,2,3,5,14,6:13,4)]

save.image(file = "Dep64Grid.RData")
