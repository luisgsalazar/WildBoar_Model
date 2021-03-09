library(sp)
library(rgeos)
library(rgdal)
library(raster)
library(tmap)
library(tmaptools)
install.packages('dplyr')
library(dplyr)
Sys.setenv(LANGUAGE="en")
file.edit('~/Renviron')

# Projection Lambert 93
CRS= "+proj=lcc +lat_1=45.898918964419 +lat_2=47.696014502038 +lat_0=46.8 +lon_0=2.337229104484 +x_0=600000 +y_0=2200000 +ellps=clrk80 +units=m +no_defs"

#Load workspace with enlarged zone data

load("X:/Luis/Model/Working_space/Data_EnlargedZone.RData")
tmap_mode('view')
tm_shape(habitats) + tm_polygons('color')
#centroids <- coordinates(habitats) %>% as.data.frame()
#inputs <- habitats@data
#input_data <- bind_cols(... = inputs, centroids)
#colnames(input_data) <- c('ID', 'Habitat', 'Color', 'Forest_cover', 'Lon', 'Lat')

 library(doParallel)
 library(plyr)
 library(foreach)

 #Restart the numeration of rownames
 rownames(habitats@data) <- 1:nrow(habitats@data)
 
 #Calculate the Neighbors
 nodes <-detectCores()
 cl <- makeCluster(nodes)
 registerDoParallel(cl)
 Neighbours <- foreach(i=1:length(habitats),.combine = 'rbind') %dopar% {
   library(rgeos)
   neighbour = rep(NA,9) #8 neighbors + 1 to get the track on the ID (first element of row)
   indexes = which(gTouches(habitats[i,], habitats, byid = T))
   neighbour[1:(length(indexes)+1)] = c(i,indexes)
   neighbour
 }
 stopCluster(cl)
 gc()

 #Add the Neighbours data set to the grid
 habitats@data[,5:12] = Neighbours[,2:9]

 tmap_mode('view')
 tm_shape(habitats)+tm_polygons('habitat', fill )
 
 #Get the coordinates
 #Get the coordinates of the polygons in the WB Grid
Habitatscoord <- coordinates(habitats)
habitats@data <- cbind(Habitatscoord, habitats@data)

#Reorder columns
colnames(habitats@data)
colnames(habitats@data) <- c("Lon", "Lat", "ID", "Habitat", 'Colors', "Forest_Cover", "Nb1", "Nb2", 'Nb3', 'Nb4', 'Nb5',
                           'Nb6', 'Nb7', 'Nb8')

habitats@data <- habitats@data[ ,c(3,1,2,6,4,7:14,5)] # get order as in the original script
#c("ID", "Lon", "Lat", "Forest_cover", 'Habitat', "Nb1", "Nb2", 'Nb3', 'Nb4', 'Nb5',
#               'Nb6', 'Nb7', 'Nb8',  'Colors')

#Covert NAs of Neighbors into 0's, if not model will have errors
habitats@data[is.na(habitats@data)]=0

write.csv(x = habitats@data, file = 'Input_enlarged.csv', sep = ',', col.names = T, row.names = T )

#Remove everything except object on the list above
rm(list=setdiff(ls(),c('habitats')))
