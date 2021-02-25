
Scenario = 1:4
#par(mfrow = c(2, 1))

#
library(rgeos)
library(magrittr)
library(ggplot2)

OutPutGraph <- NULL

for (Scen in Scenario) {
  if(Scenario %in% c(1:2)) load("habitats_Dep64.RData") else
  load("habitats_EA9.RData")
  
  #Surface in km2
  gArea(habitats)/1e6 -> surface
  
  #Load Scenario List Results
  load(paste0("./Outputs/Scenario_", Scen, "/", Scen, ".RData"))
  do.call(rbind, res) -> results
  
  ### Divide into type of results
  res_dayout <- results[ ,c("DayOutInfMat", "DayOutPopMat", "DayOutAniMat", "DayOutCarcMat", "DayOutInfAniMat")]
  res_newinf <- results[ ,c("NewInfAnimals", "NewInfCarcass", "NewInfGroups", "InfectedPerCells")]
  
  #To check if Iters are messed up
  # which(sapply(res, function(x){
  #   x <- x$NewInfAnimals %>% data.frame
  #   length(unique(x$Iter)) != 1
  # }))
  
  ### Divide Results
  apply(res_dayout, 2, identity) -> z
  Population <- do.call(cbind, z$DayOutAniMat)
  Infected   <- do.call(cbind, z$DayOutInfMat)
  carcasses  <- do.call(cbind, z$DayOutCarcMat)
  
  apply(res_newinf, 2, function(x) do.call(rbind, x)) -> z
  InfPerCells <- z$InfectedPerCells
  
  ### Calculate Epidemic Duration
  NewInfAnimals <- lapply(res, function(x){
    x= x$NewInfAnimals %>% data.frame
    c(unique(x$Iter), max(x$gTime) - 750)
  }) %>% do.call(rbind,.)

 PDensity <- plot(population[ ,1]/surface, typ = "l", main = paste0("Scenario_", Scen), xlab = "Days", ylab = "Wildboar Denisty/ km2", ylim = c(0, 5))
 apply(population[ ,2:99]/surface, 2, lines , col = rainbow(99))
  
   
}


# plot(population[ ,1]/surface, typ = "l", main = paste0("Scenario_", Scen), xlab = "Days", ylab = "Wildboar Denisty/ km2", ylim = c(0, 5))
# apply(population[ ,2:99]/surface, 2, lines , col = rainbow(99))


p2 <- plot(x = c(1, nrow(population)), y = c(0, 25000), main = paste0("Scenario_", Scen), xlab = 'Time', ylab = 'Wildboars', type ='n')
poplines <- for (i in 1:nrow(population)) {
  lines(1:nrow(population), population[ ,i], col = rainbow(i))}

p1 <- plot(Infected[ ,1], typ = "l", main = paste0("Scenario_", Scen), xlab = "Days", ylab = "Infected Wild boars",
     ylim = c(0, max(Infected) + 100))
apply(Infected[ ,2:99], 2, lines, col = rainbow(99))

# plot(c(0, 3000), c(0, 250), main = paste0("Scenario_", Scen), typ = 'n', xlab = "Days", ylab = "Incidence")
# NewInfAnimals <- res_newinf[ ,'NewInfAnimals']
# lapply(NewInfAnimals, function(x) lines(unique(x[ ,2]), table(x[ ,2]), col = rainbow(100)))

# plot(carcasses[ ,1], typ = "l", main = paste0("Scenario_", Scen), xlab = "Days", ylab = "Infected carcasses", 
#      ylim = c(0, max(carcasses) + 100))
# apply(carcasses[ ,2:99], 2, lines, col = rainbow(99))

OutPutList[[Scen]] <- list(p1, p2)
}

# voir graph 1 scenario 1
OutPutList[[1]][[1]]

plot_grid(OutPutList[[1]], tags = TRUE)

# tmp <- NULL
# 
# for(i in 1:4){
#   p1 <- qplot(c(1,4), c(2,3))
#   p2 <- qplot(c(1,2), c(2,3))
#   tmp[[i]] <- list(p1,p2)
# }
# 
# tmp[[1]][[2]]
# 
# 
# population <- population %>% data.frame
# population$time <- seq(nrow(population))
# population2 <- melt(population,id = "time")
# 
# 
# ggplot(population2, aes(x=time, y=value, color= variable))+
#   title(main=paste0("Scenario_", Scen)) +
#   geom_line()
