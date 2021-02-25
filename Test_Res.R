Scenario = 1:4
par(mfrow = c(2, 1))
library(rgeos)
library(dplyr)

for (Scen in Scenario){
  if(Scenario %in% c(1:2)) load("habitats_Dep64.RData") else
    load("habitats_EA9.RData")
  gArea(habitats)/1e6 -> surface
  
  load(paste0("./Outputs/Scenario_", Scen, "/", Scen, ".RData"))
  do.call(rbind, res) -> results
  
  res_dayout <- results[ ,c("DayOutInfMat", "DayOutPopMat", "DayOutAniMat", "DayOutCarcMat", "DayOutInfAniMat")]
  res_newinf <- results[ ,c("NewInfAnimals", "NewInfCarcass", "NewInfGroups", "InfectedPerCells")]
  
  apply(res_dayout, 2, identity) -> z
  population <- do.call(cbind, z$DayOutAniMat)
  Infected   <- do.call(cbind, z$DayOutInfMat)
  carcasses  <- do.call(cbind, z$DayOutCarcMat)
  
  apply(res_newinf, 2, function(x) do.call(rbind, x)) -> z
  InfPerCells <- z$InfectedPerCells
  
  plot(population[ ,1]/surface, typ = "l" ,xlab = "Days", ylab = "Wildboar Denisty/ km2", ylim = c(0, 5))
  apply(population[ ,2:99]/surface, 2, lines , col = rainbow(99))
  
  plot(x = c(1, nrow(population)), y = c(0, 25000), xlab = 'Time', ylab = 'Groups', type ='n')
  poplines <- for (i in 1:nrow(population)) {
    lines(1:nrow(population), population[ ,i], col = rainbow(i))
  
  }
  