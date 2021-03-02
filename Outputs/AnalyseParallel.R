library(rgeos)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)
library(survival)


Scenario = 1:2
# par(mfrow = c(2, 1))
ScenarioTable <- NULL
OutPutGraph   <- NULL

for (Scen in Scenario) {
  if(Scenario %in% c(1:2)) load("X:/Luis/Model/GIT/WildBoar_Model/Inputs/habitats_Dep64.RData") else
  load("X:/Luis/Model/GIT/WildBoar_Model/Inputs/habitats_EA9.RData")
  
  #Study Area
  if(Scenario %in% c(1:2)) Study_Area = "Department_64" else 
                           Study_Area = "France-Belgium_Border"
  
  #Details
  if(Scenario %in% c(1, 3)) Details = "Without Hunting" else Details = "Hunting Season"

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
  
  ### Extract Results from main list
  apply(res_dayout, 2, identity) -> z
  
  Population <- do.call(cbind, z$DayOutAniMat)
  Population <- Population %>% data.frame %>% mutate(Time = seq(nrow(Population)))
  Population <- melt(Population, id = "Time")
  
  Density    <- Population %>% mutate(density = value/surface)
  
  Infected   <- do.call(cbind, z$DayOutInfMat)
  Infected   <- Infected %>% data.frame %>% mutate(Time = seq(nrow(Infected)))
  Infected   <- melt(Infected, id = "Time")
  
  Carcasses  <- do.call(cbind, z$DayOutCarcMat)
  Carcasses  <- Carcasses %>% data.frame %>% mutate(Time = seq(nrow(Carcasses)))
  Carcasses  <- melt(Carcasses, id = "Time")
  
  apply(res_newinf, 2, function(x) do.call(rbind, x)) -> z
  InfPerCells <- z$InfectedPerCells
  
  ### Calculate Epidemic Duration
 EpiDuration <- lapply(res, function(x){
    x= x$NewInfAnimals %>% data.frame
    c(unique(x$Iter), max(x$gTime) - 750)
  }) %>% do.call(rbind,.)

 PPopulation <- ggplot(Population, aes(x = Time, y = value, color = variable)) +
                        ggtitle(label = Study_Area, subtitle = paste("Population", Details)) + 
                        geom_line() + theme(legend.position = "none")
 
 PDensity  <- ggplot(Density, aes(x = Time, y = density, color = variable)) +
              ggtitle(label = Study_Area, subtitle = paste("WB Density", Details)) +
              geom_line() + theme(legend.position = "none")
 
 PInfected <- ggplot(Infected, aes(x = Time, y = value, color = variable)) +
              ggtitle(label = Study_Area, subtitle = paste("Infected WB", Details)) +
              geom_line() + theme(legend.position = "none")
 
 PCarcass <- ggplot(Carcasses, aes(x = Time, y = value, color = variable)) +
             ggtitle(label = Study_Area, subtitle = paste("Infected Carcass", Details)) +
             geom_line() + theme(legend.position = "none")
 
 #nameED        <- paste0("ED", Scen)
 ScenarioTable <- rbind(ScenarioTable, c(Scenario   = Scen,
                                         Study_Area = Study_Area,
                                         Details    = Details))
 
 OutPutGraph[[Scen]] <- list(PPopulation, PDensity, PInfected, PCarcass)
   
}

# OutPutGraph[[Scenario]][Graph]
# Plot all Graphs per Scenario

OutPutGraph[[1]][[1]]


# Survival Plot
# S<-Surv(EpiDuration[ ,2],event = rep(1, 100))
# plot(survfit(S ~ 1))

ED1 <- as.data.frame(EpiDuration) %>% mutate(Scenario = 1) %>% rename(Iter = V1, EpDuration = V2)
ED2 <- as.data.frame(EpiDuration) %>% mutate(Scenario = 2) %>% rename(Iter = V1, EpDuration = V2)
ED3 <- as.data.frame(EpiDuration) %>% mutate(Scenario = 3) %>% rename(Iter = V1, EpDuration = V2)
ED4 <- as.data.frame(EpiDuration) %>% mutate(Scenario = 4) %>% rename(Iter = V1, EpDuration = V2)

EDT <- bind_rows(ED1, ED2, ED3, ED4)

survfit(Surv(EDT[ ,"EpDuration"], event = rep(1, 400)) ~Scenario, data = EDT) -> a
plot(a, col = 1:4)
legend('topright', legend = levels(factor(EDT$Scenario)), lty = 1, col = 1:4)

