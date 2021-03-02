
Scenario = 1:4
#par(mfrow = c(2, 1))

#
library(rgeos)
library(magrittr)
library(ggplot2)
library(dplyr)
library(reshape2)

OutPutGraph <- NULL

Scenario = 1:2

for (Scen in Scenario) {
  if(Scenario %in% c(1:2)) load("X:/Luis/Model/GIT/WildBoar_Model/Inputs/habitats_Dep64.RData") else
  load("X:/Luis/Model/GIT/WildBoar_Model/Inputs/habitats_EA9.RData")

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
  
  Population      <- do.call(cbind, z$DayOutAniMat)
  Population      <- Population %>% data.frame
  Population$Time <- seq(nrow(Population))
  Population      <- melt(Population, id = "Time")
  
  Infected      <- do.call(cbind, z$DayOutInfMat)
  Infected      <- Infected %>% data.frame
  Infected$Time <- seq(nrow(Infected))
  Infected      <- melt(Infected, id = "Time")
  
  Carcasses      <- do.call(cbind, z$DayOutCarcMat)
  Carcasses      <- Carcasses %>% data.frame
  Carcasses$Time <- seq(nrow(Carcasses))
  Carcasses      <- melt(Carcasses, id = "Time")
  
  apply(res_newinf, 2, function(x) do.call(rbind, x)) -> z
  InfPerCells <- z$InfectedPerCells
  
  ### Calculate Epidemic Duration
 EpiDuration <- lapply(res, function(x){
    x= x$NewInfAnimals %>% data.frame
    c(unique(x$Iter), max(x$gTime) - 750)
  }) %>% do.call(rbind,.)

 PPopulation <- ggplot(Population, aes(x = Time, y = value, color = variable)) +
                        ggplot2::ggtitle(label = paste0("Scenario_", Scen), subtitle = "Population") + 
                        geom_line()
 
 PDensity <- "Do the same stuff but before mutate the density"
 
 PInfected <- ggplot(Infected, aes(x = Time, y = value, color = variable)) +
              ggtitle(label = paste0("Scenario_", Scen), subtitle = 'Infected WB') +
              geom_line()
 
 PCarcass <- ggplot(Carcasses, aes(x = Time, y = value, color = variable)) +
             ggtitle(label = paste0("Scenario_", Scen), subtitle = 'Infected Carcasses') +
             geom_line()
   
}
# plot(population[ ,1]/surface, typ = "l", main = paste0("Scenario_", Scen), xlab = "Days", ylab = "Wildboar Denisty/ km2", ylim = c(0, 5))
# apply(population[ ,2:99]/surface, 2, lines , col = rainbow(99))

# plot(c(0, 3000), c(0, 250), main = paste0("Scenario_", Scen), typ = 'n', xlab = "Days", ylab = "Incidence")
# NewInfAnimals <- res_newinf[ ,'NewInfAnimals']
# lapply(NewInfAnimals, function(x) lines(unique(x[ ,2]), table(x[ ,2]), col = rainbow(100)))

OutPutGraph[[Scen]] <- list(PPolutation, PDensity, PInfected, PCarcass, PIncidence)


# [[Scenario]][[Graph]]
OutPutGraph[[2]][[1]]

# Plot all Graphs per Scenario
plot_grid(OutPutGraph[[2]], tags = TRUE)

# ggplot(population2, aes(x=time, y=value, color= variable))+
#   title(main=paste0("Scenario_", Scen)) +
#   geom_line()

library(survival)
S<-Surv(EpiDuration[,2],event=rep(1,100))
plot(survfit(S~1))

ED1 <- as.data.frame(EpiDuration) %>% mutate(Scenario = 1) %>% rename(Iter = V1, EpDuration = V2)
ED2 <- as.data.frame(EpiDuration) %>% mutate(Scenario = 2) %>% rename(Iter = V1, EpDuration = V2)

EDT <- bind_rows(ED1, ED2)
survfit(Surv(EDT[ ,"EpDuration"], event = rep(1, 200)) ~Scenario, data = EDT) -> a
plot(a, col = 1:2)
legend('topright', legend = levels(factor(EDT$Scenario)), lty = 1, col = 1:2)

