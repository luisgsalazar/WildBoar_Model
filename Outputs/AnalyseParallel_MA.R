library(rgeos)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)
library(survival)
library(gridExtra)
library(grid)


# List of Output Graphs ---------------------------------------------------
OutPutGraph = NULL

#par(mfrow = c(2, 1))
Scenarios <- data.frame(Sc               = c(1:8),
                        Zone             = rep(c(rep("DPA", 2), rep("FBB", 2)), 2),
                        ProbHarvest      = rep(c(0, 0.40), 4),
                        ProbHarvestAM    = rep(c(0, 0.60), 4),
                        ProbMovHunt      = rep(c(0, 0.30), 4),
                        GroupMoveHunt    = rep(c("No", "Yes"), 4),
                        TimeSeedInf      = c(rep(3650, 4), rep(730, 4)),
                        MaxIterations    = c(rep(10, 4),rep(100, 4)),
                        MaxDays          = c(rep(5*365, 4), rep(8*365, 4)),
                        MaxBCap          = rep(8, 8),
                        SurvivalProbAdF  = rep(0.85, 8),                   # Yearly Survival Prob for adults Kueling et al. (2013)    
                        SurvivalProbAdM  = rep(0.75, 8),                   # Yearly Survival Prob for adults Kueling et al. (2013) 
                        SurvivalProbSAdF = rep(0.75, 8),                   # Yearly Survival Prob for sub-adults (we use adults values as they are closer to Lange et al. 2015) Kueling et al. (2013)
                        SurvivalProbSAdM = rep(0.70, 8)
                        
)


# Data extract ------------------------------------------------------------


for (Scen in Scenarios$Sc) {
  index_Dep64 = c(1:2, 5:6, 9:10)
  if(Scen %in% index_Dep64) load("Inputs/habitats_Dep64.RData") else
    load("Inputs/habitats_EA9.RData")
  
  #Study Area
  if(Scen %in% index_Dep64) Study_Area = "Pyrénées-Atlantiques Department" else 
                            Study_Area = "France-Belgium Border"
  
  #Details
  if(Scen %in% seq(1, 7, 2)) Details = "No Hunting" else Details = "Hunting Season"
  
  #Season Introduction
  if(Scen %in% c(1:4))  Season = "None" 
  if(Scen %in% c(5:8))  Season = "Winter" 
  if(Scen %in% c(9:12)) Season = "Summer"

  #Surface in km2
  gArea(habitats)/1e6 -> surface
  
  #Load Scenario List Results
  load(paste0("Outputs/OutputApply/scenario ", Scen, ".RData"))
  do.call(rbind, res) -> results
  
  ### Divide type of results
  res_dayout <- results[ ,c("DayOutInfMat", "DayOutPopMat", "DayOutAniMat", "DayOutCarcMat", "DayOutInfAniMat")]
  res_newinf <- results[ ,c("NewInfAnimals", "NewInfCarcass", "NewInfGroups", "InfectedPerCells")]
  
  
  ### Extract Results from main list
  apply(res_dayout, 2, identity) -> z
  
  Population  <- do.call(cbind, z$DayOutAniMat)
  Popquant    <- Population %>% data.frame %>% 
                 apply(., 1, quantile, probs = c(0.5, 0.025, 0.975)) %>% t(.)
  Density     <- Population %>% data.frame %>% 
                 apply(., 1, divide_by, surface) %>% t(.)
  Denquant    <- Density %>% data.frame() %>%
                 apply(., 1, quantile, probs = c(0.5, 0.025, 0.975)) %>% t(.)
 
  
  Population  <- Population %>% data.frame %>% mutate(Time = seq(nrow(Population)))
  Population  <- melt(Population, id = "Time")
  Density <-     Population %>% mutate(density = value/surface)
  
  Cells      <- do.call(cbind, z$DayOutPopMat)
  Cellsquant <- Cells %>% data.frame %>%
    apply(., 1, quantile, probs = c(0.5, 0.025, 0.975)) %>% t(.)
  
  Cells <- Cells %>% data.frame %>% mutate(Time = seq(nrow(Cells)))
  Cells <- melt(Cells, id = "Time")
  
  Infected   <- do.call(cbind, z$DayOutInfAniMat)
  Infquant   <- Infected %>% data.frame %>%
                 apply(., 1, quantile, probs = c(0.5, 0.025, 0.975)) %>% t(.)
  InfCum     <- Infquant %>% data.frame %>% mutate(cum_median = cumsum(Infquant[ ,1]))
  Infected   <- Infected %>% data.frame %>% mutate(Time = seq(nrow(Infected)))
  Infected   <- melt(Infected, id = "Time")
  
  Carcasses  <- do.call(cbind, z$DayOutCarcMat)
  Carquant   <- Carcasses %>% data.frame %>%
                 apply(., 1, quantile, probs = c(0.5, 0.025, 0.975)) %>% t(.)
  CarcCum    <- Carquant %>% data.frame %>% mutate(cum_median = cumsum(Carquant[ ,1]))
  Carcasses  <- Carcasses %>% data.frame %>% mutate(Time = seq(nrow(Carcasses)))
  Carcasses  <- melt(Carcasses, id = "Time")
  
 if (Scen > 4){
  apply(res_newinf, 2, function(x) do.call(rbind, x)) -> z
  InfPerCells <- z$InfectedPerCells
  
  InfAnimals <- z$NewInfAnimals
  InfAnimals <- InfAnimals %>% data.frame %>% group_by(Iter, gTime) %>% summarise(n())
  by(InfAnimals$`n()`, InfAnimals$Iter, cumsum) -> tmp
  InfAnimals$Infected_Cum <- do.call(what = 'c', tmp)
  CumInfAnimals <- InfAnimals %>% group_by(., Iter) %>% summarise(median = median(Infected_Cum)) %>% summarise(quant = quantile(median, probs = c(0.5, 0.025, 0.975)))
  
  InfGroups <- z$NewInfGroups
  InfGroups <- InfGroups %>% data.frame %>% group_by(Iter, gTime) %>% summarise(n())
    by(InfGroups$`n()`, InfGroups$Iter, cumsum) -> tmp
    InfGroups$Infected_Cum <- do.call(what = 'c', tmp)
  CumInfGroups <- InfGroups %>% group_by(., Iter) %>% summarise(median = median(Infected_Cum)) %>% summarise(quant = quantile(median, probs = c(0.5, 0.025, 0.975)))
 
  InfCarcass <- z$NewInfCarcass
  InfCarcass <- InfCarcass %>% data.frame %>% group_by(Iter, gTime) %>% summarise(n())
  by(InfCarcass$`n()`, InfCarcass$Iter, cumsum) -> tmp
  InfCarcass$Infected_Cum <- do.call(what = 'c', tmp)
  CumInfCarcass <- InfCarcass %>% group_by(., Iter) %>% summarise(median = median(Infected_Cum)) %>% summarise(quant = quantile(median, probs = c(0.5, 0.025, 0.975))) 
  
  }

  ### Calculate Epidemic Duration
 # EpiDuration <- lapply(res, function(x){
 #    x= x$NewInfAnimals %>% data.frame
 #    c(unique(x$Iter), max(x$gTime) - Scenarios[Scenarios$Sc==Scen,"TimeSeedInf"])
 #    }) %>% do.call(rbind,.)


# Epidemic Duration -------------------------------------------------------

  
 Sim = NULL
 
for (j in 1:Scenarios$MaxIterations[Scenarios$Sc == Scen]) 
  Sim = c(Sim, rep(j, Scenarios$MaxDays[Scenarios$Sc == Scen]))
Carcasses$Sim = Sim
Infected$Sim  = Sim
EpiDuration   = apply(rbind(
  tapply(Carcasses$value, Carcasses$Sim, function(x) max(which(x > 0))),
  tapply(Infected$value,  Infected$Sim,  function(x) max(which(x > 0)))), 2, max) - Scenarios[Scenarios$Sc == Scen, "TimeSeedInf"]
EpiDuration = tibble(Iter = names(EpiDuration), EpDuration = as.numeric(EpiDuration))

EpiDuration <- summarise(EpiDuration, quant = quantile(EpDuration, probs = c(0.5, 0.025, 0.975)))


# Plots -------------------------------------------------------------------


 ### All Iterations Graphic 
 PPopulationA <- ggplot(Population, aes(x = Time, y = value, color = variable)) +
                ggtitle(label = Study_Area) +
                geom_line() + theme(legend.position = "none")
    
 ### Median and quantile plots        
 PPopulationQ <- ggplot(Population, aes(x = Time, y = value)) +
                ggtitle(label = Study_Area) +
                ylab("Individuals") + xlab("Time (days)") +
                geom_line(colour = "cornflowerblue", alpha = 0.3) + theme(legend.position = "none") +
                stat_summary(fun.min = function(z) { quantile(z,0.025) },
                             fun.max = function(z) { quantile(z,0.975) },
                             fun = median, geom = "smooth", colour = "darksalmon", size = 1)
 
 
 PCellsQ <- ggplot(Cells, aes(x = Time, y = value)) +
            ggtitle(label = Study_Area) +
            ylab("Individuals") + xlab("Time (days)")+
            geom_line(colour = "cornflowerblue", alpha = 0.3) + 
            theme(legend.position = "none") +
            stat_summary(fun.min = function(z) { quantile(z,0.25) },
            fun.max = function(z) { quantile(z,0.75) },
            fun = median, geom = "smooth", colour = "darksalmon", size = 1)
  
  PDensityA <- ggplot(Density, aes(x = Time, y = density, color = variable)) +
              ggtitle(label = Study_Area, subtitle = paste("WB Density", "/", Details, "/", Season)) +
              geom_line() + theme(legend.position = "none")

  PDensityQ <- ggplot(Density, aes(x = Time, y = density)) +
                 ggtitle(label = Study_Area, subtitle = paste("WB Density", "/", Details, "/", Season)) +
                 geom_line(colour = "cornflowerblue", alpha = 0.3) + theme(legend.position = "none") +
                 stat_summary(fun.min = function(z) { quantile(z,0.25) },
                              fun.max = function(z) { quantile(z,0.75) },
                              fun = median, geom = "smooth", colour = "coral", size = 1)

  PInfectedA <- ggplot(Infected, aes(x = Time, y = value, color = variable)) +
                ggtitle(label = Study_Area, subtitle = paste("Infected\n WB", "/", Details, "/", Season)) +
                geom_line() + theme(legend.position = "none")
  
  PInfectedQ <- ggplot(Infected, aes(x = Time, y = value)) +
                ggtitle(label = Study_Area) +
                ylab("Infected\n wild boars") + xlab("Time (days)")+
                geom_line(colour = "cornflowerblue", alpha = 0.3) + theme(legend.position = "none") +
                stat_summary(fun.min = function(z) { quantile(z,0.25) },
                  fun.max = function(z) { quantile(z,0.75) },
                  fun = median, geom = "smooth", colour = "brown2", size = 1)
 
  PCarcassA <- ggplot(Carcasses, aes(x = Time, y = value, color = variable)) +
               ggtitle(label = Study_Area, subtitle = paste("Infected\n Carcasses", "/", Details, "/", Season)) +
               geom_line() + theme(legend.position = "none")
  
  PCarcassQ <- ggplot(Carcasses, aes(x = Time, y = value)) +
               ggtitle(label = Study_Area) +
               ylab("Infected\n Carcasses")+xlab("Time (days)")+
               geom_line(colour = "cornflowerblue", alpha = 0.3) + theme(legend.position = "none") +
               stat_summary(fun.min = function(z) { quantile(z,0.25) },
                            fun.max = function(z) { quantile(z,0.75) },
                            fun = median, geom = "smooth", colour = "darkgoldenrod2", size = 1)
 
 # ScenarioTable <- rbind(ScenarioTable, c(Scenario   = Scen,
 #                                         Study_Area = Study_Area,
 #                                         Details    = Details,
 #                                         Season     = Season))
 if (Scen > 4) {
  
 OutPutGraph[[Scen]] <- list(PPopulationA, PPopulationQ, Popquant,
                             PDensityA,    PDensityQ,    Denquant,
                             PCellsQ,      Cellsquant,                 CumInfGroups,
                             PInfectedA,   PInfectedQ,   Infquant,     CumInfAnimals,
                             PCarcassA,    PCarcassQ,    Carquant,     CumInfCarcass,
                             EpiDuration)  } else
                    
  OutPutGraph[[Scen]] <- list(PPopulationA, PPopulationQ, Popquant,
                              PDensityA,    PDensityQ,    Denquant,
                              PCellsQ,      Cellsquant,   
                              PInfectedA,   PInfectedQ,   Infquant,    
                              PCarcassA,    PCarcassQ,    Carquant,
                              EpiDuration)
   
}

OutPutGraph[[5]][[11]]
View(OutPutGraph[[5]][[12]])
# OutPutGraph[[Scenario]][Graph/Table]
# Plot all Graphs per Scenario

# Population Dynamics Graphs

# grid.arrange(
#   grobs = list(OutPutGraph[[1]][[2]], OutPutGraph[[3]][[2]],
#                OutPutGraph[[2]][[2]], OutPutGraph[[4]][[2]]),
#   widths = c(2, 2),
#   layout_matrix = rbind(c(1, 2), c(3, 4)),
#   top = textGrob("Population Dynamics", gp = gpar(fontsize = 20, font = 3))
# )

pop <- ggarrange(OutPutGraph[[3]][[2]],  OutPutGraph[[1]][[2]],
                 OutPutGraph[[4]][[2]],  OutPutGraph[[2]][[2]]
                 + rremove("x.text"),
                 labels = c("A", "B", "C", "D"), font.label = 16,
                 ncol = 2, nrow = 2)
annotate_figure(pop,
                top = text_grob("Population Dynamics", face = "bold", size = 16))


no_hunt <- ggarrange(OutPutGraph[[7]][[2]], OutPutGraph[[5]][[2]],
                     OutPutGraph[[7]][[11]], OutPutGraph[[5]][[11]],
                     OutPutGraph[[7]][[15]], OutPutGraph[[5]][[15]]
                     + rremove("x.text"),
                     labels = c("A", "B", "C", "D", "E", "F"), font.label = 16,
                     ncol = 2, nrow = 3)
annotate_figure(no_hunt, top = text_grob("ASF Dynamics - No Hunting", face = "bold", size = 16))

hunt <- ggarrange(OutPutGraph[[8]][[2]], OutPutGraph[[6]][[2]],
                  OutPutGraph[[8]][[11]], OutPutGraph[[6]][[11]],
                  OutPutGraph[[8]][[15]], OutPutGraph[[6]][[15]]
                          + rremove("x.text"),
                          labels = c("A", "B", "C", "D", "E", "F"), font.label = 16,
                          ncol = 2, nrow = 3)
annotate_figure(no_hunt, top = text_grob("ASF Dynamics - Hunting Season", face = "bold", size = 16))

grid.arrange(
  grobs = list(OutPutGraph[[6]][[2]], OutPutGraph[[8]][[2]],
               OutPutGraph[[6]][[6]], OutPutGraph[[8]][[6]],
               OutPutGraph[[6]][[8]], OutPutGraph[[8]][[8]]),
  widths = c(2, 2),
  layout_matrix = rbind(c(1, 2),c(3, 4), c(5, 6)), 
  top = textGrob("ASF Dynamics - Hunting Season", gp = gpar(fontsize = 20, font = 3))
)

grid.arrange(
  grobs = list(OutPutGraph[[6]][[8]],OutPutGraph[[8]][[8]],
               OutPutGraph[[6]][[8]],OutPutGraph[[8]][[8]]),
  widths = c(2, 2),
  layout_matrix = rbind(c(1, 2),c(3, 4)),
  top = textGrob("ASF Dynamics" ,gp = gpar(fontsize = 20,font = 3))
)

OutPutGraph[[5]][[5]]



# Survival Plot -----------------------------------------------------------

# S<-Surv(EpiDuration[ ,2],event = rep(1, 100))
# plot(survfit(S ~ 1))
EDT=NULL
for (i in 1:4){
 
  EDT <-bind_rows(EDT, data.frame(OutPutGraph[[i+4]][[9]],Sc=i))
                          
}

survfit(Surv(EpDuration , event = rep(1, nrow(EDT))) ~Sc, data = EDT) -> a
plot(a, col = 1:i, main = "Scenarios Comparison", xlab = "Time", ylab = "Probability of infection")
legend('topright', legend = c("Pyrénées-Atlantiques/Sans Chasse",
                                   "Pyrénées-Atlantiques/Avec Chasse",
                                   "Frontière Franco-Belge/Sans Chasse",
                                   "Frontière Franco-Belge/Avec Chasse"), lty = 1, col = 1:i)

survPlot<-ggsurvplot(survfit(Surv(EpDuration,rep(1, nrow(EDT)))~Sc,
                             data=EDT),
                     cumevents=FALSE,#risk.table = "nrisk_cumevents",
                     legend.labs=c("Pyrénées-Atlantiques/Sans Chasse",
                                   "Pyrénées-Atlantiques/Avec Chasse",
                                   "Frontière Franco-Belge/Sans Chasse",
                                   "Frontière Franco-Belge/Avec Chasse"),
                     title=paste("Analyse de la Durée de l'épizootie"),
                     legend.title="",legend=c(0.8,0.8),xlab="Temps post-introduction",
                     ylab="Probabilité de persistence",
                     risk.table.y.text = FALSE,pval = TRUE,
                     events.y.text=FALSE,surv.median.line = "hv"#ncensor.plot=TRUE,  
                    )
