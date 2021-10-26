library(rgeos)
library(magrittr)
library(sp)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tibble)
library(reshape2)
library(survival)
library(grid)
library(gridExtra)
library(survminer)

# Scenarios ---------------------------------------------------------------

#The good Outcomes are the ones inside the folders with the number of the scenarios: Scenario_1/1.RData
#Scenarios 1-4 Infection in Winter (Day 750) and 5-8 Infection in summer (Day 900)
#Outcomes of the population without infections are outside file

### Table with Scenarios used for WBModel

Scenarios <- data.frame(Scenario         = c(1:12),
                        Zone             = rep(c(rep("Pyrenees-Atlantiques_Department", 2), 
                                              rep("France-Belgium_border", 2)), 3),
                        ProbHarvest      = rep(c(0, 0.40), 6),
                        ProbHarvestAM    = rep(c(0, 0.60), 6),
                        ProbMovHunt      = rep(c(0, 0.30), 6),
                        GroupMoveHunt    = rep(c("No", "Yes"),6),
                        TimeSeedInf      = c(rep(3650, 4), rep(730, 4), rep(900, 4)), #No infection/Winter/Summer
                        MaxIterations    = c(rep(10, 4), rep(100, 8)),
                        MaxDays          = c(rep(5*365, 4), rep(8*365, 8)),
                        MaxBCap          = rep(8, 12),
                        SurvivalProbAdF  = rep(0.85, 12),                   # Yearly Survival Prob for adults Kueling et al. (2013)    
                        SurvivalProbAdM  = rep(0.75, 12),                   # Yearly Survival Prob for adults Kueling et al. (2013) 
                        SurvivalProbSAdF = rep(0.75, 12),                   # Yearly Survival Prob for sub-adults (we use adults values as they are closer to Lange et al. 2015) Kueling et al. (2013)
                        SurvivalProbSAdM = rep(0.70, 12)
                        )

# List of Output Graphs ---------------------------------------------------

Scenario = 1:8
#par(mfrow = c(2, 1))
ScenarioTable <- NULL
OutPutGraph   <- NULL

for (Scen in Scenario) {
  if(Scen %in% c(1:2, 5:6, 9:10)) load("Inputs/habitats_Dep64.RData") else
  load("Inputs/habitats_EA9.RData")
  
  #Study Area
  if(Scen %in% c(1:2, 5:6, 9:10)) Study_Area = "Pyrenees-Atlantique Department" else 
                           Study_Area = "France-Belgium Border"
  
  #Details
  if(Scen %in% c(1, 3, 5, 7, 9, 11)) Details = "Without Hunting" else Details = "Hunting Season"
  
  #Season Introduction
  if(Scen %in% c(1:4)) Season = "None" 
  if(Scen %in% c(5:8)) Season = "Winter" 
  if(Scen %in% c(9:12)) Season = "Summer"

  #Surface in km2
  gArea(habitats)/1e6 -> surface
  
  #Load Scenario List Results
  #load(paste0("Outputs/scenario_", Scen, ".RData"))
  load(paste0("./Outputs/Scenario_", Scen, "/", Scen, ".RData"))
  do.call(rbind, res) -> results
  
  ### Divide type of results
  res_dayout <- results[ ,c("DayOutInfMat", "DayOutPopMat", "DayOutAniMat", "DayOutCarcMat", "DayOutInfAniMat")]
  res_newinf <- results[ ,c("NewInfAnimals", "NewInfCarcass", "NewInfGroups", "InfectedPerCells")]
  
  #To check if Iters are messed up
  # which(sapply(res, function(x){
  #   x <- x$NewInfAnimals %>% data.frame
  #   length(unique(x$Iter)) != 1
  # }))
  
  ### Extract Results from main list
  apply(res_dayout, 2, identity) -> z
  
  Groups      <- do.call(cbind, z$DayOutPopMat)
  Groupquant <- Groups %>% data.frame %>%
                apply(., 1, quantile, probs = c(0.5, 0.025, 0.975)) %>% t(.)
  
  Population  <- do.call(cbind, z$DayOutAniMat)
  Popquant    <- Population %>% data.frame %>% 
                  apply(., 1, quantile, probs = c(0.5, 0.025, 0.975)) %>% t(.)
  
  Density     <- Population %>% data.frame %>% 
                 apply(., 1, divide_by, surface) %>% t(.)
  Denquant    <- Density %>% data.frame() %>%
                 apply(., 1, quantile, probs = c(0.5, 0.025, 0.975)) %>% t(.)
  
  PopPlot <- Population %>% data.frame %>% mutate(Time = seq(nrow(Population)))
  PopPlot <- melt(PopPlot, id = "Time")
  
  DenPlot <- PopPlot %>% mutate(density = value/surface)

 GroupPlot <- Groups %>% data.frame %>% mutate(Time = seq(nrow(Groups)))
 GroupPlot <- melt(GroupPlot, id = "Time")
 
  Infected <- do.call(cbind, z$DayOutInfAniMat)
  Infquant <- Infected %>% data.frame %>%
              apply(., 1, quantile, probs = c(0.5, 0.025, 0.975)) %>% t(.)
  
  InfCum  <- Infquant %>% data.frame %>% mutate(cum_median = cumsum(Infquant[ ,1]))
  
  InfPlot   <- Infected %>% data.frame %>% mutate(Time = seq(nrow(Infected)))
  InfPlot   <- melt(InfPlot, id = "Time")
  
  Carcasses  <- do.call(cbind, z$DayOutCarcMat)
  Carquant  <- Carcasses %>% data.frame %>%
               apply(., 1, quantile, probs = c(0.5, 0.025, 0.975)) %>% t(.)
  
  CarcCum  <- Carquant %>% data.frame %>% mutate(cum_median = cumsum(Carquant[ ,1]))
  
  CarPlot   <- Carcasses %>% data.frame %>% mutate(Time = seq(nrow(Carcasses)))
  CarPlot   <- melt(CarPlot, id = "Time")
  
  # apply(res_newinf, 2, function(x) do.call(rbind, x)) -> z
  # InfPerCells <- z$InfectedPerCells
  
  if (Scen >4 ){
    apply(res_newinf, 2, function(x) do.call(rbind, x)) -> z
    InfPerCells <- z$InfectedPerCells
  }
  
  ### Calculate Epidemic Duration
  
 # EpiDuration <- lapply(res, function(x){
 #    x= x$NewInfAnimals %>% data.frame
 #    c(unique(x$Iter), max(x$gTime) - 750)
 #    }) %>% do.call(rbind,.)
 
  Sim = NULL
  
  for (j in 1:Scenarios$MaxIterations[Scenarios$Scenario == Scen]) 
  Sim = c(Sim, rep(j, Scenarios$MaxDays[Scenarios$Scenario == Scen]))
  CarPlot$Sim = Sim
  InfPlot$Sim = Sim
  EpiDuration = apply(rbind(tapply(CarPlot$value, CarPlot$Sim, function(x) max(which(x>0))),
                            tapply(InfPlot$value, InfPlot$Sim, function(x) max(which(x>0)))), 2, max) - Scenarios[Scenarios$Scenario == Scen, "TimeSeedInf"]
  EpiDuration = tibble(Iter = names(EpiDuration), EpDuration = as.numeric(EpiDuration))
  
 ### All Iterations Graphic 
 PPopulationA <- ggplot(PopPlot, aes(x = Time, y = value, color = variable)) +
                ggtitle(label = Study_Area, subtitle = "Population") +
                ylab("Number of wild boars") + xlab("Time (days)")+
                geom_line() + theme(legend.position = "none")
    
 ### Median and quantile plots        
 PPopulationQ <- ggplot(PopPlot, aes(x = Time, y = value)) +
                ggtitle(label = Study_Area, subtitle = paste("Population")) +
                ylab("Number of wild boars") + xlab("Time (days)")+
                geom_line(colour = "cornflowerblue", alpha = 0.6) + theme(legend.position = "none") +
                stat_summary(fun.min = function(z) { quantile(z, 0.25) },
                             fun.max = function(z) { quantile(z, 0.75) },
                             fun = median, geom = "smooth", colour = "darksalmon", size = 1)
  
  PDensityA <- ggplot(DenPlot, aes(x = Time, y = density, color = variable)) +
              ggtitle(label = Study_Area, subtitle = paste("WB Density")) +
              geom_line() + theme(legend.position = "none")

  PDensityQ <- ggplot(DenPlot, aes(x = Time, y = density)) +
                 ggtitle(label = Study_Area, subtitle = paste("WB Density")) +
                 geom_line(colour = "cornflowerblue", alpha = 0.6) + theme(legend.position = "none") +
                 stat_summary(fun.min = function(z) { quantile(z, 0.25) },
                              fun.max = function(z) { quantile(z, 0.75) },
                              fun = median, geom = "smooth", colour = "coral", size = 1)
  
  PGroupsA <- ggplot(GroupPlot, aes(x = Time, y = density, color = variable)) +
              ggtitle(label = Study_Area, subtitle = paste("Wild boar Groups")) +
              ylab("Number of groups") + xlab("Time (days)") +
              geom_line() + theme(legend.position = "none")
  
  PGroupsQ <- ggplot(GroupPlot, aes(x = Time, y = value)) +
                ggtitle(label = Study_Area,  paste("Wild boar Groups")) +
                ylab("Number of groups") + xlab("Time (days)") +
                geom_line(colour = "cornflowerblue", alpha = 0.6) + theme(legend.position = "none") +
                stat_summary(fun.min = function(z) { quantile(z, 0.25) },
                             fun.max = function(z) { quantile(z, 0.75) },
                             fun = median, geom = "smooth", colour = "darksalmon", size = 1)

  PInfectedA <- ggplot(InfPlot, aes(x = Time, y = value, color = variable)) +
                ggtitle(label = Study_Area, subtitle = paste("Infected wild boars")) +
                ylab("Infected wils boars") + xlab("Time (days)") +            
                geom_line() + theme(legend.position = "none")
              
  PInfectedQ <- ggplot(InfPlot, aes(x = Time, y = value)) +
                ggtitle(label = Study_Area, subtitle = paste("Infected wild boars")) +
                ylab("Infected wils boars") + xlab("Time (days)") +               
                geom_line(colour = "cornflowerblue", alpha = 0.6) + theme(legend.position = "none") +
                stat_summary(fun.min = function(z) { quantile(z,0.25) },
                  fun.max = function(z) { quantile(z, 0.75) },
                  fun = median, geom = "smooth", colour = "brown2", size = 1)
  
  PInfCum   <- ggplot(InfPlot, aes(x = Time, y = Inf_Cum)) +
    ggtitle(label = Study_Area, subtitle = paste("Cumulative infected wild boars"))
 
  PCarcassA <- ggplot(CarPlot, aes(x = Time, y = value, color = variable)) +
               ggtitle(label = Study_Area, subtitle = paste("Infected carcasses")) +
               ylab("Infected carcasses") + xlab("Time (days)") +              
               geom_line() + theme(legend.position = "none")
  
  PCarcassQ <- ggplot(CarPlot, aes(x = Time, y = value)) +
               ggtitle(label = Study_Area, subtitle = paste("Infected carcasses")) +
               ylab("Infected carcasses") + xlab("Time (days)") +            
               geom_line(colour = "cornflowerblue", alpha = 0.6) + theme(legend.position = "none") +
               stat_summary(fun.min = function(z) { quantile(z,0.25) },
                            fun.max = function(z) { quantile(z,0.75) },
                            fun = median, geom = "smooth", colour = "darkgoldenrod2", size = 1)
 
 ScenarioTable <- rbind(ScenarioTable, c(Scenario   = Scen,
                                         Study_Area = Study_Area,
                                         Details    = Details,
                                         Season     = Season,
                                         Population = Popquant[730, 1],
                                         Density    = Denquant[730, 1],
                                         Infected   = Infquant[1095, 3],
                                         InfectCum  = InfCum  [2920, 4],
                                         Carcasses  = Carquant[1095, 3],
                                         CarcassCum = CarcCum [2920, 4])
                                         )
 
 OutPutGraph[[Scen]] <- list(PPopulationA, PPopulationQ, Popquant,
                             PDensityA,    PDensityQ,    Denquant,
                             PGroupsA,     PGroupsQ,     Groupquant,
                             PInfectedA,   PInfectedQ,   Infquant,   InfCum,
                             PCarcassA,    PCarcassQ,    Carquant,   CarcCum,
                             EpiDuration)
   
}

# OutPutGraph[[Scenario]][Graph/Table]
# Plot all Graphs per Scenario

ScenarioTable

OutPutGraph[[2]][11]
view(OutPutGraph[[1]][[2]])
summary(OutPutGraph[[8]][[12]])

#No Infection Results

#Results without hunting
# grid.arrange(
#   grobs = list(OutPutGraph[[7]][[2]],  OutPutGraph[[5]][[2]],
#                OutPutGraph[[7]][[11]], OutPutGraph[[5]][[11]],
#                OutPutGraph[[7]][[15]], OutPutGraph[[5]][[15]]),
#           widths = c(2, 2),
#           layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 6)),
#           labels = c("A", "B", "C", "D", "E", "F"),
#           top = textGrob("ASF Dynamics - No Hunting" ,gp = gpar(fontsize = 20, font = 3))
# )

no_hunt <- ggarrange(OutPutGraph[[7]][[2]],  OutPutGraph[[5]][[2]],
                     OutPutGraph[[7]][[11]], OutPutGraph[[5]][[11]],
                     OutPutGraph[[7]][[15]], OutPutGraph[[5]][[15]] + rremove("x.text"),
                     labels = c("A", "B", "C", "D", "E", "F"), font.label = 16,
                     ncol = 2, nrow = 3)
            annotate_figure(no_hunt,
                            top = text_grob("ASF Dynamics - No Hunting", face = "bold", size = 16))

hunt <- ggarrange(OutPutGraph[[8]][[2]],  OutPutGraph[[6]][[2]],
                  OutPutGraph[[8]][[11]], OutPutGraph[[6]][[11]],
                  OutPutGraph[[8]][[15]], OutPutGraph[[8]][[15]] + rremove("x.text"),
                  labels = c("A", "B", "C", "D", "E", "F"), font.label = 16,
                  ncol = 2, nrow = 3)
annotate_figure(hunt,
                top = text_grob("ASF Dynamics - Hunting Season", face = "bold", size = 16))


# Survival Plot -----------------------------------------------------------

# S<-Surv(EpiDuration[ ,2],event = rep(1, 100))
# plot(survfit(S ~ 1))

EDT = NULL
for (i in 1:4){
  EDT <- bind_rows(EDT, data.frame(OutPutGraph[[i + 4]][[18]], Sc = i))
}

survfit(Surv(EpDuration, event = rep(1, nrow(EDT))) ~Sc, data = EDT) -> a
plot(a, col = 1:i, main = "Scenario Comparison", xlab = "Time", ylab = "Probability of infection")
legend('topright', legend = c("DPA / No Hunting",
                              "DPA / Hunting Season",
                              "FBB / No Hunting",
                              "FBB / Hunting Season"), lty = 1, col = 1:i)

survPlot <- ggsurvplot(survfit(Surv(EpDuration, rep(1, nrow(EDT)))~Sc,
                         data = EDT),
                         cumevents = FALSE,
                         legend.labs=c("DPA / No Hunting",
                                       "DPA / Hunting Season",
                                       "FBB / No Hunting",
                                       "FBB / Hunting Season"),
                         title = paste("Epidemic Duration Analysis"),
                         legend.title = "", legend = c(0.8, 0.8), xlab=  "Time post-introduction",
                         ylab = "Persistence Probability",
                         risk.table.y.text = FALSE, pval = TRUE,
                         events.y.text = FALSE, surv.median.line = "hv"  
)
survPlot