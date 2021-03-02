######################
#### Wild Boar Model Main Launch Panel
#### February 2021
#### Luis Salazar and Mathieu Andraud
######################

setwd("//splffiler/ueqac/Luis/Model/GIT/WildBoar_Model")

#Load libraries
library("flux")

# Load functions
source(file = "//splffiler/ueqac/Luis/Model/GIT/WildBoar_Model/WBModel.R")

# Define Scenarios
Scenarios <- data.frame(Sc = c(3:4),
                        ProbHarvest   = c(0, 0.40),   
                        ProbHarvestAM = c(0, 0.60),     
                        ProbMovHunt   = c(0, 0.30),     
                        GroupMoveHunt = c("No", "Yes"))

t0 <- Sys.time()

for (Scen in Scenarios$Sc) {
  path = paste("//splffiler/ueqac/Luis/Model/GIT/WildBoar_Model/Outputs/Scenario", Scen, sep = '_')
if  (!dir.exists(path))  dir.create(paste0("Scenario", Scen, sep = "_"))

MaxIterations = 10
ProbHarvest   = Scenarios[Scenarios$Sc == Scen, "ProbHarvest"]
ProbHarvestAM = Scenarios[Scenarios$Sc == Scen, "ProbHarvestAM"]
ProbMovHunt   = Scenarios[Scenarios$Sc == Scen, "ProbMovHunt"]
nbcores       = 34

WBModel(MaxIterations    = MaxIterations, 
                             ProbHarvest     = ProbHarvest,
                             ProbHarvestAM   = ProbHarvestAM,
                             ProbMovHunt     = ProbMovHunt,
                             savingpath      = path
                             )
 }


print(Sys.time()-t0)
