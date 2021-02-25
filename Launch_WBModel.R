######################
#### Wild Bors Model ....
#### date
#### Luis Salazar and Mathieu Andraud
######################

setwd("//splffiler/ueqac/Luis/Script_12-02-21/Dep64")
#Load libraries

library("flux")

# Load functions
source(file = "\\\\splffiler/ueqac/Luis/Model/WBModel.R")

# Define Scenarios
Scenarios <- data.frame(Sc = c(1:4),
                        ProbHarvest   = c(0, 0.40, 0, 0.40),
                        ProbHarvestAM = c(0, 0.60, 0, 0.60),
                        ProbMovHunt   = c(0, 0.30, 0, 0.30),
                        GroupMoveHunt = c("No", "Yes", "No", "Yes")
)

t0 <- Sys.time()

for (Scen in seq(nrow(Scenarios))) {path = paste("./Outputs/Scenario",Scen, sep = '_')
if  (!dir.exists(path))  dir.create(path)

MaxIterations = 100
ProbHarvest   = Scenarios[Scenarios$Sc == Scen, "ProbHarvest"]
ProbHarvestAM = Scenarios[Scenarios$Sc == Scen, "ProbHarvestAM"]
ProbMovHunt   = Scenarios[Scenarios$Sc == Scen, "ProbMovHunt"]
nbcores       = 34

if(Scen %in%c(3:4)){WBModel(MaxIterations = MaxIterations, 
                            ProbHarvest   = ProbHarvest,
                            ProbHarvestAM = ProbHarvestAM,
                            ProbMovHunt   = ProbMovHunt,
                            FileWildBoarMat   = "EA9.csv",
                            runID             = paste0("WB_Model", "_FraBel_"),
                            nbcores           = nbcores)
} else
  WBModel(MaxIterations = MaxIterations, 
          ProbHarvest   = ProbHarvest,
          ProbHarvestAM = ProbHarvestAM,
          ProbMovHunt   = ProbMovHunt,
          runID = paste0(path, "/", Scen),
          nbcores = nbcores)
}

print(Sys.time()-t0)
