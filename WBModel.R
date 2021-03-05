##########################################################################################################
####################                       Wild Boar Model         ##################

#setwd("X:/Luis/Model")


# WB Model Function -------------------------------------------------------

WBModel <- function( MaxIterations    = 100,
                     MaxDays          = 365*8,
                     HuntingSeason    = (9*7+365)-(40*7),             # Hunting season goes from October to February
                     Detailed         = TRUE,
                     InitialOccCells  = 500,
                     MaxMortProbAd    = 1 - ((0.4)^(1/365)),          # Minimum Survival Prob for adults Lange et al. (2015), we recalculate per day 
                     MaxMortProbSAd   = 1 - ((0.4)^(1/365)),          # Minimum Survival Prob for Sub-adults Lange et al. (2015) we recalculate per day 
                     MaxMortProbPig   = 1 - ((0.1)^(1/365)),          # Minimum Survival Prob for piglets Lange et al. (2015) we recalculate per day
                     SurvivalProbAdF  = 0.9986^365,                   # Yearly Survival Prob for adults Kueling et al. (2013)    
                     SurvivalProbAdM  = 0.9986^365,                   # Yearly Survival Prob for adults Kueling et al. (2013) 
                     SurvivalProbSAdF = 0.9988^365,                   # Yearly Survival Prob for sub-adults (we use adults values as they are closer to Lange et al. 2015) Kueling et al. (2013)
                     SurvivalProbSAdM = 0.9986^365,                   # Yearly Survival Prob for sub-adults(we use adults values as they are closer to Lange et al. 2015) Kueling et al. (2013)    
                     SurvivalProbPigF = 0.50,                         # Yearly Survival Prob for piglets Lange et al. (2015) 
                     SurvivalProbPigM = 0.50,                         # Yearly Survival Prob for piglets Lange et al. (2015), reduced to consider higher survival of females
                     ProbHarvest      = 0.40,                         # A wild boar has 40% chance to be harvested
                     ProbHarvestAM    = 0.60,                         # An adult male has 70% of being harvested Toigo (2008)
                     ProbMovHunt      = 0.30,                         # Probability to move due to hunting
                     AgeProbab        = c(0.39, 0.24, 0.15, 0.09, 0.06, 0.03, 0.02, 0.01, 0.01), # Age probabilities from Lange et al. (2015)
                     ProbSplitSDS     = 1/7,     #cumsum(rep(1/7,7))       # Probability of short distance split if all conditions are satisfied for a group to split per day during the week of splitting (female splitting)
                     ProbSplitMSA     = 1/(6*7), #cumsum(rep(1/(6*7),6*7)) # Probability of male splitting per group per day of the week (7days) of the total splitting period (6 weeks)
                     SettlingProb     = c(0.75, 1),                   # probability to settle for a male in a habitat cell if it was only good habitat (0.5) or good habitat with females (1.0)
                     ProbMoveHS       = 1/HuntingSeason,              # daily probability of group splitting during the whole Hunting season
                     FemPropGroup     = mean(c(14.7,18.1,12.8,14.6)), # Proportion of females within a group Merta et al. (2015)  
                     MalPropGroup     = mean(c(9.2,9,5.7,7.9)),       # Proportion of males within a group Merta et al. (2015)
                     SubAPropGroup    = mean(c(33.4,34.3,22.8,24.1)), # Proportion of sub-adults within a group Merta et al. (2015)
                     PigletPropGroup  = mean(c(42.7,38.6,58.7,53.4)), # Proportion of piglets within a group Merta et al. (2015)
                     
                     ### Distance threshold for dispersion in meters
                     
                     DistThresholdHS      = 9000,
                     DistThresholdFemales = 24000,
                     
                     ### ASF related parameters
                     
                     WGProbInf        = 1 - (1 - 0.05)^(1/7),         # Direct contact Lange et al. (2015EFSA) and Lange and Thulke (2017)
                     BGProbInf        = 0.1*(1 - (1 - 0.05)^(1/7)),
                     CarcProbInf      = 1 - (1 - 0.02) ^(1/7),        # Carcass Lange et al. (2015EFSA) and Lange and Thulke (2017)
                     ProbAccesCarc    = 0.9,                          # Lange and Thulke (2017)
                     ProbDeathInf     = 0.95,                         # Blome (2012)  
                     CarcDaySur       = 28,                           # Lange and Thulke (2017)
                     DieInConEdge     = 0,                            # Lange et al. (2015EFSA)
                     TimeSeedInf      = 750,                          # Day 730 (2nd year)
                     NumbGSeedInf     = 1,                            # Seed in one individual randomly
                     MinBCap          = 1,
                     ModeBCap         = 4,
                     MaxBCap          = 6,
                     AREA             = 1,
                     
                     #Reproduction Parameters
                     
                     # RepProbList      = list(Week = 1:26,
                     #                         Prob = c(0, 0, 0, 0.01, 0.05, 0.1, 0.2, 0.23, 0.28, 0.325, 0.35, 0.325, 0.3, 0.28, 0.25, 0.22, 0.18, 0.17, 0.16, 0.15,
                     #                                0.14, 0.13, 0.12, 0.08, 0.04, 0.0)),
                     NumOffsprProbList = list(Number = 0:8,
                                             Prob = c(0.01, 0.07, 0.16, 0.25, 0.25, 0.16, 0.07, 0.02, 0.001)),
                     DirectionList     =  list(c(6, 7, 8), c(8, 10, 13), c(11, 12, 13), c(6, 9, 11)), ## North, ## East, ## South, and ## West
                     
                     HabitatProb       = c(0.75, 0.20, 0.05),         # The probability to select a new pixel while male walking (we model that they prefer to walk in a habitat cell - Hab Cat c(1,2,3)
                     
                     #Input and ID
                     
                     FileWildBoarMat   = "Inputs/EA9.csv",
                     runID             = paste0("WB_Model", "_FraBel"),
                     nbcores           = 2,
                     savingpath        = "//splffiler/ueqac/Luis/Model/GIT/WildBoar_Model/Outputs/"
){
  
  #Source Initialization (Wild boar matrix) and Animal Processes (Ageing and natural mortality)
  
  source("//splffiler/ueqac/Luis/Model/GIT/WildBoar_Model/Initialization.R")
  source("//splffiler/ueqac/Luis/Model/GIT/WildBoar_Model/AnimalProcesses.R")

# Raster  and Distance Definition -----------------------------------------
 
  WBMat    <- DefineWbMat(FileWildBoarMat, MinBCap, ModeBCap, MaxBCap)
  coords   <- as.matrix(WBMat[ ,c('Lon', 'Lat')])
  Distance <- as.matrix(dist(coords))
  
  
  # Gamma Distribution for Reproduction Parameters --------------------------
  
  probRepro <- 0
  for(j in 1:51) probRepro <- c(probRepro, trapz(j:(j+1), dgamma(j:(j+1), shape = 6.4, scale = 2.1)))
  probRepro[1] = 1-sum(probRepro)
  RepProbList      = list(Week = 1:52,
                          Prob = probRepro)
  

# Mortality Probabilities -------------------------------------------------

  ProbMortAdF  <- 1 - (SurvivalProbAdF)^(1/365)
  ProbMortAdM  <- 1 - (SurvivalProbAdM)^(1/365)
  ProbMortSAdF <- 1 - (SurvivalProbSAdF)^(1/365)
  ProbMortSAdM <- 1 - (SurvivalProbSAdM)^(1/365)
  ProbMortPigF <- 1 - (SurvivalProbPigF)^(1/365)
  ProbMortPigM <- 1 - (SurvivalProbPigM)^(1/365)
  ProbHunted   <- 1 - (1 - ProbHarvest)^(1/HuntingSeason)
  ProbHuntedAM <- 1 - (1 - ProbHarvestAM)^(1/HuntingSeason) 
  
  # Start of Iterations -----------------------------------------------------
  
  library(doParallel)
  cl1 <- parallel::makeCluster(nbcores)
  doParallel::registerDoParallel(cl1)
  clusterEvalQ(cl1, library("flux", "dplyr"))
  res = foreach (Iter = 1:MaxIterations) %dopar% {
    
    source("//splffiler/ueqac/Luis/Model/GIT/WildBoar_Model/Initialization.R")
    source("//splffiler/ueqac/Luis/Model/GIT/WildBoar_Model/AnimalProcesses.R")
    set.seed(Iter)
    
    # Model Outcomes ----------------------------------------------------------
    
    DayOutInfMat    <- matrix(0, ncol = 1, nrow = MaxDays)
    DayOutCarcMat   <- matrix(0, ncol = 1, nrow = MaxDays)
    DayOutInfAniMat <- matrix(0, ncol = 1, nrow = MaxDays)
    DayOutPopMat    <- matrix(0, ncol = 1, nrow = MaxDays)
    DayOutAniMat    <- matrix(0, ncol = 1, nrow = MaxDays)
    NewInfGroups    <- matrix(numeric(0), ncol = 4)
    NewInfAnimals   <- matrix(numeric(0), ncol = 3)
    NewInfCarcass   <- matrix(numeric(0), ncol = 3)
    InfectedPerCells <- matrix(numeric(0), ncol = 5)
    EpDuration    <- rep(0, MaxIterations)
    FreqRelapse   <- rep(0, MaxIterations)
    cumDeath      <- rep(0, MaxIterations)
    Fpopcount     <- matrix(0, nrow = MaxDays, ncol = nrow(WBMat))
    Finfcount     <- matrix(0, nrow = MaxDays, ncol = nrow(WBMat))
    Fcarcount     <- matrix(0, nrow = MaxDays, ncol = nrow(WBMat))
    Fimmunecount  <- matrix(0, nrow = MaxDays, ncol = nrow(WBMat))
    
# Initialize Population Matrix --------------------------------------------

    PopMatWB <- InitPop(InitialOccCells,
                        Fem     = FemPropGroup,
                        Mal     = MalPropGroup,
                        SubA    = SubAPropGroup,
                        Piglet  = PigletPropGroup,
                        AgeProb = AgeProbab, WBMat, MinBCap=MinBCap, ModeBCap=ModeBCap, MaxBCap=MaxBCap)
    
    GroupsToSplit   <- matrix(numeric(0), ncol = 4)     
    ## Store the groups that have split this year
    SplittedGroups  <- numeric(0)
    MovedGroups     <- numeric(0)
    cumDeathPar     <- 0
    cumDeathInf     <- 0
    Year            <- 1
    gTime           <- 0
    Criteria        <- TRUE
    TMPOutYesterday <- FALSE
    OnlyOnce        <- TRUE

# Daily Loop --------------------------------------------------------------

    while(Criteria & (gTime < MaxDays)){
      
      gTime <- gTime + 1
      if(gTime > 365) Year <- ceiling(gTime/365)
      print(table(PopMatWB[ ,9]))
      PopMatWB <- Ageing(gTime, PopMatWB)
      
      if(gTime %in% c(1, (366 * 1:(MaxDays/365)))){

    ### Select females that will reproduce during the year
        PopMatWB        <- PopMatWB[order(PopMatWB[ ,2], PopMatWB[ ,3], PopMatWB[ ,5], decreasing = T), ]
        AnimalsWInG     <- do.call(c, sapply(unique(PopMatWB[ ,2]), function(x) 1:sum(PopMatWB[ ,2] == x)))
        BreedCapCells   <- c(WBMat[PopMatWB[ ,7], 15], rep(0, sum(PopMatWB[ ,7] == 0)))
        # Animals are allowed to breed when they are females and the group size is less than the maximum allowed size (breeding capacity)
        #(size that can be maintained by the pixel habitat)
        # and the animals above the breeding capacity based on age.
        AniAllowBreed   <- AnimalsWInG <= BreedCapCells & PopMatWB[ ,3] == 1 & PopMatWB[ ,9] < 2
        Breedingdays <- sample(0:6, sum(AniAllowBreed), replace = TRUE)
        PopMatWB[AniAllowBreed, 6] <- (sample(RepProbList$Week, sum(AniAllowBreed), rep = T, prob = RepProbList$Prob)*7)+ Breedingdays + ((Year-1)*365)
      }      
        
  ToDieNormal <- Mortality(PopMatWB, ProbMortPigF, ProbMortPigM,
                       ProbMortSAdF, ProbMortSAdM, ProbMortAdF, ProbMortAdM)

# Hunting Season ----------------------------------------------------------

  if(gTime %in% ((40*7 + ((Year - 1)*365)):((9*7 + ((Year)*365))))){
    
    #ToDieHunted  <- which(rbinom(nrow(PopMatWB), 1, ProbHunted) == 1) #Without Hunting Pro for AM
    ToDieHunted   <- which(!(PopMatWB[ ,3] ==  0 & PopMatWB[ ,4] == 3) & PopMatWB[ ,4] != 1)
    ToDieHuntedAM <- which(PopMatWB[ ,3] == 0 & PopMatWB[ ,4] == 3)
    ## Hunting without piglets and with a different hunting probability for adult males
    ToDieHunted   <- ToDieHunted[rbinom(length(ToDieHunted), 1, ProbHunted) == 1]
    ToDieHuntedAM <- ToDieHuntedAM[rbinom(length(ToDieHuntedAM), 1, ProbHuntedAM) == 1]
    Hunted        <- c(ToDieHuntedAM, ToDieHunted)
    
    #to check if there are no piglets and male adults
    #table(data.frame(PopMatWB[which(!(PopMatWB[ ,3] ==  0 & PopMatWB[ ,4] == 3)  &  PopMatWB[ ,4] != 1),c("Sex", "Age_Cat")]))
    
    ## Movement caused by hunting during Hunting Season
    ## Groups that will split due to hunting
    GroupsMoveHunt <- sort(unique(PopMatWB[Hunted ,2]))
    #GroupsMoveHunt<- GroupsMoveHunt[GroupsMoveHunt[ ,2] >= 1, , drop = FALSE]
    GroupsMoveHunt <- GroupsMoveHunt[rbinom(length(GroupsMoveHunt), 1, ProbMovHunt) == 1]
    
    GroupsMoveShD  <- sapply(PopMatWB[match(GroupsMoveHunt, PopMatWB[ ,2]), 7], function(x){
      tempo1     <- which(Distance[x, ] <= DistThresholdHS)
      any(WBMat[tempo1, 5] == 1)})
    
    GroupsMoveNumSD <- cbind(GroupsMoveHunt,  GroupsMoveShD)
    GroupsMoveNumSD <- GroupsMoveNumSD[GroupsMoveNumSD[ ,2] == 1, , drop = FALSE]
    
    ## Keep track of the groups that have moved due to hunting
    PopMatWB[PopMatWB[ ,2] %in% GroupsMoveNumSD[ ,1], 16] <- 2
    
    ## Make Short distance happened
    if(dim(GroupsMoveNumSD)[1] > 0){
      OriginalPixel     <- PopMatWB[match(GroupsMoveNumSD[ ,1], PopMatWB[ ,2]), 7]
      
      TargetPixel <- numeric(0)
      for(xx in OriginalPixel){
        tempo1     <- which(Distance[xx,] <= DistThresholdHS)
        tempo2     <- tempo1[WBMat[tempo1, 5] == 1 & !(tempo1 %in% PopMatWB[ ,7])]
        if(length(tempo2) > 1)  tempo3 <- sample(tempo2, 1)
        if(length(tempo2) == 1) tempo3 <- tempo2
        if(length(tempo2) == 0) tempo3 <- 0
        TargetPixel <- c(TargetPixel, tempo3)
      }
      
      GroupsMoveNumSD <- GroupsMoveNumSD[TargetPixel > 0, , drop = FALSE]
      OriginalPixel   <- OriginalPixel[TargetPixel > 0]
      TargetPixel     <- TargetPixel[TargetPixel > 0]
      GroupsMoveNumSD <- cbind(GroupsMoveNumSD, TargetPixel) 
      
      ## Allow Movement
      if(length(TargetPixel) > 0){
        ## Create a list to keep track of pixels where the group have been. 
        ## Exported on a daily basis, the list will be re-initiate every day.
        ## How many times does the group tries to move?
        PixelsMovedHS   <- as.list(matrix(0, ncol = length(TargetPixel)))
        CurrentPosition <- OriginalPixel
        
        for(i in 1:length(TargetPixel)){
          Trail <- 0
          previousEdge <- CurrentPosition[i]
          while(CurrentPosition[i] != TargetPixel[i] & Trail < 10){
            Trail      <- Trail + 1
            Edges      <- unlist(WBMat[CurrentPosition[i], 6:13])
            Edges      <- Edges[Edges > 0 & Edges != previousEdge]
            if(length(Edges) > 0){
              previousEdge <- CurrentPosition[i]
              NewPosition  <- Edges[which.min(Distance[TargetPixel[i], Edges])]
              if(length(NewPosition) > 1) NewPosition <- sample(NewPosition, 1)
              CurrentPosition[i] <- NewPosition
              # Here we keep track of the pixels where the group has moved to.
              PixelsMovedHS[[i]] <- c(PixelsMovedHS[[i]], NewPosition)
            }
          }
        }
        
        names(PixelsMovedHS) <- GroupsMoveNumSD[ ,1]
        MovedGroupsHS <- PopMatWB[ ,2] %in% GroupsMoveNumSD[ ,1]
        #IndexNotSplit    <- which(CurrentPosition != TargetPixel)
        #GroupsMoveNumSD <- GroupsMoveNumSD[-IndexNotSplit, , drop = FALSE]
        
        ### Update Home Pixel

        for(b in 1:dim(GroupsMoveNumSD)[1]){
          GroupsCanMove <-  PopMatWB[ ,2] %in% GroupsMoveNumSD[b, 1]
          NewHomePixel  <- rep(GroupsMoveNumSD[b, 3], sum(GroupsCanMove))
          PopMatWB[GroupsCanMove, 7]  <- NewHomePixel
          PopMatWB[GroupsCanMove, 8]  <- NewHomePixel
          
        }
      }
    }
        ToDieAll <- unique(c(ToDieNormal, ToDieHunted, ToDieHuntedAM)) 
      } else {
        ToDieAll <- ToDieNormal
      }

      ## Include death toll in the cumulative deaths (not infected)
      if(length(ToDieAll) > 0){
        PopMatWB <- PopMatWB[-ToDieAll, ]
        cumDeathPar <- cumDeathPar + length(ToDieAll)
      }
      
      ## Reset the movement due to hunting each year
      if(gTime == (9*7 + (Year*365) + 1)){
        PopMatWB[ ,16]    <- 0
        #SplittedGroups   <- numeric(0)
        PixelsMovedHS     <- as.list(matrix(0, ncol = 1))
      }

# Reproduction ------------------------------------------------------------
      
      ## On daily basis, from Jan to end June, check if there are animals to deliver and make them deliver
        if(gTime %in% PopMatWB[ ,6]){
        DelIndex        <- which(PopMatWB[ ,6] %in% gTime & PopMatWB[ ,9] != 3)
        NumOfSpring     <- sample(NumOffsprProbList$Number, length(DelIndex), rep = T, prob = NumOffsprProbList$Prob)
        DelIndex        <- DelIndex[NumOfSpring > 0]
        NumOfSpring     <- NumOfSpring[NumOfSpring > 0]
        if(length(NumOfSpring) > 0){
          newBorns <- cbind((max(PopMatWB[ ,1]) + 1):(max(PopMatWB[ ,1]) + sum(NumOfSpring)),
                          rep(PopMatWB[DelIndex, 2], NumOfSpring),
                          rbinom(sum(NumOfSpring), 1, 0.5),
                          1,
                          1,
                          0,
                          rep(PopMatWB[DelIndex, 7], NumOfSpring),
                          rep(PopMatWB[DelIndex, 7], NumOfSpring),
                          0,
                          0,
                          rep(PopMatWB[DelIndex, 1], NumOfSpring),
                          0,
                          0,
                          0,
                          0,
                          0)
          colnames(newBorns) = colnames(PopMatWB)
          PopMatWB <- rbind(PopMatWB, newBorns)         
          PopMatWB <- PopMatWB[order(PopMatWB[ ,2], PopMatWB[ ,3], PopMatWB[ ,5], decreasing = T), ]
          
        }
      }

# Female Group Splitting --------------------------------------------------
      
      ## Occurs only once in week 28 Kramer-Schadt et al. (2009)
  if(gTime %in% ((28*7 + ((Year - 1)*365)) + 0:6)){
    
    PopMatWB        <- PopMatWB[order(PopMatWB[ ,2], PopMatWB[ ,3], PopMatWB[ ,5], decreasing = T), ]
    AnimalsWInG     <- unlist(sapply(unique(PopMatWB[ ,2]), function(x) 1:sum(PopMatWB[ ,2] == x)))
    BreedCapCells   <- c(WBMat[PopMatWB[ ,7] ,15], rep(0, sum(PopMatWB[ ,7] == 0)))
    ## Identify the animals within the groups that will split
    ## Females sub-adults > breading capacity, did not split before and are not sick
    GroupSplitNum   <- cbind(sort(unique(PopMatWB[ ,2])), (tapply((AnimalsWInG > BreedCapCells &
                                                                     PopMatWB[ ,4] == 2 &
                                                                     PopMatWB[ ,3] == 1 &
                                                                     PopMatWB[,10] == 0 &
                                                                     PopMatWB[,9] < 2),  PopMatWB[ ,2], sum)))
    GroupSplitNum   <- GroupSplitNum[GroupSplitNum[ ,2] >= 2, , drop = FALSE]
    
    ## Identify only-male groups to allow females to join
    PixelMalesOnly  <- cbind(sort(unique(PopMatWB[ ,2])), (tapply((PopMatWB[ ,3] == 0), PopMatWB[ ,2], all)))
    PixelMalesOnly  <- PixelMalesOnly[PixelMalesOnly[ ,2] == 1, , drop = FALSE]
    
    GroupsSplitShD  <- sapply(PopMatWB[match(GroupSplitNum[ ,1], PopMatWB[ ,2]), 7], function(x) {
      tmp1 <-which(Distance[x,] <= DistThresholdFemales)
      any(WBMat[tmp1, 5] == 1 & !(tmp1 %in% PopMatWB[!PopMatWB[ ,2] %in% PixelMalesOnly[ ,1], 7]))})
    GroupSplitNumSD <- cbind(GroupSplitNum, GroupsSplitShD)
    GroupSplitNumSD <- GroupSplitNumSD[GroupSplitNumSD[ ,3] == 1, , drop=FALSE]
    
    ## Groups to split, may split during any day of the week. No need for a probability here the group splits only once a year.
    #if(dim(GroupSplitNumSD)[1]>0) GroupSplitNumSD <- GroupSplitNumSD[rbinom(dim(GroupSplitNumSD)[1],1,prob=ProbSplitSDS)==1,,drop=FALSE]
    PopMatWB[PopMatWB[ ,2] %in% GroupSplitNumSD[ ,1], 10] <- 1
    
    ## Make short term split to happen
    if(dim(GroupSplitNumSD)[1] > 0){
      OriginPixel     <- PopMatWB[match(GroupSplitNumSD[ ,1], PopMatWB[ ,2]), 7]
      TargetPixel <- numeric(0)
      for(xx in OriginPixel){
        tmp1 <- which(Distance[xx, ] <= DistThresholdFemales)
        tmp2      <- tmp1[WBMat[tmp1, 5] == 1 & !(tmp1 %in% PopMatWB[ ,7])]
        if(length(tmp2) > 1)  tmp3 <- sample(tmp2, 1)
        if(length(tmp2) == 1) tmp3 <- tmp2
        if(length(tmp2) == 0) tmp3 <- 0
        TargetPixel <- c(TargetPixel, tmp3)
      }
      
      GroupSplitNumSD <- GroupSplitNumSD[TargetPixel>0, , drop = FALSE]
      OriginPixel     <- OriginPixel[TargetPixel > 0]
      TargetPixel     <- TargetPixel[TargetPixel > 0]
      GroupSplitNumSD <- cbind(GroupSplitNumSD, TargetPixel) 
      
      ## This part of the code to allow movement
      if(dim(GroupSplitNumSD)[1] > 0){
        
        ## A list to keep track for the pixels where the pigs have been. 
        ## Make sure that this information is exported on daily basis, because the 
        ## list will be re-initiate every day splitting may happen.
        PixelsMoved <- as.list(matrix(0, ncol = length(TargetPixel)))
        CurrentPos <- OriginPixel
        for(i in 1:length(TargetPixel)){
          Trail    <- 0
          prevEdge <- CurrentPos[i]
          while(CurrentPos[i] != TargetPixel[i] & Trail < 10){
            Trail      <- Trail + 1
            Edges      <- unlist(WBMat[CurrentPos[i], 6:13])
            Edges      <- Edges[Edges > 0 & Edges != prevEdge]
            if(length(Edges) > 0){
              prevEdge    <- CurrentPos[i]
              NewPosition <- Edges[which.min(Distance[TargetPixel[i],Edges])]
              if(length(NewPosition) > 1) NewPosition <- sample(NewPosition, 1)
              CurrentPos[i] <- NewPosition
              # Here we keep track of the pixels where the pigs moved to
              PixelsMoved[[i]]<- c(PixelsMoved[[i]], NewPosition)
            }
          }
        }
        
        ## Here we assume that groups that did not find the way, did not actually split as the edges were not connected.
        IndexNotSplit   <- which(CurrentPos != TargetPixel)
        #GroupSplitNumSD <- GroupSplitNumSD[-IndexNotSplit, ,drop = FALSE]
        
        ## Make females from different groups moving to a new pixel to form a new group.
        if(dim(GroupSplitNumSD)[1] > 0){
          newGroupIDs <- (max(PopMatWB[,2]) + 1):(max(PopMatWB[ ,2]) + dim(GroupSplitNumSD)[1])
          if(sum(duplicated(GroupSplitNumSD[ ,4])) > 0){ #Target pixel
            for(l in unique(GroupSplitNumSD[ ,4])){
              TEMP  <- which(GroupSplitNumSD[ ,4] == l)
              if(length(TEMP) > 1){
                newGroupIDs[TEMP] <- newGroupIDs[TEMP[1]]
              }
            }
          }
          
          names(PixelsMoved) <- GroupSplitNumSD[ ,1]
          for(b in 1:dim(GroupSplitNumSD)[1]){
            
            ## If there is a male group already in the pixel and is not splitting, then make them one group with the new females.
            ## Notice the code above allows them to come into a new pixel only if it is empty or there are only males. 
            ## Don't worry tat IndexMalComb1 is matching with all PopMatWB[,7]
            IndexMalComb1  <- GroupSplitNumSD[b, 4] %in% PopMatWB[ ,7]
            if(IndexMalComb1){
              IndexMalComb2  <- unique(PopMatWB[PopMatWB[ ,7] == GroupSplitNumSD[b, 4], 2]) %in% GroupsToSplit[ ,1]
              if(!IndexMalComb2){
                AnimalsCanSplitF <- AnimalsWInG>BreedCapCells & 
                  PopMatWB[ ,3] == 1 & 
                  PopMatWB[ ,4] == 2 & 
                  PopMatWB[ ,9] < 2  & 
                  PopMatWB[ ,2] %in% GroupSplitNumSD[b, 1]
                MalesInPixel     <- PopMatWB[ ,3] == 0 & 
                  PopMatWB[ ,9] != 3 & 
                  PopMatWB[ ,7] == GroupSplitNumSD[b, 4] & 
                  !PopMatWB[ ,2] %in% GroupsToSplit[ ,1]
                newGroupIDsAn    <- rep(newGroupIDs[b], sum(AnimalsCanSplitF) + sum(MalesInPixel)) 
                NewHomePixel     <- rep(GroupSplitNumSD[b, 4], sum(AnimalsCanSplitF) + sum(MalesInPixel))
                PopMatWB[AnimalsCanSplitF|MalesInPixel, 2]  <- newGroupIDsAn
                PopMatWB[AnimalsCanSplitF|MalesInPixel, 7]  <- NewHomePixel
                PopMatWB[AnimalsCanSplitF|MalesInPixel, 8]  <- NewHomePixel
              }
              if(IndexMalComb2){
                AnimalsCanSplitF <- AnimalsWInG > BreedCapCells & 
                  PopMatWB[ ,3] == 1 & 
                  PopMatWB[ ,4] == 2 & 
                  PopMatWB[ ,9] < 2  & 
                  PopMatWB[ ,2] %in% GroupSplitNumSD[b, 1]
                newGroupIDsAn    <- rep(newGroupIDs[b], sum(AnimalsCanSplitF)) 
                NewHomePixel     <- rep(GroupSplitNumSD[b, 4], sum(AnimalsCanSplitF))
                PopMatWB[AnimalsCanSplitF, 2]  <- newGroupIDsAn
                PopMatWB[AnimalsCanSplitF, 7]  <- NewHomePixel
                PopMatWB[AnimalsCanSplitF, 8]  <- NewHomePixel
              }   
            }
            if(!IndexMalComb1){
              AnimalsCanSplitF <- AnimalsWInG > BreedCapCells & 
                PopMatWB[ ,3] == 1 & 
                PopMatWB[ ,4] == 2 & 
                PopMatWB[ ,9] < 2  & 
                PopMatWB[ ,2] %in% GroupSplitNumSD[b, 1]
              newGroupIDsAn    <- rep(newGroupIDs[b], sum(AnimalsCanSplitF)) 
              NewHomePixel     <- rep(GroupSplitNumSD[b, 4], sum(AnimalsCanSplitF))
              PopMatWB[AnimalsCanSplitF, 2]  <- newGroupIDsAn
              PopMatWB[AnimalsCanSplitF, 7]  <- NewHomePixel
              PopMatWB[AnimalsCanSplitF, 8]  <- NewHomePixel
            }   
          }
        }
      } 
    }
  }#Closes short distance splitting
  
      ## After the end of female splitting, we reset female splitting.
      if(gTime %in% ((28*7 + ((Year - 1)*365)) + 7)){
        PopMatWB[ ,10]    <- 0
        SplittedGroups   <- numeric(0)
        PixelsMoved      <- as.list(matrix(0, ncol = 1))
      }

# Male Group Splitting ----------------------------------------------------
      
      ## Period of splitting of males defined in Lange et al., (2012)
      ## Males may find a pixel to live in otherwise they will keep Wondering around 
      ## until they either die or find a place to live in
  
      if(gTime >= (25*7 + ((Year - 1)*365)) & gTime<=(30*7 + ((Year - 1)*365))){
        PopMatWB    <-  PopMatWB[order(PopMatWB[ ,2], PopMatWB[ ,3], PopMatWB[ ,5], decreasing = T), ]
       
        ## Identify the groups where sub-adult males may split: SA males from groups where males have not yet split this year
        NumbAdultMInG <- table(PopMatWB[ ,2], PopMatWB[ ,4], PopMatWB[ ,3])[ ,3 ,1]
        ## Calculate number of adult males
        NumbAdultMInG <-  cbind(sort(unique(PopMatWB[ ,2])), tapply(PopMatWB[ ,3] == 0 & 
                                                                    PopMatWB[ ,4] == 3, 
                                                                    PopMatWB[ ,2], sum))
        ## sub-adult males won't split if there is need for them in the group
        as.numeric(names(NumbAdultMInG)[which(NumbAdultMInG > 2)]) -> splitGroups
        table(PopMatWB[ ,2], PopMatWB[ ,3], PopMatWB[ ,4], PopMatWB[ ,12])[ ,1, 2, 1] -> tmp
        as.numeric(names(tmp)[tmp > 0]) -> tmp
        GroupSubAM <- splitGroups[splitGroups %in% tmp]
        rm(list = c("tmp", "NumbAdultMInG"))
        ## Groups that may split must not be group already splitting and have more than 2 adult males

          if (gTime == (25*7 + ((Year - 1)*365))) splitIndex = floor(length(GroupSubAM)*1/(5*7))
          if (gTime == (30*7 + ((Year - 1)*365))) splitIndex = length(GroupSubAM)  
          GroupSubAM <- sample(GroupSubAM, splitIndex)
          if(length(GroupSubAM) > 0){
         
        ## Identify the sub-adult males that may split.
        ## We must ensure that the groups have adult males before the sub-adults are split.
        ## No need for the sub-adults to find new areas if they can nourish and breed in their home
          tmpSAM <- which(PopMatWB[ ,2] %in% GroupSubAM & 
                            PopMatWB[ ,3] == 0 & 
                            PopMatWB[ ,4] == 2 & 
                            PopMatWB[ ,9] < 2)
          ## Check which ones that would split
          if(length(tmpSAM) > 0) {
            ## Make subadults to split based on probability from 25 to 75% (Truve et al., 2004).
            tmpSAM           <- tmpSAM[runif(length(tmpSAM)) < 0.5]#rbinom(length(tmpSAM), 1, prob = runif(length(tmpSAM), 0.25, 0.75)) == 1]
            AinmFromGToSplit <-  table(PopMatWB[tmpSAM, 2])
            AinmFromGToSplit <- AinmFromGToSplit[AinmFromGToSplit > 2]
            as.numeric(names(AinmFromGToSplit)) -> splitGroups
            tmpSAM           <- tmpSAM[PopMatWB[tmpSAM ,2] %in% splitGroups]
            
            if(length(tmpSAM) > 0) {# do not worry, if the dimension is > 0 then there will be at least 2 animals from one group to split ;-) check above.
              PopMatWB[tmpSAM,12] <- 1
              ## Store the group number of the groups that have splitted
              SplittedGroups <- c(SplittedGroups, splitGroups)
              ## Number of splitting males per group
              AinmFromGToSplit   <- cbind(splitGroups, AinmFromGToSplit, direction = round(runif(dim(AinmFromGToSplit)[1], 1, 4)), origin = unique(PopMatWB[rev(tmpSAM), 8]))
              ## The third column, informs about the direction
              ## Which should be depicted from the DirectionList. each group would walk in a specific direction
              ## If they cannot find a connecting edge in that direction, a random cell is then selected.
              ## Make them to prefer to move through habitat cells rather than fields.
              indexsNGroups  <- max(PopMatWB[ ,2]) + (1:dim(AinmFromGToSplit)[1])
              indexNGPigs    <- rep(indexsNGroups, AinmFromGToSplit[ ,2])
              indexNewGNum   <- sort(which(PopMatWB[ ,2] %in% AinmFromGToSplit[ ,1] & PopMatWB[ ,12] == 1), decreasing = T)
              PopMatWB[indexNewGNum, 2] <- indexNGPigs
              AinmFromGToSplit[ ,1] <- indexsNGroups
              GroupsToSplit <- rbind(GroupsToSplit, AinmFromGToSplit)      
            }
          }
        }
      
      
      if(dim(GroupsToSplit)[1] > 0){
        ## First we check that all groups still exist in the population
        Checktmp         <- which(!GroupsToSplit[ ,1] %in% PopMatWB[ ,2])
        if(length(Checktmp) > 0) GroupsToSplit[-Checktmp, , drop = FALSE]
        ## Determine the number of pixels they may move per day and the direction
        NumbPixMovesTod  <- round(runif(dim(GroupsToSplit)[1], 1, 3)) ## Lange et al. (2015)
        MovedPixMale <- as.list(matrix(0, ncol = length(NumbPixMovesTod > 0)))
        ToRemovejj <- numeric(0)
        if(any(NumbPixMovesTod > 0)){
          for(i in 1:max(NumbPixMovesTod)){
            if(any(PopMatWB[ ,12] == 1) & any(NumbPixMovesTod >= i)){
              IndexGroup   <- which(NumbPixMovesTod >= i)
              if(length(IndexGroup) > 0){
                for(jj in IndexGroup){
                  NeighboursIndex <- unlist(WBMat[GroupsToSplit[jj, 4], DirectionList[[GroupsToSplit[jj, 3]]]])
                  
                  #Find the directions for all neighbors
                  #Keep The ones in the selected area :exclude 0's
                  NeighboursIndex <- NeighboursIndex[NeighboursIndex > 0]
                  #Check the habitat cat
                  TMPIndPix <-sapply(NeighboursIndex, function(x) 
                                  x[which(WBMat[WBMat$ID %in% x, 5] == 1 &
                                  sum(PopMatWB[PopMatWB[ ,7] %in% x &
                                  PopMatWB[ ,9] != 3, 3] == 1) > 0 &
                                  sum(PopMatWB[PopMatWB[ ,7] %in% x, 9] != 3) > 0 &
                                  sum(PopMatWB[PopMatWB[ ,7] %in% x, 3] == 0) < 2)])
                  if (is.list(TMPIndPix)) TMPIndPix = do.call(c, TMPIndPix)
                  #Keep accessible and suitable ones only
                  Freepixel      <-  NeighboursIndex[!(NeighboursIndex %in% PopMatWB[ ,7])]
                  #good habitat cell and not occupied
                  TSuitableNb    <- c(TMPIndPix, Freepixel)
                  TSuitableNbDir <- TSuitableNb[TSuitableNb %in% WBMat[GroupsToSplit[jj, 4], 
                                                   DirectionList[[GroupsToSplit[jj, 3]]]]]                  
                  NbHabCat       <- WBMat[TSuitableNbDir, 5]
                  SuitableNbDir  <- which(NbHabCat > 0 & NbHabCat <= length(HabitatProb))
                  TSuitableNbDir <- TSuitableNbDir[SuitableNbDir]
                  
   
                  if(length(TSuitableNbDir) == 1) Destination <- TSuitableNbDir 
                  if(length(TSuitableNbDir) > 1)  Destination <- sample(TSuitableNbDir, 1 , 
                                                                       prob = HabitatProb[WBMat[WBMat$ID %in% TSuitableNbDir,5]])
                 
                  if(length(TSuitableNbDir) == 0){
                    NbHabCat <- WBMat[TSuitableNb, 5]
                    #SuitableNb <- which(NbHabCat>0 & NbHabCat<=3)
                    SuitableNb <- which(NbHabCat > 0 & NbHabCat <= length(HabitatProb))
                    TSuitableNb <- TSuitableNb[SuitableNb]
                    if(length(TSuitableNb) > 1)  Destination <- sample(TSuitableNb, 1, prob = HabitatProb[WBMat[WBMat$ID %in% TSuitableNb,5]])
                    if(length(TSuitableNb) == 1) Destination <- TSuitableNb
                    if(length(TSuitableNb) == 0) break
                  }
                  
                  ## Here we add the list with the moved pixels
                  MovedPixMale[[jj]] <- c(MovedPixMale[[jj]], Destination)
                  ## Is this only a suitable habitat pixel or a suitable habitat and has a female(s) with no or 1 male
                  index=Destination %in% Freepixel+2*(Destination %in% TMPIndPix )  
                  ## We decide whether they will settle in this pixel or not based on a random process
                  ## Notice here that the animals that die during splitting do not affect the splitting. 
                  ## We checked above that all groups in splitting matrix have animals in the Population matrix
                  Settle=FALSE
                  if (length(index)>0) {
                      if (index>0){
                    Settle <-runif(1)<=SettlingProb[index]
                    if(Settle){
                      indexspPigs               <- which(PopMatWB[ ,2] == GroupsToSplit[jj, 1] & PopMatWB[,12]==1 & PopMatWB[,9]!=3)
                      #newGroupIDM              <- rep(max(PopMatWB[,2])+1,length(indexspPigs))
                      #PopMatWB[indexspPigs,2]  <- newGroupIDM
                      PopMatWB[indexspPigs, 7]  <- Destination
                      PopMatWB[indexspPigs, 8]  <- Destination
                      PopMatWB[PopMatWB[ ,8]   == Destination, 2]
                      PopMatWB[indexspPigs, 12] <- 0
                      ToRemovejj                <- c(ToRemovejj, jj)
                      break
                    }
                    
                    if(!Settle){
                      indexspPigs             <- which(PopMatWB[,2]==GroupsToSplit[jj,1] & PopMatWB[,12]==1 & PopMatWB[,9]!=3)
                      PopMatWB[indexspPigs,8] <- Destination
                      GroupsToSplit[jj,4]     <- Destination
                     }
                    }
                  }
                }
              }
            }
          }
          names(MovedPixMale) <- GroupsToSplit[,1]
          if(length(ToRemovejj)>0) GroupsToSplit  <- GroupsToSplit[-ToRemovejj, ,drop = FALSE]
        }
      }
          tmp<-do.call(c,tapply(PopMatWB[,8],PopMatWB[,2],function(x) x[length(unique(x))>2]))
          if (length(tmp)>0) {
              for (i in 1:length(tmp)){ 
                b=min(PopMatWB[PopMatWB[,8]==tmp[i],2])
                PopMatWB[PopMatWB[,8]==tmp[i],2]=b
              }
            }
          }

# Seed Infection ----------------------------------------------------------

       if(TimeSeedInf == gTime){
       
       ## Here I seed randomly in 1 group and infect 1 animal (make it directly infectious). 
       ## We make sure that the groups have at least 2 animals.
       NumPGroup  <- tapply(PopMatWB[ ,2], PopMatWB[ ,2], length)
       temp1      <- NumPGroup[NumPGroup > 2]
       temp2      <- as.numeric(names(temp1))
       temp3      <- temp2[temp2 > 0]
       RandSeedG  <- sample(temp3, NumbGSeedInf)
       SampSeedAn <- sapply(RandSeedG, function(x) sample(which(PopMatWB[ ,2] == x), 1))
       PopMatWB[SampSeedAn, 9]  <- 2
       PopMatWB[SampSeedAn, 14] <- gTime   # the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5.3,1.3)
       PopMatWB[SampSeedAn, 15] <- PopMatWB[SampSeedAn,14] + round(rpert(length(SampSeedAn), 1, 5, 7))# the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5,0.3)
       
          }
  
# Carcass Persistence -----------------------------------------------------

  CarcDaySur <- 28
  ## Winter
  if(gTime %in% ((49*7) + ((Year - 1)*365)):((9*7 + ((Year)*365)))){
    CarcDaySur <- 90
  } 
  ## Summer
  if(gTime %in% ((22*7 + ((Year - 1)*365))):((35*7 + ((Year - 1)*365)))){
    CarcDaySur <- 5
    }
  
# ASF Dynamics ------------------------------------------------------------

      
      ## A - Here we model within group infection and infection from connecting edges via carcasses
   
      if(any(PopMatWB[ ,9] > 0)){
        # infGroups <- unique(PopMatWB[PopMatWB[ ,9] %in% 2:3, 2]) 
        NumInfPG     <- by(PopMatWB[ ,9] == 2, PopMatWB[ ,8], sum)
        InfCarcasses <- PopMatWB[, 9] == 3 & ((gTime - PopMatWB[, 15]) <= CarcDaySur)
        TotDead      <- by(InfCarcasses,PopMatWB[,8],sum)
        TotDeadPixel <- by(InfCarcasses,PopMatWB[,8],sum)
        tmpConEdgs   <- WBMat[unique(PopMatWB[ ,8]), c(1, 6:13)]
        probWGDCInf  <- 1 - (1 - WGProbInf)^NumInfPG
        ProbWGCInf   <- 1 - (1 - CarcProbInf)^TotDead
        ProbBGCInf   <- apply(tmpConEdgs,1,function(x) {
                                      x = x[x>0]
                                      1- (1 - BGProbInf)^(sum(NumInfPG[as.character(x[2:length(x)])], na.rm = TRUE))
                                })
          ## Estimate total probability of infection from contact
          ProbWGInf   <- 1 - ( 1- probWGDCInf)*(1 - ProbWGCInf)*(1 - ProbBGCInf[order(unique(PopMatWB[, 8]), decreasing = FALSE)])
          tmpSus      <-which(PopMatWB[ ,9] == 0)
          
          if(length(tmpSus) > 0){
          NewInfPG <- tmpSus[runif(length(tmpSus)) < ProbWGInf[as.character(PopMatWB[tmpSus, 8])]]
            if(length(NewInfPG) > 0){
              PopMatWB[NewInfPG,  9] <- 1 ## exposed individuals
              PopMatWB[NewInfPG, 14] <- gTime + round(rpert(length(NewInfPG), 1, 5, 9))  # the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5.3,1.3)
              PopMatWB[NewInfPG, 15] <- PopMatWB[NewInfPG, 14] + round(rpert(length(NewInfPG), 1, 5, 7))# the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5,0.3)
            }
          }

# Update ------------------------------------------------------------------
            
            # Update animals to Infectious status
            IndexToInf <- PopMatWB[ ,9] == 1 & PopMatWB[ ,14] == gTime
            if(sum(IndexToInf) > 0) {
              PopMatWB[IndexToInf, 9] <- 2
              }
            
            # Update animals to Dead/Immune state
            IndexToDOrI <- which(PopMatWB[ ,9] == 2 & PopMatWB[ ,15] == gTime)
            if(length(IndexToDOrI) > 0){
              tmpToDie <- rbinom(length(IndexToDOrI), 1, ProbDeathInf)
              if(any(tmpToDie == 1)) {
                PopMatWB[IndexToDOrI[tmpToDie == 1], 9]  <- 3
                ## Dead animals may die within their own home range 20% or any of the connecting edges 80% (splitting male groups are excluded).
                tmpDieOSG   <- IndexToDOrI[tmpToDie == 1 & (!PopMatWB[IndexToDOrI, 2] %in% GroupsToSplit[ ,1])]
                DieOutGroup <- tmpDieOSG[rbinom(length(tmpDieOSG), 1, DieInConEdge)]
                if(length(DieOutGroup) > 0){
                  EdgesDie <- sapply(DieOutGroup, function(x) {
                    tmp1 <- as.numeric(WBMat[PopMatWB[x, 8], 6:13])
                    tmp2 <- WBMat[tmp1, 5]
                    tmp3 <- tmp1[tmp2 %in% 1:2]
                    tmp3 <- tmp3[tmp3 > 0]
                    if(length(tmp3) > 1) sample(tmp3, 1) else(PopMatWB[x, 8])
                  })
                  ### If an animal will die in a connecting edge, its pixel home and current will change but not the group number as the probability
                  ### of infection from Carcasses within a home-range and the connecting edges is the same according to Lange et al. (2015,EFSA) .
                  PopMatWB[DieOutGroup, 7] <- EdgesDie
                  PopMatWB[DieOutGroup, 8] <- EdgesDie
                }
              }
              if(any(tmpToDie == 0)) PopMatWB[IndexToDOrI[tmpToDie == 0], 9]  <- 4
            }
            
            
            ## Dead animals must be removed from the splitting mechanism
            # Identify group number of animals that died today
            IndexDiedTod  <- sort(unique(PopMatWB[PopMatWB[ ,9] == 3 & PopMatWB[ ,14] == gTime, 2]))
            MIndexDiedTod <- IndexDiedTod %in% GroupsToSplit[ ,1]
            if(sum(MIndexDiedTod) > 0){
              IndexToRemSM <- IndexDiedTod[MIndexDiedTod]
              NumRemTodFSM <- tapply((PopMatWB[PopMatWB[ ,2] %in% IndexToRemSM, 9] == 3 & 
                                      PopMatWB[PopMatWB[ ,2] %in% IndexToRemSM, 15] == gTime),
                                      PopMatWB[PopMatWB[ ,2] %in% IndexToRemSM, 2], sum)
              GroupsToSplit[match(IndexToRemSM, GroupsToSplit[ ,1]), 2] <- GroupsToSplit[match(IndexToRemSM, GroupsToSplit[ ,1]), 2] - NumRemTodFSM
              if(any(GroupsToSplit[ ,2] < 1)){
                IndexDelFSM   <- which(GroupsToSplit[ ,2] < 1)
                GroupsToSplit <- GroupsToSplit[-IndexDelFSM, , drop = FALSE]
              }
            }
            
            ## Remove dead carcasses over 28 days from death from the matrix
            indexDeadOld <- PopMatWB[ ,9] == 3 & ((gTime-PopMatWB[ ,15])>= CarcDaySur)
            if (any(indexDeadOld))  {
              PopMatWB<-PopMatWB[-which(indexDeadOld),]
            }
      } #closes ASF Dynamics
      
      
      ### REINITIATE & SUMMARISE
    
      ## Here re-initiate the daily matrix for males and females movements to make sure that there is no bleed from previous days
      MovedPixMale  <- as.list(matrix(0, ncol = 1))
      PixelsMoved   <- as.list(matrix(0, ncol = 1))
      PixelsMovedHS <- as.list(matrix(0, ncol = 1))
      
      ## Make the daily summary
      DayOutInfMat[gTime] <- length(unique(PopMatWB[PopMatWB[ ,9] %in% c(2:4), 2]))
      DayOutPopMat[gTime] <- length(unique(PopMatWB[PopMatWB[,9] %in% c(0:2,4) ,2]))
      DayOutAniMat[gTime] <- nrow(PopMatWB[PopMatWB[,9] %in% c(0:2,4) ,])
      DayOutCarcMat[gTime] <- length(PopMatWB[PopMatWB[,9] == 3 ,8])
      DayOutInfAniMat[gTime] <- length(PopMatWB[PopMatWB[,9] == 2 ,8])
      
        
      tmp <- tapply(PopMatWB[,9] == 2, PopMatWB[ ,8], sum)
      PopPerCells <- tapply(PopMatWB[,9] %in% c(0:2, 4), PopMatWB[ ,8], sum)
      tmp[tmp > 0] -> tmp
      PopPerCells <- PopPerCells[names(PopPerCells) %in% names(tmp)]
      nInfectedPerCells <- cbind(names(tmp), tmp, PopPerCells)

      InfGroupsD         <- unique(PopMatWB[PopMatWB[ ,9] %in% 2:4, 2])
      InfGroupsD         <- InfGroupsD[!InfGroupsD %in% NewInfGroups[NewInfGroups[ ,1] == Iter, 3]]
      InfcellsD         <-  unique(PopMatWB[PopMatWB[ ,2] %in% InfGroupsD, 8])

      
      InfAnimals         <- which(PopMatWB[  ,9] == 2)
      InfAnimals         <- InfAnimals[!InfAnimals %in% NewInfAnimals[NewInfAnimals[ ,1] == Iter, 3]]
      InfCarcass         <- which(PopMatWB[ ,9] == 3)
      InfCarcass         <- InfCarcass[!InfCarcass %in% NewInfCarcass[NewInfCarcass[ ,1] == Iter, 3]]
      
      
      if(length(InfGroupsD) > 0){
        NewInfGroups    <- rbind(NewInfGroups, cbind(Iter, gTime, InfGroupsD, InfcellsD))
      }
      if(length(InfAnimals) > 0){
        NewInfAnimals   <- rbind(NewInfAnimals, cbind(Iter, gTime, InfAnimals))
      }
      if(length(InfCarcass) > 0){
        NewInfCarcass <- rbind(NewInfCarcass, cbind(Iter, gTime, InfCarcass))
      }
      if(length(nInfectedPerCells) > 0){
        InfectedPerCells <- rbind(InfectedPerCells, cbind(Iter, gTime, nInfectedPerCells))
        }
      ## A variable to count whether the disease faded out this iteration and started up again 
      ## This is counted once/iteration despite in while loop
      TMPOutToday <- TMPOutYesterday & sum(PopMatWB[ ,9] %in% 1:2) > 0
      if(TMPOutToday & OnlyOnce){
        FreqRelapse[Iter] <- 1
        OnlyOnce <- FALSE
      }
      
      TMPOutYesterday <- gTime > TimeSeedInf & sum(PopMatWB[ ,9] %in% 1:2) == 0 & sum(PopMatWB[ ,9] == 3) > 0
     
     
  } #While(criteria...)
    
    to.res = list(NewInfAnimals    = NewInfAnimals,
                  NewInfCarcass    = NewInfCarcass,
                  NewInfGroups     = NewInfGroups,
                  DayOutInfMat     = DayOutInfMat,
                  DayOutPopMat     = DayOutPopMat,
                  DayOutAniMat     = DayOutAniMat,
                  DayOutCarcMat    = DayOutCarcMat,
                  DayOutInfAniMat  = DayOutInfAniMat,
                  InfectedPerCells = InfectedPerCells)
  } #Close for(Iter in...)
  
  parallel::stopCluster(cl1)
  save(res, file =  paste0(savingpath, runID, '.RData'))
}
 