rpert <- function(n, a, l, b) {
  mu <- (a+4*l+b)/6
  if (mu==l) v <- w <- 3 else {
    v <- (mu-a)*(2*l-a-b)/(l-mu)/(b-a)
    w <- v*(b-mu)/(mu-a)
  }
  a+(b-a)*rbeta(n, v, w)
}

# Definition of WB Matrix --------------------------------------------------

  DefineWbMat <- function(file, MinBCap, ModeBCap, MaxBCap){
  WBMat <- read.table(file, sep = ",", header = T)
  WBMat <- cbind(WBMat, 0)
  WBMat <- WBMat[ ,-1]
  colnames(WBMat) <- c("ID", "Lon", "Lat", "Forest_cover", "Habitat", "Nb1", "Nb2",
                       "Nb3", "Nb4", "Nb5", "Nb6", "Nb7", "Nb8", "Colors", "Breeding")
  WBMat[WBMat[ ,5] == 1, 15] <- ceiling(rpert(sum(WBMat[ ,5] == 1), MinBCap, ModeBCap, MaxBCap))
  return(WBMat)  
}


 InitPop <- function(InitialOccCells, 
                    Mat     = WBMat,
                    Fem     = FemPropGroup,
                    Mal     = MalPropGroup,
                    SubA    = SubAPropGroup,
                    Piglet  = PigletPropGroup,
                    AgeProb = AgeProbab, MinBCap, ModeBCap, MaxBCap){
  TMPHP <- Mat[ ,5] == 1
  homePixelsAll <- sample(Mat[TMPHP, 1], InitialOccCells)
  PopMatWB = NULL
  size = 0
  for(i in 1:InitialOccCells){
    ## The initial distribution of the groups was based on Merta et al. (2015). They were calculated relative to the average
    ## Female percentage. for instance Males were average % Males divided by average % Females. 
    Females       <- round(runif(MinBCap, ModeBCap, MaxBCap))
    Males         <- round(Females*Mal/Fem)
    SubAdults     <- round(Females*SubA/Fem)
    Piglets       <- round(Females*Piglet/Fem)
    GroupID       <- rep(i, sum(Females, Males, SubAdults, Piglets))
    Sex           <- c(rep(1, Females), rep(0, Males), rbinom(SubAdults, 1, prob = 0.5), rbinom(Piglets, 1, prob = 0.5))
    AgeCat        <- c(rep(3, Females),rep(3, Males), rep(2, SubAdults), rep(1, Piglets))
    Dam           <- c(rep(0, sum(Females, Males)), sample(1:Females, SubAdults, rep = T), sample(1:Females, Piglets, rep = T))
    # Females deliver between Jan and June, so the age of the Piglets will be between 183 and 365 days.
    # The same thing for sub-adults but with 365 days extra.
    tmpAgeSubA    <- sample((183:365) + 365, Females, rep = T)
    tmpAgePig     <- sample((183:365), Females, rep = T)
    Age           <- c(sample(2:10, size = sum(Females, Males), T, prob = AgeProb)*365, tmpAgeSubA[Dam[(Females + Males + 1):(Females + Males + SubAdults)]],
                       tmpAgePig[Dam[(Females + Males + SubAdults + 1):(Females + Males + SubAdults + Piglets)]])
    # Adjust the Dam number to fit the actual ID of the DAMs 
    Dam[(Females + Males + 1):(Females + Males + SubAdults)] <- Dam[(Females + Males + 1):(Females + Males + SubAdults)] + size# -1 because we have extra row at start
    Dam[(Females + Males + SubAdults + 1):(Females + Males + SubAdults + Piglets)] <- Dam[(Females + Males + SubAdults + 1):(Females + Males + SubAdults + Piglets)] + size
    Breed         <- rep(0, sum(Females, Males, SubAdults, Piglets))
    HomePixel     <- rep(homePixelsAll[i], sum(Females, Males, SubAdults, Piglets))
    CurrPixel     <- HomePixel
    infectStatus  <- rep(0, sum(Females, Males, SubAdults, Piglets))
    SplitStatus   <- rep(0, sum(Females, Males, SubAdults, Piglets))
    SplitMale     <- rep(0, sum(Females, Males, SubAdults, Piglets))
    TimeDeath     <- rep(0, sum(Females, Males, SubAdults, Piglets))
    TimeToInfect  <- rep(0, sum(Females, Males, SubAdults, Piglets))
    TimeToDeath   <- rep(0, sum(Females, Males, SubAdults, Piglets))
    GroupMove     <- rep(0, sum(Females, Males, SubAdults, Piglets))
    IDs           <- (size+ 1):(size + sum(Females, Males, SubAdults, Piglets))
    
    InitMatWBPop  <- cbind(IDs, GroupID, Sex, AgeCat, Age, Breed, HomePixel, CurrPixel, infectStatus, 
                           SplitStatus, Dam, SplitMale, TimeDeath, TimeToInfect, TimeToDeath, GroupMove)
    PopMatWB   <- rbind(PopMatWB, InitMatWBPop)
    # Mat[unique(HomePixel), 15] <- Females
    size = nrow(PopMatWB)
    
  }
  
  colnames(PopMatWB) <- c("IDs", "Group_ID", "Sex", "Age_Cat", "Age_days", "Breed",
                          "Home_pixel", "Current_pixel", "Infect_status", "Split_status", 
                          "Dam", "Split_male", "Time_death","Time_to_infect", "Time_to_death", "Group_move")
  return(PopMatWB)
}
