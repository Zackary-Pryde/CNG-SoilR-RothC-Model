# 0. File Information ----

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 
# This R script is the first attempt at RothC model CALIBRATION using SoilR.
# The code below follows the code set out at the following resource:
# https://fao-gsp.github.io/GSOCseq/software-environment.html
# 
# GOAL OF THE EXERCISE:
#   
#   Here, we want to run the MODEL SPIN UP using a specific value for Cinput.
#   
#   We then take the result, compare the SOC stock resulting from spin up and
#   compare it to the aerial average SOC stock obtained via a sampling campaign.
#   
#   If the difference between SOC_spinup and SOC_sampling is greater than 1%, 
#   then try another Cinput value for spin up. 
#   
#   Continue this process until the condition is met, at which point we save the
#   result as our calibrated model (i.e., initial SOC pool inputs).
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# 1. Required Packages/Dependencies ----

library(pacman)
p_load(raster, rgdal, ncdf4, SoilR, abind, soilassessment, Formula)

# 2. Defining the RothC function ----

# ROTH C MODEL FUNCTION .
Roth_C<-function(Cinputs,years,
                 DPMptf, RPMptf, BIOptf, HUMptf, FallIOM,
                 Temp,Precip,Evp,Cov2,soil.thick,SOC,clay,DR,bare1){
  #Temperature factor per month
  fT=fT.RothC(Temp[,2])
  
  #Moisture effects per month . 
  fw1func<-function(P, E, S.Thick = 30, pClay = 32.0213, pE = 1, bare) {
    M = P - E * pE
    Acc.TSMD = NULL
    for (i in 2:length(M)) {
      B = ifelse(bare[i] == FALSE, 1, 1.8)
      Max.TSMD = -(20 + 1.3 * pClay - 0.01 * (pClay^2)) * (S.Thick/23) * (1/B)
      Acc.TSMD[1] = ifelse(M[1] > 0, 0, M[1])
      if (Acc.TSMD[i - 1] + M[i] < 0) {
        Acc.TSMD[i] = Acc.TSMD[i - 1] + M[i]
      } else 
        (Acc.TSMD[i] = 0)
      if (Acc.TSMD[i] <= Max.TSMD) {
        Acc.TSMD[i] = Max.TSMD
      }
    }
    b = ifelse(Acc.TSMD > 0.444 * Max.TSMD, 1, (0.2 + 0.8 * ((Max.TSMD - Acc.TSMD)/(Max.TSMD - 0.444 * Max.TSMD))))
    b<-clamp(b,lower=0.2)
    return(data.frame(Acc.TSMD, b, Max.TSMD))
  }
  
  fW_2<- fw1func(P=(Precip[,2]), E=(Evp[,2]), S.Thick = soil.thick, pClay = clay, pE = 1, bare=bare1$Soil_Cover)$b
  
  #Vegetation Cover effects 
  fC <- Cov2[,2]
  
  # Set the factors frame for Model calculations
  xi.frame=data.frame(years,rep(fT*fW_2*fC*fPR,length.out=length(years)))
  
  # RUN THE MODEL FROM SOILR
  Model3_spin=RothCModel(t=years,C0=c(DPMptf, RPMptf, BIOptf, HUMptf, FallIOM),In=Cinputs,DR=DR,clay=clay,xi=xi.frame, pass=TRUE) 
  Ct3_spin=getC(Model3_spin)
  
  # Get the final pools of the time series
  poolSize3_spin=as.numeric(tail(Ct3_spin,1))
  return(poolSize3_spin)
}