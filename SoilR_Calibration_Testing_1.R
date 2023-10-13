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
p_load(raster, rgdal, ncdf4, SoilR, abind, soilassessment, Formula, tidyverse)

# 2. Defining the RothC function ----

Roth_C<-function(Cinputs,years,
                 DPMptf, RPMptf, BIOptf, HUMptf, FallIOM,
                 Temp,Precip,Evp,soil.thick,SOC,clay,DR,bare1){
  
  fT.RothC_CNG = function (Temp_input) {
    47.91/(1 + exp(106.06/(ifelse(Temp_input >= -18.27, Temp_input, NA) + 18.27)))
  }
  
  fT = fT.RothC_CNG(Temp[,2])
  
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
  
  fW_2<- fw1func(P=(Precip[,2]), E=(Evp[,2]), S.Thick = soil.thick, pClay = clay, pE = 1, bare=bare1$Bare)$b
  
  Cov2 = bc %>%
    mutate(Soil_Cover = ifelse(Bare == TRUE, 1,0.6))
  
  fC <- Cov2[,2]
  
  xi.frame=data.frame(years,rep(fT*fW_2*fC,length.out=length(years)))
  
  Model3_spin=RothCModel(t=years,
                         C0=c(DPMptf, RPMptf, BIOptf, HUMptf, FallIOM),
                         ks = c(k.DPM = 10, k.RPM = 0.3, k.BIO = 0.66, k.HUM = 0.02, k.IOM = 0),
                         In=Cinputs,
                         DR=DR,
                         clay=clay,
                         xi=xi.frame, 
                         pass=TRUE, 
                         solver = euler) 
  Ct3_spin=getC(Model3_spin)
  
  poolSize3_spin=as.numeric(tail(Ct3_spin,1))
  return(poolSize3_spin)
}

# 3. Defining Calibration Input Conditions ----

# Climate Data

Temp = data.frame("Month" = 1:12, 
                  "Temp" = c(24.8, 24.8, 22.8, 19.7, 15.6, 12.8, 12.2, 14.4, 16.9, 19.1, 21.4, 23.3))

Precip = data.frame("Month" = 1:12, 
                    "Precip" = c(25.8, 27.6, 44.3, 27.6, 11.8, 7.9, 13.7, 16.5, 13.7, 21.2, 29.5, 25.8))

Evp = data.frame("Month" = 1:12, 
                 "Evp" = c(265.1, 223.1, 191.4, 142.8, 112.9, 89.3, 97.2, 125.9, 158.6, 200.3, 248.6, 274.9))

# Edaphic Data

Soil_Depth = 30
Clay_Stratum = 21.48
SOC_Stratum = 97.4651688163155

# RothC Simulation Duration

years = seq(1/12,1000,by=1/12)

# Bare/Cover Data

bc = data.frame("Month" = 1:12,
                "Bare" = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))

# Other SOC related calibration inputs

DPMi = 0
RPMi = 0
BIOi = 0
HUMi = 0
IOM = 0.049*SOC_Stratum^(1.139)

DPM_to_RPM_Ratio = 1.44
Cinputs_Initial = 1

Roth_C(Cinputs = Cinputs_Initial, 
       years = years, 
       DPMptf = DPMi,
       RPMptf = RPMi,
       BIOptf = BIOi,
       HUMptf = HUMi, 
       FallIOM = IOM, 
       Temp = Temp, Precip = Precip, Evp = Evp, 
       soil.thick = Soil_Depth, 
       SOC = SOC_Stratum, 
       clay = Clay_Stratum, 
       DR = DPM_to_RPM_Ratio, 
       bare1 = bc)


fb<-Roth_C(Cinputs=0,years=years,DPMptf=0, RPMptf=0, BIOptf=0, HUMptf=0, FallIOM=IOM,Temp=Temp,Precip=Precip,Evp=Evp,soil.thick=Soil_Depth,SOC=SOC_Stratum,clay=Clay_Stratum,DR=DPM_to_RPM_Ratio,bare1=bc)
fb_t<-fb[1]+fb[2]+fb[3]+fb[4]+fb[5]
m<-(fb_t-FallIOM)/(Cinputs_Initial)
Ceq<-(SOC_Stratum-FallIOM)/m
