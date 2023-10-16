# 0. File Information ----

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 
# This R script is the first attempt at RothC model SIMULATION using SoilR.
# The code below follows the code set out at the following resource:
# https://fao-gsp.github.io/GSOCseq/software-environment.html
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# 1. Required Packages/Dependencies ----

library(pacman)
p_load(raster, rgdal, ncdf4, SoilR, abind, soilassessment, Formula)

# GSOC METHOD ----

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# NB I am setting up these data to track that used for the PD1 MR1 Farm: Storms River.
#   Weather File = JANSENVILLE
#   Stratum = Knysna-Amatole montane forests / Cfb : Temperate, no dry season, warm summer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# ROTH C MODEL FUNCTION .
Roth_C<-function(Cinputs,years,
                 DPMptf, RPMptf, BIOptf, HUMptf, FallIOM,
                 Temp,Precip,Evp,soil.thick,clay,DR,bare1){
  
  fT.RothC_CNG = function (Temp) {
    47.91/(1 + exp(106.06/(ifelse(Temp >= -18.27, Temp, NA) + 18.27)))
  }
  
  fT = fT.RothC_CNG(Weather_File[,2]) # Temperature effects per month
  
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
  
  Cov2 = bc %>%
    mutate(Soil_Cover = ifelse(Soil_Cover == TRUE, 1,0.6))
  
  fC <- Cov2[,2]
  
  # Set the factors frame for Model calculations
  xi.frame=data.frame(years,rep(fT*fW_2*fC,length.out=length(years)))
  
  # RUN THE MODEL FROM SOILR
  Model3_spin=RothCModel(t=years,C0=c(DPMptf, RPMptf, BIOptf, HUMptf, FallIOM),In=Cinputs,DR=1.44,clay=clay,xi=xi.frame, pass=TRUE, solver = deSolve.lsoda.wrapper) 
  Ct3_spin=getC(Model3_spin)
  
  # Get the final pools of the time series
  poolSize3_spin=as.numeric(tail(Ct3_spin,1))
  return(poolSize3_spin)
}

# Defining the input variables BASELINE ----

JANSENVILLE_Weather_File = data.frame("Month" = 1:12,
                                      "Temp" = c(24.75, 24.75, 22.8, 19.7, 15.55, 12.75, 12.2, 14.45, 16.95, 19.15, 21.35, 23.3),
                                      "Precip" = c(25.83, 27.65, 44.31, 27.65, 11.77, 7.9, 13.69, 16.54, 13.69, 21.23, 29.46, 25.83),
                                      "Evp" = c(12, 18, 35, 58, 82, 90, 97, 84, 54, 31, 14, 10))

STRATUM_Edaphic_File = data.frame("Soil_Depth" = 30,
                                  "SOC_Stratum" = 97.4651,
                                  "ClayPerc_Stratum" = 21.4800)

ALMP_BL = data.frame("Month" = 1:12,
                     "Bare" = c(FALSE,FALSE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
                     "Cinput" = rep(1.18131300047788, 12),
                     "FYM" = c(0,0,0,0,0,0,0,0,0,0,0,0),
                     "Irrigation" = c(0,0,0,0,0,0,0,0,0,0,0,0))

ALMP_PR = data.frame("Month" = 1:12,
                     "Bare" = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
                     "Cinput" = rep(1.56466235908391, 12),
                     "FYM" = c(0,0,0,0,0,0,0,0,0,0,0,0),
                     "Irrigation" = c(0,0,0,0,0,0,0,0,0,0,0,0))

# Initial conditions from spin up
DPMptf = 0.4371314
RPMptf = 13.8658136
BIOptf = 1.9966239
HUMptf = 71.2542223
FallIOM = 9.0259982

years = seq(1/12,2,by=1/12) 

#SOC_BL = 

