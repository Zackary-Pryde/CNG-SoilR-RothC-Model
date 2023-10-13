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

Temp = data.frame("Month" = 1:12, 
                  "Temp" = c(24.75, 24.75, 22.8, 19.7, 15.55, 12.75, 12.2, 14.45, 16.95, 19.15, 21.35, 23.3))

Precip = data.frame("Month" = 1:12, 
                    "Precip" = c(25.83, 27.65, 44.31, 27.65, 11.77, 7.9, 13.69, 16.54, 13.69, 21.23, 29.46, 25.83)) # NOTE that these values are excluding any irrigation

Evp=data.frame(Month=1:12, Evp=c(12, 18, 35, 58, 82, 90, 97, 84, 54, 31,
                                 14, 10))

soil.thick = 30 # Soil thickness (organic layer topsoil) (cm)
SOC = 97.4651 # Soil Organic Carbon Stock (Mg/ha). NB: Mg refers to Megagram = metric Tonne.
clay = 21.4800 # Percent clay (%)
Cinputs = 0.5 # Annual C inputs to soil (Mg/ha/yr). NB: This is the value that we need to check on. For now, Ill assume its the same as that used to calibrate

years = seq(1/12,1000,by=1/12) 

fT = fT.RothC(Temp[,2]) # Temperature effects per month

fW = fW.RothC(P=(Precip[,2]), E=(Evp[,2]), 
              S.Thick = soil.thick, 
              pClay = clay, 
              pE = 1, bare = FALSE)$b # Moisture effects per month

xi.frame = data.frame(years,
                      rep(fT*fW, length.out = length(years)))

FallIOM=0.049*SOC^(1.139) #IOM using Falloon method


# 2. Defining the RothC CALIBRATION function ----


Cinputs = 1
Perc_Diff = 1

while (Perc_Diff > 0.01) {
  
  Model1=RothCModel(t = years,
                    C0=c(DPM = 0, RPM = 0, BIO = 0, HUM = 0, IOM = FallIOM),
                    In = Cinputs, 
                    clay = clay, 
                    xi = xi.frame, pass = TRUE, solver = deSolve.lsoda.wrapper) # Loads the model
  
  Ct1 = getC(Model1) # Calculates stocks for each pool per month
  
  poolSize1=as.numeric(tail(Ct1,1))
  
  SOC_Cali = sum(poolSize1)
  
  Perc_Diff = (SOC - SOC_Cali)/SOC
  
  print(paste0("Difference = ",SOC - SOC_Cali," (",round(Perc_Diff*100, 2)," %)"))
  
  Cinputs = Cinputs + Perc_Diff*Cinputs}


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
