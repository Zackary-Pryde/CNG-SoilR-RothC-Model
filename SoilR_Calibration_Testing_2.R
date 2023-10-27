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
p_load(raster, ncdf4, SoilR, abind, soilassessment, Formula, tidyverse)

# 2. Calibration Input Data (NB Changes to be made) ----

JANSENVILLE_Weather_File = data.frame("Month" = 1:12,
                                      "Temp" = c(24.75, 24.75, 22.8, 19.7, 15.55, 12.75, 12.2, 14.45, 16.95, 19.15, 21.35, 23.3),
                                      "Precip" = c(25.83, 27.65, 44.31, 27.65, 11.77, 7.9, 13.69, 16.54, 13.69, 21.23, 29.46, 25.83),
                                      "Evp" = c(8.55, 7.97, 6.17, 4.76, 3.64, 2.98, 3.13, 4.06, 5.29, 6.46, 8.29, 8.87))

STRATUM_Edaphic_File = data.frame("Soil_Depth" = 30,
                                  "SOC_Stratum" = 97.4651,
                                  "ClayPerc_Stratum" = 21.4800)

# 3. Defining the RothC CALIBRATION function ----

RothC_Calibration_CNG = function(Weather_File,Edaphic_File) {
  
  years = seq(1/12,1000,by=1/12) 
  
  fT.RothC_CNG = function (Temp) {
    47.91/(1 + exp(106.06/(ifelse(Temp >= -18.27, Temp, NA) + 18.27)))
  }
  
  fT = fT.RothC_CNG(Weather_File[,2]) # Temperature effects per month
  
  fW = fW.RothC(P=(Weather_File[,3]), E=(Weather_File[,4]), 
                S.Thick = Edaphic_File$Soil_Depth, 
                pClay = Edaphic_File$ClayPerc_Stratum, 
                pE = 1, bare = FALSE)$b
  
  xi.frame = data.frame(years,
                        rep(fT*fW, length.out = length(years)))
  
  FallIOM=0.049*Edaphic_File$SOC_Stratum^(1.139)
  
  Cinputs = 1
  Perc_Diff = 1
  
  while (Perc_Diff > 0.01) {
    
    Model1=RothCModel(t = years,
                      C0=c(DPM = 0, RPM = 0, BIO = 0, HUM = 0, IOM = FallIOM),
                      In = Cinputs, 
                      clay = Edaphic_File$ClayPerc_Stratum, 
                      DR = 1.44,
                      xi = xi.frame, pass = TRUE, solver = deSolve.lsoda.wrapper) # Loads the model
    
    Ct1 = getC(Model1) # Calculates stocks for each pool per month
    
    poolSize1=as.numeric(tail(Ct1,1))
    names(poolSize1)<-c("DPM", "RPM", "BIO", "HUM", "IOM")
    
    SOC_Cali = sum(poolSize1)
    
    Perc_Diff = (Edaphic_File$SOC_Stratum - SOC_Cali)/Edaphic_File$SOC_Stratum
    
    print(paste0("Difference = ",Edaphic_File$SOC_Stratum - SOC_Cali," (",round(Perc_Diff*100, 2)," %)"))
    
    Cinputs = Cinputs + Perc_Diff*Cinputs
    
    if (Perc_Diff <= 0.01) {
      cat("\n","CALIBRATION COMPLETE:", "\n","\n")
    }
  }
  return(poolSize1)
}

# 4. Usage ----

RothC_Calibration_CNG(Weather_File = JANSENVILLE_Weather_File,
                      Edaphic_File = STRATUM_Edaphic_File)
