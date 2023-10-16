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
Roth_C<-function(years,
                 weather_file,
                 edaphic_file,
                 ALMP_file, 
                 calibrated_model){
  
  fT.RothC_CNG = function (Temp) {
    47.91/(1 + exp(106.06/(ifelse(Temp >= -18.27, Temp, NA) + 18.27)))
  }
  
  fT = fT.RothC_CNG(weather_file[,2]) # Temperature effects per month
  
  fw1func<-function(P, E, S.Thick, pClay, pE = 1, bare) {
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
  
  fW_2<- fw1func(P=(weather_file[,3]), E=(weather_file[,4]), S.Thick = edaphic_file$Soil_Depth, pClay = edaphic_file$ClayPerc_Stratum, pE = 1, bare=ALMP_file$Bare)$b
  
  Cov2 = ALMP_file[,c("Month","Bare")] %>%
    mutate(Soil_Cover = ifelse(Bare == TRUE, 1,0.6))
  
  fC <- Cov2[,2]
  
  # Set the factors frame for Model calculations
  xi.frame=data.frame(years,rep(fT*fW_2*fC,length.out=length(years)))
  
  # RUN THE MODEL FROM SOILR
  Model3_spin = RothCModel(t = years,
                           C0 = c(calibrated_model[1,2], 
                                  calibrated_model[2,2], 
                                  calibrated_model[3,2], 
                                  calibrated_model[4,2], 
                                  calibrated_model[5,2]),
                           ks = c(k.DPM = 10, k.RPM = 0.3, k.BIO = 0.66, k.HUM = 0.02, k.IOM = 0),
                           In = ALMP_file$Cinput[1],
                           FYM = ALMP_file$FYM[1],
                           DR = 1.44,
                           clay = edaphic_file$ClayPerc_Stratum,
                           xi = xi.frame, 
                           pass = TRUE, 
                           solver = deSolve.lsoda.wrapper)
  
  Ct3_spin=getC(Model3_spin)
  
  # Get the final pools of the time series
  poolSize3_spin=as.numeric(tail(Ct3_spin,1))
  names(poolSize3_spin)<-c("DPM", "RPM", "BIO", "HUM", "IOM")
  
  SOC_Result = sum(poolSize3_spin)
  
  cat("\n","SIMULATION COMPLETE:", "\n","\n")
  
  print(poolSize3_spin)
  
  cat("\n","Modelled SOC after ",length(years)/12," years = ",SOC_Result, " (t/ha)")
  
  #print(paste0("Modelled SOC after ",length(years)/12," years = ",SOC_Result, " (t/ha)"))
  
  #return(poolSize3_spin)
}

# Defining the input variables BASELINE ----

JANSENVILLE_Weather_File = data.frame("Month" = 1:12,
                                      "Temp" = c(24.75, 24.75, 22.8, 19.7, 15.55, 12.75, 12.2, 14.45, 16.95, 19.15, 21.35, 23.3),
                                      "Precip" = c(25.83, 27.65, 44.31, 27.65, 11.77, 7.9, 13.69, 16.54, 13.69, 21.23, 29.46, 25.83),
                                      "Evp" = c(8.55, 7.97, 6.17, 4.76, 3.64, 2.98, 3.13, 4.06, 5.29, 6.46, 8.29, 8.87))

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

calibrated_model = data.frame("Soil_Carbon_Pool" = c("DPMptf", "RPMptf", "BIOptf", "HUMptf", "FallIOM"),
                              "Value" = c(0.4381499, 13.8652529, 1.9964292, 71.2632778, 9.0259982))

years = seq(1/12,2,by=1/12) 

Roth_C(years = years,
       weather_file = JANSENVILLE_Weather_File, 
       edaphic_file = STRATUM_Edaphic_File, 
       ALMP_file = ALMP_BL, 
       calibrated_model = calibrated_model)

Roth_C(years = years,
       weather_file = JANSENVILLE_Weather_File, 
       edaphic_file = STRATUM_Edaphic_File, 
       ALMP_file = ALMP_PR, 
       calibrated_model = calibrated_model)
