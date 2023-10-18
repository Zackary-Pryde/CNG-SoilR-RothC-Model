# Required Packages/Dependencies TTT

library(pacman)
p_load(raster, rgdal, ncdf4, SoilR, abind, soilassessment, Formula, ggplot2, tidyverse)


JANSENVILLE_Weather_File = data.frame("Month" = 1:12,
                                      "Temp" = c(24.75, 24.75, 22.8, 19.7, 15.55, 12.75, 12.2, 14.45, 16.95, 19.15, 21.35, 23.3),
                                      "Precip" = c(25.83, 27.65, 44.31, 27.65, 11.77, 7.9, 13.69, 16.54, 13.69, 21.23, 29.46, 25.83),
                                      "Evp" = c(8.55, 7.97, 6.17, 4.76, 3.64, 2.98, 3.13, 4.06, 5.29, 6.46, 8.29, 8.87))

STRATUM_Edaphic_File = data.frame("Soil_Depth" = 30,
                                  "SOC_Stratum" = 97.4651,
                                  "ClayPerc_Stratum" = 21.4800)

ALMP_BL = data.frame("Month" = 1:12,
                     "Bare" = c(FALSE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
                     "Cinput" = rep(1.18131300047788, 12),
                     "FYM" = c(0,0,0,0,0,0,0,0,0,0,0,0),
                     "Irrigation" = c(0,0,0,0,0,0,0,0,0,0,0,0))

ALMP_PR = data.frame("Month" = 1:24,
                     "Bare" = c(FALSE,FALSE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,
                                FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
                     "Cinput" = c(rep(1.56466235908391, 12), 
                                  rep(1.96466235908391, 12)),
                     "FYM" = c(0,0,0,0,0,0,0,0,0,0,0,0,
                               0,0.5,0.5,0.5,0.5,0.5,0,0,0,0,0,0),
                     "Irrigation" = c(0,0,0,0,0,0,0,0,0,0,0,0,
                                      0,0,0,0,0,0,0,0,0,0,0,0))

calibrated_model_Input = data.frame("Soil_Carbon_Pool" = c("DPMptf", "RPMptf", "BIOptf", "HUMptf", "FallIOM"),
                                    "Value" = c(0.4510887, 13.3954103, 1.9193236, 72.2360421, 9.0259982))

Get_Delta_SOC_RothC = function(Weather_File, Edaphic_File, ALMP_File_BL, ALMP_File_PR, Calibrated_Model) {
  
  Years = seq(1/12,1,1/12) # Always 1 year
  num_years = ceiling(nrow(ALMP_File_PR)/12)
  month_ranges <- split(ALMP_File_PR, factor(rep(1:num_years, each = 12)))
  
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
      mutate(Bare = ifelse(Bare == TRUE, 1,0.6))
    
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
                             In = sum(ALMP_file$Cinput),
                             FYM = sum(ALMP_file$FYM),
                             DR = 1.44,
                             clay = edaphic_file$ClayPerc_Stratum,
                             xi = xi.frame, 
                             pass = TRUE, 
                             solver = euler)
    
    Ct3_spin=getC(Model3_spin)
    colnames(Ct3_spin)<-c("DPM", "RPM", "BIO", "HUM", "IOM")
    Model_Result = as.data.frame(Ct3_spin)
    Model_Result$SOC_Stock = Model_Result$DPM + Model_Result$RPM + Model_Result$BIO + Model_Result$HUM + Model_Result$IOM
    return(Model_Result)
  }
  
  Calibrated_Model_i_pr = Calibrated_Model
  Calibrated_Model_i_bl = Calibrated_Model
  
  SOC_MODEL_RESULT = list()
  
  for (year_index in names(month_ranges)) {
    
    cat("Project Year = ",year_index)
    
    ALMP_PR_Year_i = month_ranges[[year_index]]
    #year_number <- seq_along(names(month_ranges))
    
    SOC_BL = Roth_C(years = Years, weather_file = Weather_File, edaphic_file = Edaphic_File, ALMP_file = ALMP_File_BL, calibrated_model = Calibrated_Model_i_bl)
    names(SOC_BL) = paste0(names(SOC_BL), "_BL")
    
    cat("\n","BASELINE:", "\n","\n")
    print(SOC_BL)
    
    SOC_PR = Roth_C(years = Years, weather_file = Weather_File, edaphic_file = Edaphic_File, ALMP_file = ALMP_PR_Year_i, calibrated_model = Calibrated_Model_i_pr)
    names(SOC_PR) = paste0(names(SOC_PR), "_PR")
    
    cat("\n","PROJECT", "\n","\n")
    print(SOC_PR)
    
    SOC_MODEL_RESULT[[year_index]] = cbind(SOC_BL, SOC_PR)
    
    Calibrated_Model_i_pr = data.frame("Soil_Carbon_Pool" = c("DPMptf", "RPMptf", "BIOptf", "HUMptf", "FallIOM"),
                                       "Value" = c(SOC_PR[12,1:5]$DPM_PR, SOC_PR[12,1:5]$RPM_PR, SOC_PR[12,1:5]$BIO_PR, SOC_PR[12,1:5]$HUM_PR, SOC_PR[12,1:5]$IOM_PR))
    
    Calibrated_Model_i_bl = data.frame("Soil_Carbon_Pool" = c("DPMptf", "RPMptf", "BIOptf", "HUMptf", "FallIOM"),
                                       "Value" = c(SOC_BL[12,1:5]$DPM_BL, SOC_BL[12,1:5]$RPM_BL, SOC_BL[12,1:5]$BIO_BL, SOC_BL[12,1:5]$HUM_BL, SOC_BL[12,1:5]$IOM_BL))
  }
  
  return(SOC_MODEL_RESULT)
}

TEST = Get_Delta_SOC_RothC(Weather_File = JANSENVILLE_Weather_File,
                           Edaphic_File = STRATUM_Edaphic_File, 
                           ALMP_File_BL = ALMP_BL, 
                           ALMP_File_PR = ALMP_PR, 
                           Calibrated_Model = calibrated_model_Input)

RESULT = do.call(rbind,TEST)

RESULT$Month = 1:24

ggplot(data = RESULT, aes(x = Month)) + theme_minimal() + 
  geom_line(aes(y = SOC_Stock_BL), color = "darkblue", linetype = "dashed", lwd = 0.75) + 
  geom_line(aes(y = SOC_Stock_PR), color = "darkred", linetype = "dashed", lwd = 0.75) +
  labs(x = "Month",
       y = "SOC (t/Ha)",
       title = "Storms River Over Two Project Years") + 
  geom_vline(xintercept = c(1,13), linetype = "dashed", color = "darkgray")
