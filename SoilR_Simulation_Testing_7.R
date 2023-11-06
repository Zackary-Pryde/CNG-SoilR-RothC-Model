library(pacman)
p_load(raster, ncdf4, SoilR, abind, soilassessment, Formula, ggplot2, tidyverse, odbc, tidyverse, RODBC, DBI, dplyr, keyring)

keyring_unlock("cng_SQL_Credentials")
connection <- odbc::dbConnect(odbc(),Driver = "ODBC Driver 17 for SQL Server",
                              Server = key_get("SQL_Server", keyring ="cng_SQL_Credentials"),
                              Database = key_get("Database", keyring ="cng_SQL_Credentials"),
                              UID = key_get("Username", keyring ="cng_SQL_Credentials"),
                              PWD = key_get("Password", keyring ="cng_SQL_Credentials"))
keyring_lock("cng_SQL_Credentials")

# Define Simulation Function
Get_Delta_SOC_RothC = function(Paddock_UID_Input) {
  
  FFM_Information <- dbGetQuery(connection, paste0("SELECT * FROM dbo.SoilR_Farm_Field_Master WHERE Paddock_UID = '", Paddock_UID_Input,"'"))
  
  Weather_File = dbGetQuery(connection, paste0("SELECT * FROM dbo.SoilR_Weather_File WHERE Weather_Station = '", FFM_Information$Weather_Station,"'")) %>%
    dplyr::select(Month, "Temp" = Temperature, "Precip" = Precipitation, "Evp" = Evapotranspiration)
  
  Edaphic_File = dbGetQuery(connection, paste0("SELECT * FROM dbo.SoilR_Stratum_File WHERE Stratum = '", FFM_Information$Stratum,"'")) %>%
    dplyr::select(Soil_Depth, "SOC_Stratum" = SOC_Stock, "ClayPerc_Stratum" = Clay_Percentage)
  
  ALMP_File_BL = dbGetQuery(connection, paste0("SELECT * FROM dbo.SoilR_ALM_File WHERE Paddock_UID = '", FFM_Information$Paddock_UID,"'", "AND Scenario = 'Baseline'")) %>%
    dplyr::select(Month, Bare, Cinput, FYM, Irrigation)
  
  ALMP_File_PR = dbGetQuery(connection, paste0("SELECT * FROM dbo.SoilR_ALM_File WHERE Paddock_UID = '", FFM_Information$Paddock_UID,"'", "AND Scenario = 'Project'")) %>%
    dplyr::select(Month, Bare, Cinput, FYM, Irrigation)
  
  Calibrated_Model = dbGetQuery(connection, paste0("SELECT * FROM dbo.SoilR_Calibrated_Model WHERE Model_Name = '", FFM_Information$Calibrated_Model,"'")) %>%
    dplyr::select(DPM,RPM,BIO,HUM,IOM)
  
  Calibrated_Model = data.frame("Soil_Carbon_Pool" = c("DPMptf", "RPMptf", "BIOptf", "HUMptf", "FallIOM"),
                                "Value" = c(Calibrated_Model$DPM, Calibrated_Model$RPM, Calibrated_Model$BIO, Calibrated_Model$HUM, Calibrated_Model$IOM))
  
  Years = seq(1/12,2,1/12) # Always 1 year
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
    
    fW_2<- fw1func(P=(weather_file[,3] + ALMP_file[,5]), E=(weather_file[,4]), S.Thick = edaphic_file$Soil_Depth, pClay = edaphic_file$ClayPerc_Stratum, pE = 1, bare=ALMP_file$Bare)$b
    
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
    
    ALMP_PR_Year_i = month_ranges[[year_index]]
    
    SOC_BL = Roth_C(years = Years, weather_file = Weather_File, edaphic_file = Edaphic_File, ALMP_file = ALMP_File_BL, calibrated_model = Calibrated_Model_i_bl)
    names(SOC_BL) = paste0(names(SOC_BL), "_BL")
    
    SOC_PR = Roth_C(years = Years, weather_file = Weather_File, edaphic_file = Edaphic_File, ALMP_file = ALMP_PR_Year_i, calibrated_model = Calibrated_Model_i_pr)
    names(SOC_PR) = paste0(names(SOC_PR), "_PR")
    
    Calibrated_Model_i_pr = data.frame("Soil_Carbon_Pool" = c("DPMptf", "RPMptf", "BIOptf", "HUMptf", "FallIOM"),
                                       "Value" = c(SOC_PR[13,1:5]$DPM_PR, SOC_PR[13,1:5]$RPM_PR, SOC_PR[13,1:5]$BIO_PR, SOC_PR[13,1:5]$HUM_PR, SOC_PR[13,1:5]$IOM_PR))
    
    Calibrated_Model_i_bl = data.frame("Soil_Carbon_Pool" = c("DPMptf", "RPMptf", "BIOptf", "HUMptf", "FallIOM"),
                                       "Value" = c(SOC_BL[13,1:5]$DPM_BL, SOC_BL[13,1:5]$RPM_BL, SOC_BL[13,1:5]$BIO_BL, SOC_BL[13,1:5]$HUM_BL, SOC_BL[13,1:5]$IOM_BL))
    
    SOC_BL = SOC_BL[1:12,]
    SOC_PR = SOC_PR[1:12,]
    
    SOC_MODEL_RESULT[[year_index]] = cbind(SOC_BL, SOC_PR)
  }
  
  SOC_MODEL_RESULT_DF = do.call(rbind,SOC_MODEL_RESULT)
  SOC_MODEL_RESULT_DF$Month = 1:nrow(SOC_MODEL_RESULT_DF)
  SOC_MODEL_RESULT_DF$Delta_SOC_Stock = SOC_MODEL_RESULT_DF$SOC_Stock_PR - SOC_MODEL_RESULT_DF$SOC_Stock_BL
  SOC_MODEL_RESULT_DF$Paddock_UID = FFM_Information$Paddock_UID
  
  SOC_MODEL_RESULT_DF = SOC_MODEL_RESULT_DF[,c("Paddock_UID", "Month", "DPM_BL", "RPM_BL", "BIO_BL", "HUM_BL", "IOM_BL", "SOC_Stock_BL", "DPM_PR", "RPM_PR", "BIO_PR", "HUM_PR", "IOM_PR", "SOC_Stock_PR", "Delta_SOC_Stock")]
  rownames(SOC_MODEL_RESULT_DF) = NULL
  
  DF = SOC_MODEL_RESULT_DF
  DB <- dbGetQuery(connection, "SELECT * FROM dbo.SoilR_Simulation_Raw_Output")
  
  if (any(DB$Paddock_UID == levels(as.factor(DF$Paddock_UID)))) {
    DB = DB[-which(DB$Paddock_UID == levels(as.factor(DF$Paddock_UID))),]
    DB = rbind(DB, DF)
    rownames(DB) = NULL
  } else {
    DB = rbind(DB, DF)
    rownames(DB) = NULL
  }
  
  dbWriteTable(connection, name = "SoilR_Simulation_Raw_Output", DB, overwrite = TRUE)
  
  Years_to_Report = seq(12,nrow(DF),12)
  DF = DF[Years_to_Report,c("Paddock_UID","SOC_Stock_BL","SOC_Stock_PR","Delta_SOC_Stock")]
  DF$Year = 1:nrow(DF)
  
  DF <- DF %>%
    mutate(Delta_SOC_Stock_Claimable = ifelse(Year == 1, Delta_SOC_Stock, Delta_SOC_Stock - lag(Delta_SOC_Stock)))
  
  DF$Field_Size = FFM_Information$Field_Size # NB Need to change code to account for year-specific field size
  DF$Emissions_Output = (DF$Delta_SOC_Stock_Claimable * DF$Field_Size) * (44/12)
  
  DF = DF[,c("Paddock_UID", "Year", "Field_Size", "SOC_Stock_BL", "SOC_Stock_PR","Delta_SOC_Stock", "Delta_SOC_Stock_Claimable", "Emissions_Output")]
  
  DB <- dbGetQuery(connection, "SELECT * FROM dbo.SoilR_Simulation_Summary_Output")
  
  if (any(DB$Paddock_UID == levels(as.factor(DF$Paddock_UID)))) {
    DB = DB[-which(DB$Paddock_UID == levels(as.factor(DF$Paddock_UID))),]
    DB = rbind(DB, DF)
    rownames(DB) = NULL
  } else {
    DB = rbind(DB, DF)
    rownames(DB) = NULL
  }
  
  dbWriteTable(connection, name = "SoilR_Simulation_Summary_Output", DB, overwrite = TRUE)
  
}

# Define Iteration Function (NB Takes a vector Input)
Get_Delta_SOC_RothC_For_Many = function(Input_Paddock_UIDs) {
  
  for (Paddock in levels(as.factor(Input_Paddock_UIDs))) {
    Paddock_UID_Input = Paddock
    Get_Delta_SOC_RothC(Paddock_UID_Input = Paddock)
    print(paste0((which(levels(as.factor(Input_Paddock_UIDs)) == Paddock)/length(levels(as.factor(Input_Paddock_UIDs))))*100, " %"))
  }
  
  dbDisconnect(connection)
}

# USAGE

# Obtain the Paddock_UIDs to Run

#Paddock_UIDs_To_Run <- dbGetQuery(connection, paste0("SELECT DISTINCT [Paddock_UID] FROM [dbo].[SoilR_Farm_Field_Master] WHERE Paddock_UID != 'TestF1'"))
Paddock_UIDs_To_Run <- dbGetQuery(connection, paste0("SELECT DISTINCT [Paddock_UID] FROM [dbo].[SoilR_Farm_Field_Master]"))
Paddock_UIDs_To_Run = Paddock_UIDs_To_Run$Paddock_UID

# Run the new function with the Paddock_UID vector as its input
Get_Delta_SOC_RothC_For_Many(Input_Paddock_UIDs = Paddock_UIDs_To_Run)
