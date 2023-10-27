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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#   
#   Additions to SoilR_Calibration_Testing_1 Workflow:
#     - Extending Calibration function to use SQL databases
#     - User input of Model name (run for loop on model names if need be)
#     - Sink to SQL table
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# 1. Required Packages/Dependencies ----

library(pacman)
p_load(raster, ncdf4, SoilR, abind, soilassessment, Formula, ggplot2, tidyverse, odbc, tidyverse, RODBC, DBI, dplyr, keyring)

# SQL Database connection
keyring_unlock("cng_SQL_Credentials")
connection <- odbc::dbConnect(odbc(),Driver = "ODBC Driver 17 for SQL Server",
                              Server = key_get("SQL_Server", keyring ="cng_SQL_Credentials"),
                              Database = key_get("Database", keyring ="cng_SQL_Credentials"),
                              UID = key_get("Username", keyring ="cng_SQL_Credentials"),
                              PWD = key_get("Password", keyring ="cng_SQL_Credentials"))
keyring_lock("cng_SQL_Credentials")

# 3. Defining the RothC CALIBRATION function ----

RothC_Calibration_CNG = function(Model_Name_Input) {
  
  # Obtain Data Associated with Input
  
  # Get the Model that has been specified
  FFM_Information <- dbGetQuery(connection, paste0("SELECT DISTINCT Calibrated_Model, Stratum, Weather_Station FROM dbo.SoilR_Farm_Field_Master WHERE Calibrated_Model = '", Model_Name_Input,"'"))
  
  # Get the Info for specified STRATUM
  Edaphic_File <- dbGetQuery(connection, paste0("SELECT * FROM dbo.SoilR_Stratum_File WHERE Stratum = '", FFM_Information$Stratum,"'")) %>%
    dplyr::select(Soil_Depth, "SOC_Stratum" = SOC_Stock, "ClayPerc_Stratum" = Clay_Percentage)
  
  # Get the Info for specified WEATHER STATION
  Weather_File = dbGetQuery(connection, paste0("SELECT * FROM dbo.SoilR_Weather_File WHERE Weather_Station = '", FFM_Information$Weather_Station,"'")) %>%
    dplyr::select(Month, "Temp" = Temperature, "Precip" = Precipitation, "Evp" = Evapotranspiration)
  
  # Calibration on Input data
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
    
    SOC_Cali = sum(poolSize1)
    
    Perc_Diff = (Edaphic_File$SOC_Stratum - SOC_Cali)/Edaphic_File$SOC_Stratum
    
    Cinputs = Cinputs + Perc_Diff*Cinputs
  }
  
  Calibrated_Model = data.frame("Model_Name" = FFM_Information$Calibrated_Model,
                                "Stratum" = FFM_Information$Stratum,
                                "Weather_Station" = FFM_Information$Weather_Station,
                                "DPM" = poolSize1[1],
                                "RPM" = poolSize1[2],
                                "BIO" = poolSize1[3],
                                "HUM" = poolSize1[4],
                                "IOM" = poolSize1[5])
  
  DF = Calibrated_Model
  DB <- dbGetQuery(connection, "SELECT * FROM dbo.SoilR_Calibrated_Model")
  
  if (any(DB$Model_Name == Calibrated_Model$Model_Name)) {
    DB[which(DB$Model_Name == Calibrated_Model$Model_Name),] = Calibrated_Model[1,]
  } else {
    DB = rbind(DB, DF)
  }
  
  dbWriteTable(connection, name = "SoilR_Calibrated_Model", DB, overwrite = TRUE)
  
  dbDisconnect(connection)
  return(Calibrated_Model)
}

# 4. Usage ----

Testing = RothC_Calibration_CNG(Model_Name_Input = "Test - JANSENVILLE")


