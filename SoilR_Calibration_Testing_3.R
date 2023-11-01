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
#   Additions to SoilR_Calibration_Testing_2 Workflow:
#     - Extending Calibration function to use a specified set of calibrated models
#       which may or may not exist in the FFM table. 
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

RothC_Calibration_CNG_Many_Models = function(Models_To_Calibrate) {
  
  for (Model_Input in levels(as.factor(Models_To_Calibrate))) {
    
    # Get the Model that has been specified
    Model_i <- Model_Input
    
    # Get Stratum and weather station from model name... NB - Note that in this case model name MUST be in the format of: "Stratum - Weather_Station"
    Stratum_i = unlist(strsplit(Model_Input, " - "))[1]
    Weather_Station_i = unlist(strsplit(Model_Input, " - "))[2]
    
    # Get the Info for specified STRATUM
    Edaphic_File <- dbGetQuery(connection, paste0("SELECT * FROM dbo.SoilR_Stratum_File WHERE Stratum = '", Stratum_i,"'")) %>%
      dplyr::select(Soil_Depth, "SOC_Stratum" = SOC_Stock, "ClayPerc_Stratum" = Clay_Percentage)
    
    # Get the Info for specified WEATHER STATION
    Weather_File = dbGetQuery(connection, paste0("SELECT * FROM dbo.SoilR_Weather_File WHERE Weather_Station = '", Weather_Station_i,"'")) %>%
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
    
    Calibrated_Model = data.frame("Model_Name" = Model_i,
                                  "Stratum" = Stratum_i,
                                  "Weather_Station" = Weather_Station_i,
                                  "DPM" = poolSize1[1],
                                  "RPM" = poolSize1[2],
                                  "BIO" = poolSize1[3],
                                  "HUM" = poolSize1[4],
                                  "IOM" = poolSize1[5])
    
    DF = Calibrated_Model
    DB <- dbGetQuery(connection, "SELECT * FROM dbo.SoilR_Calibrated_Model")
    
    if (any(DB$Model_Name == DF$Model_Name)) {
      DB[which(DB$Model_Name == DF$Model_Name),] = DF[1,]
    } else {
      DB = rbind(DB, DF)
    }
    
    dbWriteTable(connection, name = "SoilR_Calibrated_Model", DB, overwrite = TRUE)
    print(paste0("Calibration Complete: ", Model_i))
  }
  dbDisconnect(connection)
}

# 4. Usage ----

Models = c("Lowland fynbos and renosterveld / Cfb : Temperate, no dry season, warm summer - PORT-ELIZABETH",
           "Drakensberg montane grasslands, woodlands and forests / Cwb : Temperate, dry winter, warm summer - KOKSTAD",
           "Drakensberg montane grasslands, woodlands and forests / Cwb : Temperate, dry winter, warm summer - ESTCOURT",
           "Drakensberg montane grasslands, woodlands and forests / Cfa : Temperate, no dry season, hot summer - PORT-ELIZABETH",
           "Drakensberg montane grasslands, woodlands and forests / BSk : Arid, steppe, cold - PORT-ELIZABETH",
           "Knysna-Amatole montane forests / Cfb : Temperate, no dry season, warm summer - JANSENVILLE",
           "Knysna-Amatole montane forests / Cfa : Temperate, no dry season, hot summer - JANSENVILLE",
           "Lowland fynbos and renosterveld / Cfb : Temperate, no dry season, warm summer - JANSENVILLE",
           "KwaZulu-Cape coastal forest mosaic / BSh : Arid, steppe, hot - PORT-ELIZABETH",
           "Knysna-Amatole montane forests / Cfb : Temperate, no dry season, warm summer - GEORGE-AIRPORT",
           "Drakensberg montane grasslands, woodlands and forests / Cfb : Temperate, no dry season, warm summer - QUEENSTOWN",
           "Nama Karoo / BSk : Arid, steppe, cold - PORT-ELIZABETH")

Testing = RothC_Calibration_CNG_Many_Models(Models_To_Calibrate = Models)
