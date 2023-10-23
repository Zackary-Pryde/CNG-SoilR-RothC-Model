# 0. File Information ----

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 
# The purpose of this R script is to INSERT SoilR data into SQL tables in our database. 
# These data are used in:
#   a) RothC Model Calibration
#   b) RothC Model Simulation
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# 1. Required Packages/Dependencies ----

library(pacman)
p_load(odbc,tidyverse,RODBC,DBI,dplyr, keyring)

# 2. Connecting to the SQL database securely ----

# KEYRING METHOD

keyring_unlock("cng_SQL_Credentials")
connection <- odbc::dbConnect(odbc(),Driver = "ODBC Driver 17 for SQL Server",
                              Server = key_get("SQL_Server", keyring ="cng_SQL_Credentials"),
                              Database = key_get("Database", keyring ="cng_SQL_Credentials"),
                              UID = key_get("Username", keyring ="cng_SQL_Credentials"),
                              PWD = key_get("Password", keyring ="cng_SQL_Credentials"))
keyring_lock("cng_SQL_Credentials")

# Weather File to be inserted - - - - - - - - - - - - - -

JANSENVILLE_Weather_File = data.frame("Month" = 1:12,
                                      "Temperature" = c(24.75, 24.75, 22.8, 19.7, 15.55, 12.75, 12.2, 14.45, 16.95, 19.15, 21.35, 23.3),
                                      "Precipitation" = c(25.83, 27.65, 44.31, 27.65, 11.77, 7.9, 13.69, 16.54, 13.69, 21.23, 29.46, 25.83),
                                      "Evapotranspiration" = c(8.55, 7.97, 6.17, 4.76, 3.64, 2.98, 3.13, 4.06, 5.29, 6.46, 8.29, 8.87),
                                      "Weather_Station" = rep("JANSENVILLE", 12))

DF = JANSENVILLE_Weather_File
DB <- dbGetQuery(connection, "SELECT * FROM dbo.SoilR_Weather_File")
DB = rbind(DB, DF)
dbWriteTable(connection, name = "SoilR_Weather_File", DB, overwrite = TRUE)

# Stratum File to be inserted - - - - - - - - - - - - - -

STRATUM_Edaphic_File = data.frame("Stratum" = "Storms River Testing",
                                  "Clay_Percentage" = 21.4800,
                                  "SOC_Percentage" = NA,
                                  "Bulk_Density" = NA,
                                  "Soil_Depth" = 30.0000,
                                  "SOC_Stock" = 97.4651)

DF = STRATUM_Edaphic_File
DB <- dbGetQuery(connection, "SELECT * FROM dbo.SoilR_Stratum_File")
DB = rbind(DB, DF)
dbWriteTable(connection, name = "SoilR_Stratum_File", DB, overwrite = TRUE)

# Calibrated Model File to be inserted - - - - - - - - - - - - - -

Calibrated_Model = data.frame("Model_Name" = "Test - JANSENVILLE",
                              "Stratum" = "Test",
                              "Weather_Station" = "JANSENVILLE",
                              "DPM" = 0.4510887,
                              "RPM" = 13.3954103,
                              "BIO" = 1.9193236,
                              "HUM" = 72.2360421,
                              "IOM" = 9.0259982)

DF = Calibrated_Model
DB <- dbGetQuery(connection, "SELECT * FROM dbo.SoilR_Calibrated_Model")
DB = rbind(DB, DF)
dbWriteTable(connection, name = "SoilR_Calibrated_Model", DB, overwrite = TRUE)
