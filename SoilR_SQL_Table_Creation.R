# 0. File Information ----

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 
# The purpose of this R script is to create SQL tables in our database which will store
# the relevant INPUT data which is used in:
#   a) RothC Model Calibration
#   b) RothC Model Simulation
# 
# The code below is meant for one-time usage where SQL tables are created.
# SQL table updates and data upserts will be written into a separate script
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# 1. Required Packages/Dependencies ----

library(pacman)
p_load(odbc,tidyverse,RODBC,DBI,dplyr, keyring)

# 2. Connecting to the SQL database securely ----
# 
# We need to form a SQL database connection without sharing credentials within the script.
# Two methods may be applicable here:
#   1. Use a secrets keyring
#   2. Use interactive prompts

# KEYRING METHOD
# 
# NB Keyring has been created locally and is password protected. 

keyring_unlock("cng_SQL_Credentials")
connection <- odbc::dbConnect(odbc(),Driver = "ODBC Driver 17 for SQL Server",
                              Server = key_get("SQL_Server", keyring ="cng_SQL_Credentials"),
                              Database = key_get("Database", keyring ="cng_SQL_Credentials"),
                              UID = key_get("Username", keyring ="cng_SQL_Credentials"),
                              PWD = key_get("Password", keyring ="cng_SQL_Credentials"))
keyring_lock("cng_SQL_Credentials")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# INTERACTIVE PROMPTING
# connection <- odbc::dbConnect(odbc(),Driver = "ODBC Driver 18 for SQL Server",
#                               Server = rstudioapi::showPrompt(title = "SQL Server", message = "Specify Server", default = ""),
#                               Database = rstudioapi::showPrompt(title = "Database", message = "Specify Database", default = ""),
#                               UID = rstudioapi::showPrompt(title = "Username", message = "Enter Username", default = ""),
#                               PWD = rstudioapi::showPrompt(title = "Password", message = "Enter Password", default = ""))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# 3. Create SQL tables for SoilR Input Data

# 3.1 Weather File Input Table

create_table_sql <- "
  CREATE TABLE SoilR_Weather_File (
    Month INT,
    Temperature FLOAT,
    Precipitation FLOAT,
    Evapotranspiration FLOAT,
    Weather_Station VARCHAR(MAX)
  )
"

# Execute the SQL statement to create the table
dbExecute(connection, create_table_sql)

# 3.2 Stratum File Input Table

create_table_sql <- "
  CREATE TABLE SoilR_Stratum_File (
    Stratum VARCHAR(MAX),
    Clay_Percentage FLOAT,
    SOC_Percentage FLOAT,
    Bulk_Density FLOAT,
    Soil_Depth FLOAT,
    SOC_Stock FLOAT
  )
"

# Execute the SQL statement to create the table
dbExecute(connection, create_table_sql)

# 3.3 ALM File Input Table

create_table_sql <- "
  CREATE TABLE SoilR_ALM_File (
    FarmID INT,
    FieldID INT,
    Scenario VARCHAR(MAX),
    Month INT,
    Bare VARCHAR(MAX),
    Cinput FLOAT,
    FYM FLOAT,
    Irrigation FLOAT
  )
"

# Execute the SQL statement to create the table
dbExecute(connection, create_table_sql)

# 3.4 Calibrated Model Input Table

create_table_sql <- "
  CREATE TABLE SoilR_Calibrated_Model (
    Model_Name VARCHAR(MAX),
    Stratum VARCHAR(MAX),
    Weather_Station VARCHAR(MAX),
    Cinput FLOAT,
    DPM FLOAT,
    RPM FLOAT,
    BIO FLOAT,
    HUM FLOAT,
    IOM FLOAT
  )
"

# Execute the SQL statement to create the table
dbExecute(connection, create_table_sql)

dbDisconnect(connection)
