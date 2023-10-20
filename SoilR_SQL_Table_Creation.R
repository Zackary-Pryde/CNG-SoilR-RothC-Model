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
connection <- odbc::dbConnect(odbc(),Driver = "ODBC Driver 18 for SQL Server",
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
