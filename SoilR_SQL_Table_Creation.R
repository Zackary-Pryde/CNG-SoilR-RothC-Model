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
p_load(raster, rgdal, ncdf4, SoilR, abind, soilassessment, Formula)
