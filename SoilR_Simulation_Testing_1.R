# 0. File Information ----

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 
# This R script is the first attempt at RothC model SIMULATION using SoilR.
# The code below follows a brief tutorial set out at the following resource:
# https://www.bgc-jena.mpg.de/TEE/basics/2015/11/19/RothC/
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# 1. Dependencies ----

library(pacman)
p_load(SoilR)

# 2. Clay content, litter inputs, and climate data ----

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# NB I am setting up these data to track that used for the PD1 MR1 Farm: Storms River.
#   Weather File = JANSENVILLE
#   Stratum = Knysna-Amatole montane forests / Cfb : Temperate, no dry season, warm summer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Weather File Data

Temp = data.frame("Month" = 1:12, 
                  "Temp" = c(24.75, 24.75, 22.8, 19.7, 15.55, 12.75, 12.2, 14.45, 16.95, 19.15, 21.35, 23.3))

Precip = data.frame("Month" = 1:12, 
                    "Precip" = c(25.83, 27.65, 44.31, 27.65, 11.77, 7.9, 13.69, 16.54, 13.69, 21.23, 29.46, 25.83)) # NOTE that these values are excluding any irrigation

Evp = data.frame("Month" = 1:12, 
                 "Evp" = c(265.1, 223.1, 191.4, 142.8, 112.9, 89.29, 97.17, 125.9, 158.6, 200.3, 248.6, 274.9))
