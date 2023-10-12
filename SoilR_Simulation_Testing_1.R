# 0. File Information ----

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 
# This R script is the first attempt at RothC model SIMULATION using SoilR.
# The code below follows a brief tutorial set out at the following resource:
# https://www.bgc-jena.mpg.de/TEE/basics/2015/11/19/RothC/
# 
# NB: Tutorial covers Method 1: "Using clay content, litter inputs, and climate data"
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

# 2. Weather Data ----

Temp = data.frame("Month" = 1:12, 
                  "Temp" = c(24.75, 24.75, 22.8, 19.7, 15.55, 12.75, 12.2, 14.45, 16.95, 19.15, 21.35, 23.3))

Precip = data.frame("Month" = 1:12, 
                    "Precip" = c(25.83, 27.65, 44.31, 27.65, 11.77, 7.9, 13.69, 16.54, 13.69, 21.23, 29.46, 25.83)) # NOTE that these values are excluding any irrigation

<<<<<<< HEAD
# Evp = data.frame("Month" = 1:12, 
#                 "Evp" = c(265.1, 223.1, 191.4, 142.8, 112.9, 89.29, 97.17, 125.9, 158.6, 200.3, 248.6, 274.9))

 Evp = data.frame("Month" = 1:12, 
                  "Evp" = c(8.55, 7.97, 6.17, 4.76, 3.64, 2.98, 3.13, 4.06, 5.29, 6.46, 8.29, 8.87))
 
# Evp=data.frame(Month=1:12, Evp=c(12, 18, 35, 58, 82, 90, 97, 84, 54, 31,
#                                 14, 10))
=======
Evp = data.frame("Month" = 1:12,
                 "Evp" = c(265.1, 223.1, 191.4, 142.8, 112.9, 89.29, 97.17, 125.9, 158.6, 200.3, 248.6, 274.9))

# Evp=data.frame(Month=1:12, Evp=c(12, 18, 35, 58, 82, 90, 97, 84, 54, 31,
#                                  14, 10))
>>>>>>> 108d16104d9768d8361b233583ece7b75b04b52c

# 3. Edaphic Data ----

soil.thick = 30 # Soil thickness (organic layer topsoil) (cm)
SOC = 97.4651688163155 # Soil Organic Carbon Stock (Mg/ha). NB: Mg refers to Megagram = metric Tonne.
clay = 21.4800892252279 # Percent clay (%)
Cinputs = 0.5 # Annual C inputs to soil (Mg/ha/yr). NB: This is the value that we need to check on. For now, Ill assume its the same as that used to calibrate

# 4. RothC Simulation Duration ----

years = seq(1/12,500,by=1/12) 

# 5. Effects of Climate on Decomposition ----

fT = fT.RothC(Temp[,2]) # Temperature effects per month

fW = fW.RothC(P=(Precip[,2]), E=(Evp[,2]), 
              S.Thick = soil.thick, 
              pClay = clay, 
              pE = 0.75, bare = FALSE)$b # Moisture effects per month

xi.frame = data.frame(years,
                      rep(fT*fW, length.out = length(years)))

# 6. Size of the Inert Organic Matter pool (IOM) ----

FallIOM=0.049*SOC^(1.139) #IOM using Falloon method

# 7. Initializing the model and solving C stocks for each pool (Not sure what to call this section) ----

Model1=RothCModel(t = years,
                  C0=c(DPM = 0, RPM = 0, BIO = 0, HUM = 0, IOM = FallIOM),
                  In = Cinputs, 
                  clay = clay, 
                  xi = xi.frame) # Loads the model

Ct1 = getC(Model1) # Calculates stocks for each pool per month

# Plotting the results...

matplot(years, 
        Ct1, 
        type="l", 
        lty=1, 
        col=1:5,
        xlab = "Time (years)", 
        ylab="C stocks (Mg/ha)")

legend("topleft", 
       c("DPM", "RPM", "BIO", "HUM", "IOM"),
       lty=1, 
       col=1:5, 
       bty="n")

# 8. Final pool sizes after 500 year spinup ----

poolSize1=as.numeric(tail(Ct1,1))
names(poolSize1)<-c("DPM", "RPM", "BIO", "HUM", "IOM")
poolSize1
