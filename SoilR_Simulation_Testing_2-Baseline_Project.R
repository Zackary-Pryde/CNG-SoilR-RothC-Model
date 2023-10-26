# 0. File Information ----

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 
# This R script is the first attempt at RothC model SIMULATION using SoilR.

# Use the existing RothC files to attempt to reach the same solution for Field 1  
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

### 1. Dependencies ----

library(pacman)
p_load(raster, rgdal, ncdf4, SoilR, abind, soilassessment, Formula)

# The code below follows a brief tutorial set out at the following resource:
# https://www.bgc-jena.mpg.de/TEE/basics/2015/11/19/RothC/
# 
# NB: Tutorial covers Method 1: "Using clay content, litter inputs, and climate data"
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# 2. Clay content, litter inputs, and climate data ----

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# This script will use data inputs and outputs for the PD1 MR1 Farm: Storms River - Field 1.

# NB I am setting up these data to track that used for the PD1 MR1 Farm: Storms River.
#   Weather File = JANSENVILLE
#   Stratum = Knysna-Amatole montane forests / Cfb : Temperate, no dry season, warm summer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


### 2. Weather Data ----

Temp = data.frame("Month" = 1:12, 
                  "Temp" = c(24.8, 24.8, 22.8, 19.7, 15.6, 12.8, 12.2, 14.4, 16.9, 19.1, 21.4, 23.3))

Precip = data.frame("Month" = 1:12, 
                    "Precip" = c(25.8, 27.6, 44.3, 27.6, 11.8, 7.9, 13.7, 16.5, 13.7, 21.2, 29.5, 25.8)) # NOTE that these values are excluding any irrigation

Evp = data.frame("Month" = 1:12, 
                 "Evp" = c(265.1, 223.1, 191.4, 142.8, 112.9, 89.3, 97.2, 125.9, 158.6, 200.3, 248.6, 274.9))

# Evp = data.frame("Month" = 1:12, 
#                  "Evp" = c(8.55, 7.97, 6.17, 4.76, 3.64, 2.98, 3.13, 4.06, 5.29, 6.46, 8.29, 8.87))

### 3. Edaphic Data ----

soil.thick = 30 # Soil thickness (organic layer topsoil) (cm)
clay = 21.48 # Percent clay (%)

### 4. RothC Simulation Duration ----

# set to 2 years, 2020 and 2021

years = seq(1/12,2,by=1/12)

### 5. Effects of Climate on Decomposition ----

# temperature is the same in the baseline and project scenarios

fT.RothC_CNG = function (Temp) {
  47.91/(1 + exp(106.06/(ifelse(Temp >= -18.27, Temp, NA) + 18.27)))
}

fT = fT.RothC_CNG(Temp[,2])

# moisture effects are different between the baseline and project scenarios due to change in bare vs. covered

fW.base = fW.RothC(P=(Precip[,2]), E=(Evp[,2]),
                   S.Thick = soil.thick,
                   pClay = clay,
                   pE = 1, bare = FALSE)$b # logical value assumes same bare vs. covered in all months

fW.proj = fW.RothC(P=(Precip[,2]), E=(Evp[,2]),
                   S.Thick = soil.thick,
                   pClay = clay,
                   pE = 1, bare = FALSE)$b # logical value assumes same bare vs. covered in all months


# adjustments to the fW equation to allow for monthly bare vs cover - code from FAO GSOC
# https://fao-gsp.github.io/GSOCseq/stage-2-running-the-model-1.html#overview-of-the-main-commands-to-perform-the-rothc-calculations

base.soil.cover = data.frame("Month" = 1:12,
                        "Cov" = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))

proj.soil.cover = data.frame("Month" = 1:12,
                             "Cov" = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))

# bare1 = data.frame("Month" = 1:12,
#                    "Soil_cover" = c(TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE))

#Moisture effects per month . 
fw1func<-function(P, E, S.Thick = 30, pClay = 32.0213, pE = 1, bare) 
{
  
  M = P - E * pE
  Acc.TSMD = NULL
  for (i in 2:length(M)) {
    B = ifelse(bare[i] == FALSE, 1, 1.8)
    Max.TSMD = -(20 + 1.3 * pClay - 0.01 * (pClay^2)) * (S.Thick/23) * (1/B)
    Acc.TSMD[1] = ifelse(M[1] > 0, 0, M[1])
    if (Acc.TSMD[i - 1] + M[i] < 0) {
      Acc.TSMD[i] = Acc.TSMD[i - 1] + M[i]
    }
    else (Acc.TSMD[i] = 0)
    if (Acc.TSMD[i] <= Max.TSMD) {
      Acc.TSMD[i] = Max.TSMD
    }
  }
  b = ifelse(Acc.TSMD > 0.444 * Max.TSMD, 1, (0.2 + 0.8 * ((Max.TSMD - 
                                                              Acc.TSMD)/(Max.TSMD - 0.444 * Max.TSMD))))
  b<-clamp(b,lower=0.2)
  return(data.frame(Acc.TSMD, b, Max.TSMD))
}

fW.base_2<- fw1func(P=(Precip[,2]), E=(Evp[,2]), S.Thick = soil.thick, pClay = clay, pE = 1, bare = base.soil.cover$Cov)$b 

fW.proj_2<- fw1func(P=(Precip[,2]), E=(Evp[,2]), S.Thick = soil.thick, pClay = clay, pE = 1, bare = proj.soil.cover$Cov)$b 

### 6. Combining climate and edaphic effects ----
# set the soil cover factor (RothC documentation)
# if soil is vegetated c = 0.6
# if soil is bare c = 1.0 

base.cover.factor <- base.soil.cover[,2]
base.cover.factor[which(base.cover.factor == TRUE)] = 1
base.cover.factor[which(base.cover.factor == 0)] = 0.6
base.cover.factor

proj.cover.factor <- proj.soil.cover[,2]
proj.cover.factor[which(proj.cover.factor == TRUE)] = 1
proj.cover.factor[which(proj.cover.factor == 0)] = 0.6
proj.cover.factor


xi.frame.base = data.frame(years,
                           rep(fT*fW.base_2*base.cover.factor, length.out = length(years)))

xi.frame.proj = data.frame(years,
                           rep(fT*fW.proj_2*proj.cover.factor, length.out = length(years)))

### 7. Set up the baseline and project simulations ----

# ks: default decomposition rate constants for different pools are used (the same in baseline and project)
# C0: initial amount of carbon in each pool obtained from the spin-up / initialised model (always the same in baseline and project)
# In: plant residue inputs (different between baseline and project in this case)
# FYM: farm yard manure input (the same between baseline and project in this case)
# DR: ratio of decomposable plant material to resistant plant material (default value)

base.inputs = data.frame("Month" = 1:12,
                         "In" = rep(1.18, 12))

proj.inputs = data.frame("Month" = 1:12,
                         "In" = rep(1.56, 12))


base.model <- RothCModel(t = years, ks = c(k.DPM = 10, k.RPM = 0.3, k.BIO = 0.66, k.HUM = 0.02, k.IOM = 0),
                         C0 = c(0.6396, 13.2247, 1.8943, 72.5269, 9.0260), In = as.data.frame(base.inputs$In), FYM = as.data.frame(rep(0,12)), DR = 1.44,
                         clay = clay, xi = xi.frame.base, solver = euler ,pass = TRUE)

proj.model <- RothCModel(t= years, ks = c(k.DPM = 10, k.RPM = 0.3, k.BIO = 0.66, k.HUM = 0.02, k.IOM = 0),
                         C0 = c(0.6396, 13.2247, 1.8943, 72.5269, 9.0260), In = as.data.frame(proj.inputs$In), FYM = as.data.frame(rep(0,12)), DR = 1.44,
                         clay = clay, xi = xi.frame.proj, solver = euler ,pass = TRUE)  
?RothCModel()

base.model.output = getC(base.model)
poolsize.base = as.numeric(tail(base.model.output,1))
names(poolsize.base) <- c("DPM", "RPM", "BIO", "HUM", "IOM")
poolsize.base
sum(poolsize.base)

proj.model.output = getC(proj.model)
poolsize.proj = as.numeric(tail(proj.model.output,1))
names(poolsize.proj) <- c("DPM", "RPM", "BIO", "HUM", "IOM")
poolsize.proj
sum(poolsize.proj)