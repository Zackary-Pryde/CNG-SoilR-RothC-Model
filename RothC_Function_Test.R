# Required Packages/Dependencies ----
library(pacman)
p_load(raster, ncdf4, SoilR, abind, soilassessment, Formula, ggplot2, tidyverse, readr, odbc, RODBC, DBI, dplyr)

options(digits = 8)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - -
# NOTES ----
# Farm:               Basson Family Trust
# Paddock:            Field 2
# Area:               662.34904676
# Stratum:            Lowland fynbos and renosterveld / Cfb : Temperate, no dry season, warm summer
# Weather Station:    PORT-ELIZABETH
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - -

# Input Data ----

# SQL Credentials ----
connection <- odbc::dbConnect(odbc(),Driver = "ODBC Driver 17 for SQL Server",
                              Server = "sql-agricarbon-dev.database.windows.net",
                              Database = "sqldb-agricarbon-dev",
                              UID = "agricarbon_admin",
                              PWD = "Cl1mat3_N3utral@2022")

# Data acquisition and preparation ----

Paddock_UID_Input = "Vernon Scheepers/Storms River/Field2"

FFM_Information <- dbGetQuery(connection, paste0("SELECT * FROM dbo.RothC_DOS_Farm_Field_Master WHERE Paddock_UID = '", Paddock_UID_Input,"'"))

Weather_File_Input = dbGetQuery(connection, paste0("SELECT * FROM dbo.RothC_DOS_Weather_File WHERE Weather_Station = '", FFM_Information$Weather_Station,"'")) %>%
  dplyr::select(Month, "Temp" = Temperature, "Precip" = Precipitation, "Evp" = Evapotranspiration)

Edaphic_File_Input = dbGetQuery(connection, paste0("SELECT * FROM dbo.RothC_DOS_Stratum_File WHERE Stratum = '", FFM_Information$Stratum,"'")) %>%
  dplyr::select(Soil_Depth, "SOC_Stratum" = SOC_Stock, "ClayPerc_Stratum" = Clay_Percentage)

ALMP_BL = dbGetQuery(connection, paste0("SELECT * FROM dbo.RothC_DOS_ALM_File WHERE Paddock_UID = '", FFM_Information$Paddock_UID,"'", "AND Scenario = 'Baseline'")) %>%
  dplyr::select(Month, Bare, Cinput, FYM, Irrigation)

ALMP_PR = dbGetQuery(connection, paste0("SELECT * FROM dbo.RothC_DOS_ALM_File WHERE Paddock_UID = '", FFM_Information$Paddock_UID,"'", "AND Scenario = 'Project'")) %>%
  dplyr::select(Month, Bare, Cinput, FYM, Irrigation)

Calibrated_Model_DOS = dbGetQuery(connection, paste0("SELECT * FROM dbo.RothC_DOS_Calibrated_Model WHERE Model_Name = '", FFM_Information$Calibrated_Model,"'")) %>%
  dplyr::select(DPM,RPM,BIO,HUM,IOM)

Weather_File = Weather_File_Input
Edaphic_File = Edaphic_File_Input
ALMP_File_BL = ALMP_BL
ALMP_File_PR = ALMP_PR
Calibrated_Model = Calibrated_Model_DOS

# TEMPERATURE

#Temp = weather_file$Temp
Temp = Weather_File$Temp
fT = 47.91/(1 + exp(106.06/(ifelse(Temp > -18.27, Temp, NA) + 18.27))) # R documentation incorrect. Cant have zero denomenator -> Must be strictly greater than...

# WEATHERING

P = Weather_File$Precip + ALMP_File_BL$Irrigation
E = Weather_File$Evp
S.Thick = Edaphic_File$Soil_Depth
pClay = Edaphic_File$ClayPerc_Stratum
pE = 1
bare = ALMP_File_BL$Bare

M = P - E * pE

Acc.TSMD = NULL
for (i in 2:length(M)) {
  
  B = ifelse(bare == FALSE, 1, 1.8)
  Max.TSMD = min(-(20 + (1.3 * pClay) - (0.01 * (pClay^2))) * (S.Thick/23) * (1/B))
  Acc.TSMD[1] = ifelse(M[1] > 0, 0, M[1])
  
  if (Acc.TSMD[i - 1] + M[i] < 0) {
    Acc.TSMD[i] = Acc.TSMD[i - 1] + M[i]
  } else 
    (Acc.TSMD[i] = 0)
  
  if (Acc.TSMD[i] <= Max.TSMD) {
    Acc.TSMD[i] = Max.TSMD
  }
}

b = ifelse(Acc.TSMD > 0.444 * Max.TSMD, 1, (0.2 + (1 - 0.2) * ((Max.TSMD - Acc.TSMD)/(Max.TSMD - 0.444 * Max.TSMD)))) # Changed to match the documentation
b<-clamp(b,lower=0.2)
data.frame(Acc.TSMD, b, Max.TSMD)

# Bare Cover modifyer (i.e., Crop Retainment)

Cov2 = ALMP_BL[,c("Month","Bare")] %>%
  mutate(Bare = ifelse(Bare == TRUE, 1,0.6))

fC <- Cov2[,2]

(a = fT)
(b = b)
(c = fC)
(abc = a*b*c)


years = seq(0,1,1/12) # Always 1 year

# (1 - e^(-abckt))
Effects = data.frame("a" = fT,
                     "b" = b,
                     "c" = fC,
                     "k_DPM" = 10/12,
                     "k_RPM" = 0.3/12,
                     "k_BIO" = 0.66/12,
                     "k_HUM" = 0.02/12,
                     "k_IOM" = 0/12)

#xi.frame = data.frame(years[2:length(years)], rep(fT*b*fC,length.out=length(years[2:length(years)])))
xi.frame = data.frame(seq(0,1,1/12)[2:13],fT*b*fC)

calibrated_model = Calibrated_Model_DOS

# RUN THE MODEL FROM SOILR
Model3_spin = RothCModel(t = years,
                         C0 = c(calibrated_model$DPM, 
                                calibrated_model$RPM, 
                                calibrated_model$BIO, 
                                calibrated_model$HUM, 
                                calibrated_model$IOM),
                         ks = c(k.DPM = 10.0/12, k.RPM = 0.3/12, k.BIO = 0.66/12, k.HUM = 0.02/12, k.IOM = 0/12),
                         In = data.frame(ALMP_BL$Month,ALMP_BL$Cinput),
                         #In = sum(ALMP_BL$Cinput),
                         FYM = data.frame(ALMP_BL$Month,ALMP_BL$FYM),
                         #FYM = sum(ALMP_BL$FYM),
                         DR = 1.44, #TRY 0.59
                         clay = Edaphic_File_Input$ClayPerc_Stratum,
                         xi = xi.frame, 
                         pass = TRUE, 
                         solver = euler)

Ct3_spin=getC(Model3_spin)
#getC14(object = Model3_spin)
#getF14C(object = Model3_spin)
#getF14(object = Model3_spin)
#getTimes(Model3_spin)
getReleaseFlux(Model3_spin)

Model3_spin[c("time", "C")]


colnames(Ct3_spin)<-c("DPM", "RPM", "BIO", "HUM", "IOM")
Model_Result = as.data.frame(Ct3_spin)
Model_Result$SOC_Stock = Model_Result$DPM + Model_Result$RPM + Model_Result$BIO + Model_Result$HUM + Model_Result$IOM
Model_Result

calibrated_model = Model_Result[13,1:5]


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Test theory from RothC application documentation



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Test theory from RothC in Sierra et al

# C input scalar
I = sum(ALMP_File_BL$Cinput)

# Matrix 1
g = 0.59
M1 = matrix(c(g, 
              1-g, 
              0, 
              0, 
              0),
            nrow = 5, 
            ncol = 1,
            byrow = TRUE)

# Matrix 2
k1 = 10
k2 = 0.3
k3 = 0.66
k4 = 0.02

x_prop = 1.67*(1.85 + 1.60*(exp(-0.0786*Edaphic_File_Input$ClayPerc_Stratum)))

a_31 = k1*(0.46/(x_prop+1))
a_32 = k2*(0.46/(x_prop+1))
a_33 = k3*(0.46/(x_prop+1))
a_34 = k4*(0.46/(x_prop+1))

a_41 = k1*(0.54/(x_prop+1))
a_42 = k2*(0.54/(x_prop+1))
a_43 = k3*(0.54/(x_prop+1))
a_44 = k4*(0.54/(x_prop+1))

M2 = matrix(c(-k1,0,0,0,0,
              0,-k2,0,0,0,
              a_31,a_32,-k3+a_33,a_34,0,
              a_41,a_42,a_43,-k4+a_44,0,
              0,0,0,0,0),
            nrow = 5,
            ncol = 5,
            byrow = TRUE)

# Matrix 3

C1 = Calibrated_Model_DOS$DPM
C2 = Calibrated_Model_DOS$RPM
C3 = Calibrated_Model_DOS$BIO
C4 = Calibrated_Model_DOS$HUM
C5 = Calibrated_Model_DOS$IOM

M3 = matrix(c(C1,
              C2,
              C3,
              C4,
              C5),
            nrow = 5,
            ncol = 1,
            byrow = TRUE)

# Matrix Algebra for result at one timestep

I
M1
M2
M3

dC_by_dt = I*M1 + M2 %*% M3
dC_by_dt/12

t2 = M3 + (dC_by_dt/12)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# C input scalar
I = sum(ALMP_File_BL$Cinput)

# Matrix 1
g = 0.59
M1 = matrix(c(g, 
              1-g, 
              0, 
              0, 
              0),
            nrow = 5, 
            ncol = 1,
            byrow = TRUE)

# Matrix 2
k1 = 10
k2 = 0.3
k3 = 0.66
k4 = 0.02

x_prop = 1.67*(1.85 + 1.60*(exp(-0.0786*Edaphic_File_Input$ClayPerc_Stratum)))

a_31 = k1*(0.46/(x_prop+1))
a_32 = k2*(0.46/(x_prop+1))
a_33 = k3*(0.46/(x_prop+1))
a_34 = k4*(0.46/(x_prop+1))

a_41 = k1*(0.54/(x_prop+1))
a_42 = k2*(0.54/(x_prop+1))
a_43 = k3*(0.54/(x_prop+1))
a_44 = k4*(0.54/(x_prop+1))

M2 = matrix(c(-k1,0,0,0,0,
              0,-k2,0,0,0,
              a_31,a_32,-k3+a_33,a_34,0,
              a_41,a_42,a_43,-k4+a_44,0,
              0,0,0,0,0),
            nrow = 5,
            ncol = 5,
            byrow = TRUE)

# Matrix 3

M3 = t2

# Matrix Algebra for result at one timestep

I
M1
M2
M3

dC_by_dt = I*M1 + M2 %*% M3
dC_by_dt/12

t3 = M3 + (dC_by_dt/12)
