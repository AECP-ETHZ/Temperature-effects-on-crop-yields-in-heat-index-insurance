# =========================================================================================================
# ---------------------------------------------------------------------------------------------------------
#
# Temperature Effects on Crop Yields in Heat Index Insurance
#
# Supplementary R code for main results
# 
# Authors: Janic Bucheli, Tobias Dalhaus, Robert Finger
#
# Corresponding author: Janic Bucheli (jbucheli@ethz.ch)
#
# Please check our github repository for more codes and information
# https://github.com/AECP-ETHZ/Temperature-effects-on-crop-yields-in-heat-index-insurance
#
# Citation: Bucheli, J., Dalhaus T., & Finger R. 2021. Temperature Effects on Crop Yields in Heat Index Insurance.
#           Food Policy.
# 
# ---------------------------------------------------------------------------------------------------------
# =========================================================================================================

# This file shows the codes to calculate our key results. Please find additional codes for this publication
# on our github repository: https://github.com/AECP-ETHZ/Temperature-effects-on-crop-yields-in-heat-index-insurance 
#
# Content:
# 1) Loading packages
# 2) Specifications, loading and preparing data
# 3) Calculation of hourly temperature exposure during the risk period (defined by phenology)
# 4) Restructuring datasets
# 5) Expected utility model: Functions
# 6) Expected utility model: Uninsured status (base-line scenario)
# 7) Winter wheat: out-of-sample calibration (with restricted cubic splines) and testing for three knots
# 8) Winter rapeseed: out-of-sample calibration (with restricted cubic splines) and testing for three knots
# 9) Winter wheat: out-of-sample calibration (with restricted cubic splines) and testing for five knots
# 10) Winter rapeseed: out-of-sample calibration (with restricted cubic splines) and testing for five knots

# =========================================================================================================
# 1. Loading packages & increasing memory size
# =========================================================================================================
'
install.packages("raster")
install.packages("dplyr")
install.packages("data.table")
install.packages("reshape2")
install.packages("Hmisc")
install.packages("lubridate")
install.packages("plyr")
install.packages("esmisc")
install.packages("sandwich")
install.packages("MASS")
install.packages("e1071")
'
library(raster)
library(dplyr)
library(data.table)
library(reshape2)
library(Hmisc)
library(lubridate)
library(plyr)
library(esmisc)
library(sandwich)
library(MASS)
library(e1071)

# =========================================================================================================
# 2. Specifications, loading and preparing data
# =========================================================================================================

# -----------------------------------------
# Yield data & farm location
# -----------------------------------------

# Define crop and last year
crops     <- c("wheat","rapeseed") 
last_year <- 2018

# Read yield data: [[1]] == winter wheat; [[2]] == winter rapeseed
farm_yields <- list()
for (k in 1:length(crops)){
  farm_yields[[k]] <- read.csv(paste("Yield/",crops[k],"_untrended18.csv", sep=""))
}

# Assign colnames to yield data (longitude, latitude and years)
years_panel                <- seq(last_year,1995,-1)
colnames(farm_yields[[1]]) <- c("farm.lon", "farm.lat", years_panel)
colnames(farm_yields[[2]]) <- c("farm.lon", "farm.lat", years_panel)

# Create spatial points raster files with farm locations
farm_coordinates <- list()
for (k in 1:length(crops)){
  temp1 <- farm_yields[[k]][,1:2]
  coordinates(temp1)    <- c("farm.lon", "farm.lat")
  proj4string(temp1)    <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  farm_coordinates[[k]] <- temp1
  rm(temp1)
}

# -----------------------------------------
# Temperature data
# -----------------------------------------

# Define start & end date of panel (pick first and las day of year, respectively)
start_panel <- as.Date("1995-01-01", format ="%Y-%m-%d")
end_panel   <- as.Date("2018-12-31", format="%Y-%m-%d")

# Read the Tmin and Tmax data with stack (faster than brick)
Tmin <- stack("Meteo/TMin/tn_ens_mean_0.1deg_reg_v20.0e.nc")
Tmax <- stack("Meteo/TMax/tx_ens_mean_0.1deg_reg_v20.0e.nc")

# Subset the historical temperature data w.r.t. yield panel
# Each layer in Tmin & Tmax is a year.
first_layer_date <- as.Date(substring(Tmin[[1]]@data@names,2), format = "%Y.%m.%d")
last_layer_date  <- as.Date(substring(Tmin[[dim(Tmin)[3]]]@data@names,2), format = "%Y.%m.%d")
seq_layers_dates <- seq(as.Date(first_layer_date), as.Date(last_layer_date), "days")
sub_layers_panel <- which(seq_layers_dates >= start_panel & seq_layers_dates <= end_panel)

sub_Tmin <- subset(Tmin, sub_layers_panel)
sub_Tmax <- subset(Tmax, sub_layers_panel)

# Remove large dataset from global environment
rm(Tmin,Tmax, first_layer_date,last_layer_date)

# -----------------------------------------
# Match farm locations with weather data
# -----------------------------------------

# Tmin
Tmin_farm <- list()
for (k in 1:length(crops)){
  temp           <- as.data.frame(matrix(NA,ncol= length(seq_layers_dates[sub_layers_panel]), nrow=nrow(farm_yields[[k]])))
  colnames(temp) <- seq_layers_dates[sub_layers_panel]
  
  for (t in 1:length(sub_layers_panel)){
    temp[1:nrow(farm_yields[[k]]),t] <- raster::extract(sub_Tmin[[t]], farm_coordinates[[k]])
  }
  Tmin_farm[[k]] <- temp 
  rm(temp)
}

# Tmax
Tmax_farm <- list()
for (k in 1:length(crops)){
  temp           <- as.data.frame(matrix(NA,ncol= length(seq_layers_dates[sub_layers_panel]), nrow=nrow(farm_yields[[k]])))
  colnames(temp) <- seq_layers_dates[sub_layers_panel]
  
  for (t in 1:length(sub_layers_panel)){
    temp[1:nrow(farm_yields[[k]]),t] <- raster::extract(sub_Tmax[[t]], farm_coordinates[[k]])
  }
  Tmax_farm[[k]] <- temp 
  rm(temp)
}

# Remove large dataset from global environment
rm(sub_Tmin,sub_Tmax)

# -----------------------------------------
# Load phenology data (for risk period)
# -----------------------------------------

# Wheat: start = stem elongation, end = milk ripeness
PhenoFarm_wheat_start <- read.csv("Phenology/PhenoFarm/wheat_stemelong_new.csv")[,-c(1)]
PhenoFarm_wheat_end   <- read.csv("Phenology/PhenoFarm/wheat_milkripe_new.csv") [,-c(1)]
colnames(PhenoFarm_wheat_start) <- seq(1995,last_year,1)
colnames(PhenoFarm_wheat_end)   <- seq(1995,last_year,1)

# rapeseed: start = bud formation, end = ripeness
PhenoFarm_rapeseed_start <-  read.csv("Phenology/PhenoFarm/rapeseed_budformation_new.csv")[,-c(1)]
PhenoFarm_rapeseed_end   <-  read.csv("Phenology/PhenoFarm/rapeseed_ripeness_to18.csv")[,-c(1,2)]
colnames(PhenoFarm_rapeseed_start) <- seq(1995,last_year,1)
colnames(PhenoFarm_rapeseed_end)   <- seq(1995,last_year,1)

# ==============================================================================================================
# 3. Calculation of hourly temperature exposure during the risk period 
# ==============================================================================================================

# -----------------------------------------
# Winter wheat
# -----------------------------------------

# Hourly steps in unit radian
# By 2pi/24 to set the hour steps, and add/substract 2*pi*1/48 to estimate at the middle of the hour
# First seq() is for first half-day; second seq() for second half-day
radian_observations <- c(seq(((-0.5*pi)+(2*pi/48)),(0.5*pi-(2*pi/48)), by=(2*pi/24)), seq((( 0.5*pi)+(2*pi/48)),(1.5*pi-(2*pi/48)), by=(2*pi/24)))

# List for each farm: farm_wheat_exposure [[i]] [[t]] with i == farm, t == year
farm_wheat_exposure <- list()

for (i in 1:nrow(farm_yields[[1]])){
  
  # Temperature exposure of farm i n years t; each element of in "years_exposure" is a year of farm i.
  years_exposure <- list()
  
  for (y in 1:(last_year - 1995 + 1)){
    
    start_risk <- as.Date(PhenoFarm_wheat_start[i,y], formate = "%Y-%m-%d")
    # End the risk period before entrance into new growth phase (-1)
    end_risk   <- as.Date(PhenoFarm_wheat_end [i,y], formate = "%Y-%m-%d")-1
    
    # Create data frame for hourly temperature readings
    daily_temperature_readings            <- as.data.frame(matrix(NA, nrow= as.numeric(end_risk - start_risk +1), ncol=24))
    row.names(daily_temperature_readings) <- seq(start_risk, end_risk,1)
    colnames(daily_temperature_readings)  <- seq(0,23,1)
    
    for (d in 1:nrow(daily_temperature_readings)){
      
      # Get the relevant temperature data of this day
      Tmin_t  <- Tmin_farm[[1]][i,as.character(row.names(daily_temperature_readings)[d])]
      Tmin_t1 <- Tmin_farm[[1]][i,as.character(as.Date(row.names(daily_temperature_readings)[d])+1)]
      Tmax_t  <- Tmax_farm[[1]][i,as.character(row.names(daily_temperature_readings)[d])]
      
      # Get amplitude of first half (fh) and second half (sh) of the day
      amplitude_fh <- (Tmax_t - Tmin_t) / 2
      amplitude_sh <- (Tmax_t - Tmin_t1)/ 2
      
      # Get daily average of first half and second half of the day
      Tmean_fh <- (Tmax_t + Tmin_t) / 2
      Tmean_sh <- (Tmax_t + Tmin_t1)/ 2
      
      # Calculate the hourly temperature of the first half of the day (first sine curve)
      # Daily observation starts 30' after Tmin
      
      # Get hourly temperatures for the first half of the day (first sine curve)
      for (fh in 1:12){
        daily_temperature_readings[d,fh] <- Tmean_fh + amplitude_fh * sin(radian_observations[fh])
      }
      
      # Get hourly temperatures for the second half of the day (second sine curve)
      for (sh in 1:12){
        daily_temperature_readings[d,12+sh] <- Tmean_sh + amplitude_sh * sin(radian_observations[12+sh])
      }
    }
    
    years_exposure[[y]] <- daily_temperature_readings
    rm(daily_temperature_readings)
    years_exposure[[y]]$date<-seq(start_risk, end_risk,1)
  }
  farm_wheat_exposure[[i]] <- years_exposure 
  rm(years_exposure)
}

# -----------------------------------------
# Winter rapeseed
# -----------------------------------------

# List for each farm: farm_rapeseed_exposure [[i]] [[t]] with i == farm, t == year
farm_rapeseed_exposure <- list()

for (i in 1:nrow(farm_yields[[2]])){
  
  # Temperature exposure of farm i n years t; each element of in "years_exposure" is a year of farm i.
  years_exposure <- list()
  
  for (y in 1:(last_year - 1995 + 1)){
    
    start_risk <- as.Date(PhenoFarm_rapeseed_start[i,y], formate = "%Y-%m-%d")
    # End the risk period before entrance into new growth phase (-1)
    end_risk   <- as.Date(PhenoFarm_rapeseed_end [i,y], formate = "%Y-%m-%d")-1
    
    # Create data frame for hourly temperature readings
    daily_temperature_readings <- as.data.frame(matrix(NA, nrow= as.numeric(end_risk - start_risk +1), ncol=24))
    row.names(daily_temperature_readings) <- seq(start_risk, end_risk,1)
    colnames(daily_temperature_readings) <- seq(0,23,1)
    
    for (d in 1:nrow(daily_temperature_readings)){
      
      # Get the relevant temperature data of this day
      Tmin_t  <- Tmin_farm[[2]][i,as.character(row.names(daily_temperature_readings)[d])]
      Tmin_t1 <- Tmin_farm[[2]][i,as.character(as.Date(row.names(daily_temperature_readings)[d])+1)]
      Tmax_t  <- Tmax_farm[[2]][i,as.character(row.names(daily_temperature_readings)[d])]
      
      # Get amplitude of first half (fh) and second half (sh) of the day (amplitude of half-daily sine function)
      amplitude_fh <- (Tmax_t - Tmin_t) / 2
      amplitude_sh <- (Tmax_t - Tmin_t1)/ 2
      
      # Get daily average of first half and second half of the day
      Tmean_fh <- (Tmax_t + Tmin_t) / 2
      Tmean_sh <- (Tmax_t + Tmin_t1)/ 2
      
      # Calculate the hourly temperature of the first half of the day (first sine curve)
      # Daily observation starts at Tmin (radian = -pi(12)) 
  
      # Get hourly temperatures for the first half of the day (first sine curve)
      for (fh in 1:12){
        daily_temperature_readings[d,fh] <- Tmean_fh + amplitude_fh * sin(radian_observations[fh])
      }
      
      # Get hourly temperatures for the second half of the day (second sine curve)
      for (sh in 1:12){
        daily_temperature_readings[d,12+sh] <- Tmean_sh + amplitude_sh * sin(radian_observations[12+sh])
      }
    }
    
    years_exposure[[y]]      <- daily_temperature_readings
    years_exposure[[y]]$date <-seq(start_risk, end_risk,1)
    rm(daily_temperature_readings)
  }
  farm_rapeseed_exposure[[i]]  <- years_exposure 
  rm(years_exposure)
  print(i / nrow(farm_yields[[2]]))
}

# Remove large dataset from global environment
rm(Tmin_farm,Tmax_farm, amplitude_fh, amplitude_sh,fh,sh, Tmax_t,Tmin_t,Tmin_t1,Tmean_fh,Tmean_sh)

# ==============================================================================================================
# 4. Restructuring data
# ==============================================================================================================

# The following steps facilitate data handling in insurance calibration

# -----------------------------------------
# Winter wheat
# -----------------------------------------

# Wheat yields
wheat_yield_melted           <- melt(cbind(farm_yields[[1]][,-c(1:2)],id=1:nrow(farm_yields[[1]])),id="id")
colnames(wheat_yield_melted) <- c("Farm","year","yield")
wheat_yield_melted$year      <- as.numeric(as.character(wheat_yield_melted$year))

# Temperature exposure wheat
temp_wheat_exposure1 <-rbindlist(lapply(farm_wheat_exposure, rbindlist),idcol=T)
temp_wheat_exposure2 <- melt(temp_wheat_exposure1, id=c(".id","date"))
colnames(temp_wheat_exposure2) <- c("Farm","Date","Hour", "Temperature")
temp_wheat_exposure2$year      <-year(temp_wheat_exposure2$Date)

# -----------------------------------------
# rapeseed
# -----------------------------------------

# rapeseed yields
rapeseed_yield_melted <-melt(cbind(farm_yields[[2]][,-c(1:2)],id=1:nrow(farm_yields[[2]])),id="id")
colnames(rapeseed_yield_melted) <-c("Farm","year","yield")
rapeseed_yield_melted$year <-as.numeric(as.character(rapeseed_yield_melted$year))

# Temperature exposure rapeseed
temp_rapeseed_exposure1 <-rbindlist(lapply(farm_rapeseed_exposure, rbindlist),idcol=T)
temp_rapeseed_exposure2 <- melt(temp_rapeseed_exposure1, id=c(".id","date"))
colnames(temp_rapeseed_exposure2) <- c("Farm","Date","Hour", "Temperature")
temp_rapeseed_exposure2$year      <-year(temp_rapeseed_exposure2$Date)

# ==============================================================================================================
# 5. Expected utility model
# ==============================================================================================================

# Coefficients of constant relative risk aversion
# 0.5 == slightly risk-averse, 2 == moderately risk-averse, 4 == extremely risk-averse
alpha <- 2

# Power utility function
utility_function <- function(alpha,revenue){
  
  if(alpha == 1){
    log(revenue)
  } else {
    ((revenue^(1 - alpha)) / (1 - alpha))
  }
}

# Inverse utility function (to derive certainty equivalents)
inverse_utility_function <- function(alpha,eu){
  
  if(alpha == 1){
    exp(eu)
  } else {
    (eu*(1-alpha))^(1/(1-alpha))
  }
}

# ==============================================================================================================
# 6. Expected utility model: Uninsured status (base-line scenario)
# ==============================================================================================================

# Step 1: Detrending yields 
# Step 2: Calculate expected utility
# Note that we report revenue in dt/ha and not on monetary terms.

# -----------------------------------------
# Winter wheat
# -----------------------------------------

lin_trend_wheat     <- rlm(yield ~ year , data = wheat_yield_melted, method="M")
quadr_trend__wheat  <- rlm(yield ~ year + I(year^2) , data=wheat_yield_melted, method="M")

# Use Akaike information criterion (AIC) to identify more appropriate trend
AIC(lin_trend_wheat)
AIC(quadr_trend__wheat)

# AIC favors quadratic trend for farm-individual winter wheat yields (lower value)
# Calculation of detrended yields using a quadratic trend
detr_wheat_yields <- as.data.frame(matrix(NA, nrow=nrow(farm_yields[[1]]), ncol=length(seq(1995,last_year,1))))
colnames(detr_wheat_yields) <- seq(1995,last_year,1)

for (i in 1:nrow(detr_wheat_yields)){  
  for (t in 1:length(seq(1995,last_year,1))){
    
    detr_wheat_yields[i,t] <- ((coef(quadr_trend__wheat)[1] +last_year * coef(quadr_trend__wheat)[2] + (last_year^2) * coef(quadr_trend__wheat)[3])) - ((coef(quadr_trend__wheat)[1] + as.numeric(colnames(detr_wheat_yields)[t]) * coef(quadr_trend__wheat)[2] + (as.numeric(colnames(detr_wheat_yields)[t]))^2 * coef(quadr_trend__wheat)[3])) +  farm_yields[[1]][i,colnames(detr_wheat_yields)[t]]
    
  }
}

# Utility of uninsured status
utility_wheat_uninsured <- as.data.frame(matrix(NA, nrow=nrow(farm_yields[[1]]), ncol = length(seq(1995,last_year,1))))
for (i in 1:nrow(farm_yields[[1]])){
  for (t in 1:length(seq(1995,last_year,1))){
    utility_wheat_uninsured[i,t] <- utility_function(alpha = alpha, revenue = detr_wheat_yields[i,t])
  }
}

# Expected utility statitsics
summary_EU_wheat_uninsured <- as.data.frame(matrix(NA, nrow=nrow(farm_yields[[1]]),ncol=4))
colnames(summary_EU_wheat_uninsured) <- c("EU", "Expected Revenue", "CE", "Risk Premium")

for (i in 1:nrow(farm_yields[[1]])){
  
  # Expected utility
  summary_EU_wheat_uninsured [i,1] <- mean(as.numeric(utility_wheat_uninsured[i,]), na.rm=T)
  
  # Expected Revenue
  summary_EU_wheat_uninsured [i,2] <- mean(as.numeric(detr_wheat_yields[i,]), na.rm=T)
  
  # Certainty equivalent
  summary_EU_wheat_uninsured [i,3] <- inverse_utility_function(alpha = alpha, eu = as.numeric(summary_EU_wheat_uninsured[i,1])) 
  
  # Risk Premium
  summary_EU_wheat_uninsured [i,4] <- summary_EU_wheat_uninsured [i,2] - summary_EU_wheat_uninsured [i,3] 
}

# -----------------------------------------
# rapeseed
# -----------------------------------------

lin_trend_rapeseed   <- rlm(yield ~ year , data = rapeseed_yield_melted, method="M")
quadr_trend_rapeseed <- rlm(yield ~ year + I(year^2) , data=rapeseed_yield_melted, method="M")

# Use Akaike information criterion (AIC) to identify more appropriate trend
AIC(lin_trend_rapeseed)
AIC(quadr_trend_rapeseed)

# AIC favors quadratic trend for farm-individual rapeseed yields (lower value)
# Calculation of detrended yields using a quadratic trend

detr_rapeseed_yields <- as.data.frame(matrix(NA, nrow=nrow(farm_yields[[2]]), ncol=length(seq(1995,last_year,1))))
colnames(detr_rapeseed_yields) <- seq(1995,last_year,1)

for (i in 1:nrow(detr_rapeseed_yields)){  
  for (t in 1:length(seq(1995,last_year,1))){
    
    detr_rapeseed_yields[i,t] <- detr_rapeseed_yields[i,t] <- ((coef(quadr_trend_rapeseed)[1] +last_year * coef(quadr_trend_rapeseed)[2] + (last_year^2) * coef(quadr_trend_rapeseed)[3]))  - ((coef(quadr_trend_rapeseed)[1] + as.numeric(colnames(detr_rapeseed_yields)[t]) * coef(quadr_trend_rapeseed)[2] + (as.numeric(colnames(detr_rapeseed_yields)[t]))^2 * coef(quadr_trend_rapeseed)[3])) +  farm_yields[[2]][i,colnames(detr_rapeseed_yields)[t]]
    
  }
}

# Utility of uninsured status
utility_rapeseed_uninsured <- as.data.frame(matrix(NA, nrow=nrow(farm_yields[[2]]), ncol = length(seq(1995,last_year,1))))
for (i in 1:nrow(farm_yields[[2]])){
  for (t in 1:length(seq(1995,last_year,1))){
    utility_rapeseed_uninsured[i,t] <- utility_function(alpha = alpha, revenue = detr_rapeseed_yields[i,t])
  }
}

# Expected utility statitsics
summary_EU_rapeseed_uninsured <- as.data.frame(matrix(NA, nrow=nrow(farm_yields[[2]]),ncol=4))
colnames(summary_EU_rapeseed_uninsured) <- c("EU", "Expected Revenue", "CE", "Risk Premium")

for (i in 1:nrow(farm_yields[[2]])){
  
  # Expected utility
  summary_EU_rapeseed_uninsured [i,1] <- mean(as.numeric(utility_rapeseed_uninsured[i,]), na.rm=T)
  
  # Expected Revenue
  summary_EU_rapeseed_uninsured [i,2] <- mean(as.numeric(detr_rapeseed_yields[i,]), na.rm=T)
  
  # Certainty equivalent
  summary_EU_rapeseed_uninsured [i,3] <- inverse_utility_function(alpha = alpha, eu = as.numeric(summary_EU_rapeseed_uninsured[i,1])) 
  
  # Risk Premium
  summary_EU_rapeseed_uninsured [i,4] <- summary_EU_rapeseed_uninsured [i,2] - summary_EU_rapeseed_uninsured [i,3] 
}

# ==============================================================================================================
# 7. Winter wheat: out-of-sample calibration and testing for three knots
# ==============================================================================================================

# -------------------------------------------
# Definition of knots (out-of-sample)
# -------------------------------------------

# List with knots for each farm and model: [[i]] == farm; [knot,model]
list_knot_locations_kn3_wheat <- list()

for (i in 1:nrow(farm_yields[[1]])){
  
  # Leave-out data from farm i
  temp_sample <- temp_wheat_exposure2[which(temp_wheat_exposure2$Farm != i),]
  temp_knot_df <- as.data.frame(matrix(NA, nrow=3, ncol=3))
  row.names(temp_knot_df) <- c("knot1", "knot2", "knot3")
  colnames(temp_knot_df)  <- c("Best fit","Equal", "Quantile") 
  
  # Model 1: Best fit
  
  # lower and upper bound define lowest and highest knot location
  lower_bound <- ceiling(quantile(temp_sample$Temperature,0.05, type=1))
  upper_bound <- floor(quantile(temp_sample$Temperature,0.95, type=1))
  temp_range <- seq(lower_bound, upper_bound,1) 
  kn3_combi <- combn(temp_range,3, simplify = T)
  
  # Minimum space between 2 knots
  required_space <- 5
  
  # Get all knot combinations
  differences <- matrix(NA, nrow=3, ncol=ncol(kn3_combi))
  differences[1,] <- abs(kn3_combi[1,] - kn3_combi[2,])
  differences[2,] <- abs(kn3_combi[1,] - kn3_combi[3,])
  differences[3,] <- abs(kn3_combi[2,] - kn3_combi[3,])
  knots_3_combi_good <- kn3_combi[,which(differences[1,] >= required_space & differences[2,] >= required_space & differences[3,] >= required_space)]
  
  # vector containing residual sum of squares (RSS)
  RSS_3knots <- vector(length=ncol(knots_3_combi_good ))
  
  # Get the RSS for each model
  for (c in 1:ncol(knots_3_combi_good)){
 
  # Get new time series and aggregate hourly values  
  temp1 <- as.matrix(rcspline.eval(temp_sample$Temperature,knots= knots_3_combi_good[,c], inclx=T))
  temp2 <- cbind(temp_sample,temp1)
  
  temp3<-temp2 %>%
    group_by(Farm,year)%>%
    summarise_at(vars(x,V2),sum)
  
  # Match yearly values with yearly yields
  sub_wheat_yield_melted <- wheat_yield_melted[which(wheat_yield_melted$Farm != i),]
  temp4 <- join(sub_wheat_yield_melted,temp3, by=c("Farm", "year"))
  
  # Run regression
  temp_reg <- lm(yield ~ x + V2 + year + I(year^2) + as.factor(Farm)-1, data = temp4)
  RSS_3knots[c] <- sum(resid(temp_reg)^2)
  
  rm(temp_reg,temp4,sub_wheat_yield_melted,temp3,temp2,temp1)
  
  print(paste(c/ncol(knots_3_combi_good),"for farm i", i/nrow(farm_yields[[1]])))
  
  }
  
  # Best number of knots
  temp_knot_df[,1] <- knots_3_combi_good[,which.min(RSS_3knots)]
  rm(RSS_3knots, knots_3_combi_good, differences, kn3_combi,temp_range,upper_bound, lower_bound)
  
  
  # Model 2: equally spaced
  temp_knot_df[,2] <- c(min(temp_sample$Temperature) + ((max(temp_sample$Temperature)-min(temp_sample$Temperature)) / 4),
                        min(temp_sample$Temperature) + 2* ((max(temp_sample$Temperature)-min(temp_sample$Temperature)) / 4),
                        min(temp_sample$Temperature) + 3* ((max(temp_sample$Temperature)-min(temp_sample$Temperature)) / 4))
  
  
  
  # Model 3: 10%, 50%, 90% quantile
  temp_knot_df[,3] <- c(quantile(temp_sample$Temperature,0.1, type=1),
                                     quantile(temp_sample$Temperature,0.5, type=1),
                                     quantile(temp_sample$Temperature,0.9, type=1))
  
  
  
  list_knot_locations_kn3_wheat[[i]] <-  temp_knot_df
  rm(temp_sample, temp_knot_df)
  print(i / nrow(farm_yields[[1]]))
}

# Testing whether model with best fit is better than linear model
# Leave-out data from farm i

# 1 = cubic spline is superior to linear model
comparison_vector <- vector()

for (i in 1:nrow(farm_yields[[1]])){
  temp_sample <- temp_wheat_exposure2[which(temp_wheat_exposure2$Farm != i),]
  
  # The RCS model
  # Get new time series and aggregate hourly values  
  temp1 <- as.matrix(rcspline.eval(temp_sample$Temperature,knots= list_knot_locations_kn3_wheat[[i]][,1], inclx=T))
  temp2 <- cbind(temp_sample,temp1)
  
  temp3<-temp2 %>%
    group_by(Farm,year)%>%
    summarise_at(vars(x,V2),sum)
  
  # Match yearly values with yearly yields
  sub_wheat_yield_melted <- wheat_yield_melted[which(wheat_yield_melted$Farm != i),]
  temp4 <- join(sub_wheat_yield_melted,temp3, by=c("Farm", "year"))
  
  # Run cubic model (restricted cubic spline model)
  temp_cub_reg <- lm(yield ~ x + V2 + year + I(year^2) + as.factor(Farm)-1, data = temp4)
  
  # Run linear model
  temp_lin_reg <- lm(yield ~ x + year + I(year^2) + as.factor(Farm)-1, data = temp4)
  
  AIC_lin <- AIC(temp_lin_reg)
  AIC_cub <- AIC(temp_cub_reg)
  
  if (AIC_cub < AIC_lin){
    comparison_vector[i] <- 1
  } else {comparison_vector[i] <- 0}
  
  rm(temp_sample, temp1, temp2, temp3,sub_wheat_yield_melted,temp_cub_reg, temp_lin_reg,AIC_lin, AIC_cub)
  print(i / nrow(farm_yields[[1]]))
}

# Check 
comparison_vector
rm(comparison_vector)


# -------------------------------------------
# Out-of-sample calibration: leave-out-farm i
# -------------------------------------------

# Yields together with new hourly temperature time series: [[i]] == farm; [[m]] == new time series for knot specification m
list_wheat_matrix_modelkn3_aggregated <- list() # [[i]][[m]]

# New time series of hourly temperature exposures: [[i]] == farm; [[m]] == new time series for knot specification m
list_wheat_matrix_modelkn3_hourly <- list() # [[i]][[m]]

for (i in 1:nrow(farm_yields[[1]])){
  temp_hourly <- list() #[[m]]
  temp_aggregated <- list() # [[m]]
  
  for (m in 1:ncol(list_knot_locations_kn3_wheat[[i]])){
    
    # Get new time series for RCS   
    temp1 <- as.matrix(rcspline.eval(temp_wheat_exposure2$Temperature,knots= list_knot_locations_kn3_wheat[[i]][,m], inclx=T))
    temp2 <- cbind(temp_wheat_exposure2, temp1)
    
    # Aggregate to yearly values
    temp3<-temp2 %>%
      group_by(Farm,year)%>%
      summarise_at(vars(x,V2),sum)
    
    temp_hourly[[m]]     <- temp2
    temp_aggregated[[m]] <- join(wheat_yield_melted,temp3, by=c("Farm", "year"))
    
    rm(temp1,temp2,temp3)
  }
list_wheat_matrix_modelkn3_aggregated[[i]] <- temp_aggregated
list_wheat_matrix_modelkn3_hourly[[i]] <- temp_hourly

rm(temp_aggregated, temp_hourly)
}

# Now, we have all data together (yields and temperature time series)
# We continue with out-of-sample regression. See the paper for more details.

# List with model outputs. [[i]] == farm i; [[m]] == knot specification m
list_farm_wheat_panelmodelkn3_modeloutputs <- list()

for (i in 1:nrow(farm_yields[[1]])){
  
  # Contains the m models of farm i.
  list_farm_wheat_panelmodelkn3_models.of.farm <- list() 
  for (m in 1:ncol(list_knot_locations_kn3_wheat[[i]])){
    
    # Aggregated data without farm i (leave-out-farm i for out-of-sample training)
    temp1 <- subset(list_wheat_matrix_modelkn3_aggregated[[i]][[m]], list_wheat_matrix_modelkn3_aggregated[[i]][[m]] [,"Farm"] != i)
    
    # Run panel regression
    temp_reg <- lm(yield ~ x + V2 + year + I(year^2) + as.factor(Farm)-1, data = temp1)
    list_farm_wheat_panelmodelkn3_models.of.farm[[m]] <- temp_reg
    
    rm(temp1,temp_reg)
    
  }
  
  list_farm_wheat_panelmodelkn3_modeloutputs [[i]] <- list_farm_wheat_panelmodelkn3_models.of.farm   
  rm(list_farm_wheat_panelmodelkn3_models.of.farm)
}


# Calculation of hourly effects (marginal temperature effect estimated at farm i)
list_farm_wheat_hourly_effect_kn3 <- list() # [[farm]] [[model]][]

for (i in 1:nrow(farm_yields[[1]])){
  
  temp_list <- list()
  for (m in 1:ncol(list_knot_locations_kn3_wheat[[i]])){
    
    # Subset of farm i
    temp1 <- which(list_wheat_matrix_modelkn3_hourly[[i]][[m]][,"Farm"] == i) 
    temp2 <- list_wheat_matrix_modelkn3_hourly[[i]][[m]][temp1,]  
    
    # Calculate temperature effect
    temp2$effect <- as.matrix(temp2[,6:7]) %*% coef(list_farm_wheat_panelmodelkn3_modeloutputs [[i]][[m]])[1:2]
    
    temp_list[[m]] <- temp2  
    rm(temp1,temp2)  
  }
  list_farm_wheat_hourly_effect_kn3[[i]] <- temp_list
  rm(temp_list)
}

# -------------------------------------------
# Historical payouts & premium
# -------------------------------------------

# Different strike level temperatures
strike_temperature <- c(seq(13,36,by=1))

# Historical payouts with following structure: [[i]] == of farm i, [[m]] with model specification m, 
# [[m]] contains a data.frame with [strike level, year]
list_farm_wheat_payouts_kn3 <- list() 
for (i in 1:nrow(farm_yields[[1]])){
  
  # temp list for specification[[m]]
  temp_list <- list() 
  for (m in 1:ncol(list_knot_locations_kn3_wheat[[i]])){
    
    temp_model <- list_farm_wheat_hourly_effect_kn3 [[i]] [[m]]
    
    # DF for payouts
    temp_df <- as.data.frame(matrix(NA, nrow=length(strike_temperature), ncol=length(seq(1995,last_year,1))))
    row.names(temp_df) <- strike_temperature 
    colnames(temp_df) <- seq(1995,last_year,1)
    
    for (t in 1:ncol(temp_df)){
      # Subset of year
      sub_model_df <- subset(temp_model, temp_model$year == as.numeric(colnames(temp_df)[t]))
      # Subset only negative expected impacts trigger a payout 
      sub2_model_df <- subset(sub_model_df, sub_model_df$effect < 0)
      
      for (strike in 1:length(strike_temperature)){
        
        # Subset only temperature above the strike level trigger a payout
        sub3_model_df <- subset(sub2_model_df, sub2_model_df$Temperature >= strike_temperature[strike]) 
        
        # Calculate the cumulative payout
        temp_df[strike,t] <- sum(sub3_model_df$effect) * (-1)
        rm(sub3_model_df)
      }
      
      rm(sub2_model_df, sub_model_df)  
    }
    temp_list[[m]] <- temp_df
    rm(temp_df,temp_model)
  }
  list_farm_wheat_payouts_kn3[[i]] <- temp_list
  rm(temp_list)
} 

# Calculate the premium
# Structure: [farm i, model m, strike level]
array_premium_wheat_kn3 <- array(dim=c(nrow(farm_yields[[1]]), ncol(list_knot_locations_kn3_wheat[[i]]),length(strike_temperature)))
dimnames(array_premium_wheat_kn3)[[3]] <- strike_temperature

for (i in 1:nrow(farm_yields[[1]])){
  for (m in 1:ncol(list_knot_locations_kn3_wheat[[i]])){
    for (strike in 1:length(strike_temperature)){
      temp1 <- which(!is.na(detr_wheat_yields[i,]))
      array_premium_wheat_kn3[i,m,strike] <- mean(as.numeric(list_farm_wheat_payouts_kn3[[i]][[m]][strike, temp1]), na.rm=T)
      rm(temp1)
    }  
  }
}

# -------------------------------------------
# Insured revenue and expected utility
# -------------------------------------------

# Calculation of revenue [[strike]][m,i,t]
list_farm_wheat_wealth_modelskn3_strike <- list()

for (strike in 1:length(strike_temperature)){
  
  array_farm_wheat_wealth_modelskn3 <- array(dim=c(ncol(list_knot_locations_kn3_wheat[[1]]), nrow(farm_yields[[1]]), length(seq(1995,last_year,1))))
  dimnames(array_farm_wheat_wealth_modelskn3) [[3]] <-seq(1995,last_year,1) 
  
  for (m in 1:ncol(list_knot_locations_kn3_wheat[[1]])){
    for (i in 1:nrow(farm_yields[[1]])){
      for (t in 1: length(seq(1995,last_year,1))){
        
        array_farm_wheat_wealth_modelskn3 [m,i,t] <- as.numeric(detr_wheat_yields[i,t]) + as.numeric(list_farm_wheat_payouts_kn3 [[i]][[m]] [strike,t]) - as.numeric(array_premium_wheat_kn3[i,m,strike]) 
        
      }
    }
  }
  
  list_farm_wheat_wealth_modelskn3_strike [[strike]] <- array_farm_wheat_wealth_modelskn3 
  rm(array_farm_wheat_wealth_modelskn3)
}

# Utility of being insured with specification m.
# List has structure: [[strike]][mi,i,t]
list_farm_wheat_wealth_modelkn3_strike_utility <- list()

for (strike in 1:length(strike_temperature)){
  
  temp_array <- array(dim=c(ncol(list_knot_locations_kn3_wheat[[i]]), nrow(farm_yields[[1]]), length(seq(1995,last_year,1))))
  dimnames(temp_array) [[3]] <-seq(1995,last_year,1) 
  
  for (m in 1:ncol(list_knot_locations_kn3_wheat[[i]])){
    for (i in 1:nrow(farm_yields[[1]])){
      for (t in 1: length(seq(1995,last_year,1))){
        
        temp_array [m,i,t] <- utility_function(alpha = alpha, revenue = list_farm_wheat_wealth_modelskn3_strike [[strike]][m,i,t])
        
      }
    }
  }
  
  list_farm_wheat_wealth_modelkn3_strike_utility  [[strike]] <- temp_array 
  rm(temp_array)
}

# Expected utility statistics of being insured
# Structure of list: 
summary_EU_wheat_insured_kn3 <- list() # [[strike]] [m,i,c(EU,E,CE,RP)]

for (strike in 1:length(strike_temperature)){
  temp_array <- array(dim=c(ncol(list_knot_locations_kn3_wheat[[1]]), nrow(farm_yields[[1]]), 4))
  dimnames(temp_array) [[3]] <- c("EU", "Expected Revenue", "CE", "RP")
  
  for (m in 1:ncol(list_knot_locations_kn3_wheat[[1]])){
    for (i in 1:nrow(farm_yields[[1]])){
      
      # Expected utility
      temp_array  [m,i,1] <- mean(as.numeric(list_farm_wheat_wealth_modelkn3_strike_utility[[strike]][m,i,]), na.rm=T)
      
      # Expected Revenue
      temp_array  [m,i,2] <- mean(as.numeric(list_farm_wheat_wealth_modelskn3_strike[[strike]][m,i,]), na.rm=T)
      
      # Certainty Equivalent
      temp_array  [m,i,3] <- inverse_utility_function(alpha = alpha, eu = as.numeric(temp_array[m,i,1]))
      
      # Risk Premium
      temp_array  [m,i,4] <- temp_array  [m,i,2] -  temp_array  [m,i,3]
    }
  }
  summary_EU_wheat_insured_kn3 [[strike]] <- temp_array
  rm(temp_array)
}

# -------------------------------------------
# Testing
# -------------------------------------------

average_rel_change_wheat_knot3_uninsured <- as.data.frame(matrix(nrow=ncol(list_knot_locations_kn3_wheat[[1]]), ncol=length(strike_temperature)))
colnames(average_rel_change_wheat_knot3_uninsured) <- strike_temperature

for (m in 1:ncol(list_knot_locations_kn3_wheat[[1]])){
  for (strike in 1:length(strike_temperature)){
    
    average_rel_change_wheat_knot3_uninsured [m,strike] <- mean((summary_EU_wheat_insured_kn3[[strike]][m,,4] - summary_EU_wheat_uninsured[,4]) / summary_EU_wheat_uninsured[,4])
    
  }
}

average_rel_change_wheat_knot3_uninsured_test <- as.data.frame(matrix(nrow=ncol(list_knot_locations_kn3_wheat[[1]]), ncol=length(strike_temperature)))
colnames(average_rel_change_wheat_knot3_uninsured_test) <- strike_temperature

for (m in 1:ncol(list_knot_locations_kn3_wheat[[1]])){
  for (strike in 1:length(strike_temperature)){
    
    average_rel_change_wheat_knot3_uninsured_test [m,strike] <- wilcox.test(summary_EU_wheat_insured_kn3[[strike]][m,,4], summary_EU_wheat_uninsured[,4], paired=T, alternative = "l")$p.value
    
  }
}

# ==============================================================================================================
# 8. Rapeseed: out-of-sample calibration and testing for three knots
# ==============================================================================================================

# -------------------------------------------
# Definition of knots (out-of-sample)
# -------------------------------------------

# List with knots for each farm and model: [[i]] == farm; [knot,model]
list_knot_locations_kn3_rapeseed <- list()

for (i in 1:nrow(farm_yields[[2]])){
  
  # Leave-out data from farm i
  temp_sample <- temp_rapeseed_exposure2[which(temp_rapeseed_exposure2$Farm != i),]
  temp_knot_df <- as.data.frame(matrix(NA, nrow=3, ncol=3))
  row.names(temp_knot_df) <- c("knot1", "knot2", "knot3")
  colnames(temp_knot_df)  <- c("Best fit","Equal", "Quantile") 
  
  # Model 1: Best fit
  
  # lower and upper bound define lowest and highest knot location
  lower_bound <- ceiling(quantile(temp_sample$Temperature,0.05, type=1))
  upper_bound <- floor(quantile(temp_sample$Temperature,0.95, type=1))
  temp_range <- seq(lower_bound, upper_bound,1) 
  kn3_combi <- combn(temp_range,3, simplify = T)
  
  # Minimum space between 2 knots in degree-Celsius
  required_space <- 5
  
  # Get all knot combinations
  differences <- matrix(NA, nrow=3, ncol=ncol(kn3_combi))
  differences[1,] <- abs(kn3_combi[1,] - kn3_combi[2,])
  differences[2,] <- abs(kn3_combi[1,] - kn3_combi[3,])
  differences[3,] <- abs(kn3_combi[2,] - kn3_combi[3,])
  knots_3_combi_good <- kn3_combi[,which(differences[1,] >= required_space & differences[2,] >= required_space & differences[3,] >= required_space)]
  
  # vector containing residual sum of squares (RSS)
  RSS_3knots <- vector(length=ncol(knots_3_combi_good ))
  
  # Get the RSS for each model
  for (c in 1:ncol(knots_3_combi_good)){
    
    # Get new time series and aggregate hourly values  
    temp1 <- as.matrix(rcspline.eval(temp_sample$Temperature,knots= knots_3_combi_good[,c], inclx=T))
    temp2 <- cbind(temp_sample,temp1)
    
    temp3<-temp2 %>%
      group_by(Farm,year)%>%
      summarise_at(vars(x,V2),sum)
    
    # Match yearly values with yearly yields
    sub_rapeseed_yield_melted <- rapeseed_yield_melted[which(rapeseed_yield_melted$Farm != i),]
    temp4 <- join(sub_rapeseed_yield_melted,temp3, by=c("Farm", "year"))
    
    # Run regression
    temp_reg <- lm(yield ~ x + V2 + year + I(year^2) + as.factor(Farm)-1, data = temp4)
    RSS_3knots[c] <- sum(resid(temp_reg)^2)
    
    rm(temp_reg,temp4,sub_rapeseed_yield_melted,temp3,temp2,temp1)
    
    print(paste(c/ncol(knots_3_combi_good),"for farm i", i/nrow(farm_yields[[2]])))
    
  }
  
  # Best number of knots
  temp_knot_df[,1] <- knots_3_combi_good[,which.min(RSS_3knots)]
  rm(RSS_3knots, knots_3_combi_good, differences, kn3_combi,temp_range,upper_bound, lower_bound)
  
  
  # Model 2: equally spaced
  temp_knot_df[,2] <- c(min(temp_sample$Temperature) + ((max(temp_sample$Temperature)-min(temp_sample$Temperature)) / 4),
                        min(temp_sample$Temperature) + 2* ((max(temp_sample$Temperature)-min(temp_sample$Temperature)) / 4),
                        min(temp_sample$Temperature) + 3* ((max(temp_sample$Temperature)-min(temp_sample$Temperature)) / 4))
  
  
  
  # Model 3: 10%, 50%, 90% quantile
  temp_knot_df[,3] <- c(quantile(temp_sample$Temperature,0.1, type=1),
                        quantile(temp_sample$Temperature,0.5, type=1),
                        quantile(temp_sample$Temperature,0.9, type=1))
  
  
  
  list_knot_locations_kn3_rapeseed[[i]] <-  temp_knot_df
  rm(temp_sample, temp_knot_df)
  print(i / nrow(farm_yields[[2]]))
}


# Testing whether model with best fit is better than linear model
# Leave-out data from farm i

# 1 = cubic spline is superior to linear model
comparison_vector <- vector()

for (i in 1:nrow(farm_yields[[2]])){
  temp_sample <- temp_rapeseed_exposure2[which(temp_rapeseed_exposure2$Farm != i),]
  
  # The RCS model
  # Get new time series and aggregate hourly values  
  temp1 <- as.matrix(rcspline.eval(temp_sample$Temperature,knots= list_knot_locations_kn3_rapeseed[[i]][,1], inclx=T))
  temp2 <- cbind(temp_sample,temp1)
  
  temp3<-temp2 %>%
    group_by(Farm,year)%>%
    summarise_at(vars(x,V2),sum)
  
  # Match yearly values with yearly yields
  sub_rapeseed_yield_melted <- rapeseed_yield_melted[which(rapeseed_yield_melted$Farm != i),]
  temp4 <- join(sub_rapeseed_yield_melted,temp3, by=c("Farm", "year"))
  
  # Run cubic model (restricted cubic spline model)
  temp_cub_reg <- lm(yield ~ x + V2 + year + I(year^2) + as.factor(Farm)-1, data = temp4)
  
  # Run linear model
  temp_lin_reg <- lm(yield ~ x + year + I(year^2) + as.factor(Farm)-1, data = temp4)
  
  AIC_lin <- AIC(temp_lin_reg)
  AIC_cub <- AIC(temp_cub_reg)
  
  if (AIC_cub < AIC_lin){
    comparison_vector[i] <- 1
  } else {comparison_vector[i] <- 0}
  
  rm(temp_sample, temp1, temp2, temp3,sub_rapeseed_yield_melted,temp_cub_reg, temp_lin_reg,AIC_lin, AIC_cub)
  print(i / nrow(farm_yields[[2]]))
}

# Check 
comparison_vector
rm(comparison_vector)

# -------------------------------------------
# Out-of-sample calibration: leave-out-farm i
# -------------------------------------------

# Yields together with new hourly temperature time series: [[i]] == farm, [[m]] == new time series for knot specification m
list_rapeseed_matrix_modelkn3_aggregated <- list() 

# New time series of hourly temperature exposures: [[i]] [[m]] == new time series for knot specification m
list_rapeseed_matrix_modelkn3_hourly <- list()
memory.size(max=T)


for (i in 1:nrow(farm_yields[[2]])){
  temp_hourly <- list() #[[m]]
  temp_aggregated <- list() #[[m]]
  
  for (m in 1:ncol(list_knot_locations_kn3_rapeseed[[i]])){
    
    # Get new time series for RCS   
    temp1 <- as.matrix(rcspline.eval(temp_rapeseed_exposure2$Temperature,knots= list_knot_locations_kn3_rapeseed[[i]][,m], inclx=T))
    temp2 <- cbind(temp_rapeseed_exposure2, temp1)
    
    # Aggregate to yearly values
    temp3<-temp2 %>%
      group_by(Farm,year)%>%
      summarise_at(vars(x,V2),sum)
    
    temp_hourly[[m]] <- temp2
    temp_aggregated[[m]] <- join(rapeseed_yield_melted,temp3, by=c("Farm","year"))
    
    rm(temp1,temp2,temp3)
  }
  
list_rapeseed_matrix_modelkn3_hourly[[i]]     <- temp_hourly
list_rapeseed_matrix_modelkn3_aggregated[[i]] <- temp_aggregated
    
rm(temp_hourly, temp_aggregated)   
  
}

# Now, we have all data together (yields and temperature time series)
# We continue with regression. See the paper for more details.

# List with model outputs. [[i]] == farm i; [[m]] == knot specification m
list_farm_rapeseed_panelmodelkn3_modeloutputs <- list()

for (i in 1:nrow(farm_yields[[2]])){
  
  # Contains the m models of farm i.
  list_farm_rapeseed_panelmodelkn3_models.of.farm <- list() 
  for (m in 1:ncol(list_knot_locations_kn3_rapeseed[[2]])){
    
    # Aggregated data without farm i (leave-out-farm i for out-of-sample training)
    temp1 <- subset(list_rapeseed_matrix_modelkn3_aggregated[[i]][[m]], list_rapeseed_matrix_modelkn3_aggregated[[i]][[m]] [,"Farm"] != i)
    
    # Run panel regression
    temp_reg <- lm(yield ~ x + V2 + year + I(year^2) + as.factor(Farm)-1, data = temp1)
    list_farm_rapeseed_panelmodelkn3_models.of.farm[[m]] <- temp_reg
    
    rm(temp1,temp_reg)
    
  }
  
  list_farm_rapeseed_panelmodelkn3_modeloutputs [[i]] <- list_farm_rapeseed_panelmodelkn3_models.of.farm   
  rm(list_farm_rapeseed_panelmodelkn3_models.of.farm)
}


# Calculation of hourly effects (marginal temperature effect)
list_farm_rapeseed_hourly_effect_kn3 <- list() # [[farm]] [[model]][]

for (i in 1:nrow(farm_yields[[2]])){
  
  temp_list <- list()
  for (m in 1:ncol(list_knot_locations_kn3_rapeseed[[i]])){
    
    # Subset of farm i
    temp1 <- which(list_rapeseed_matrix_modelkn3_hourly[[i]][[m]][,"Farm"] == i) 
    temp2 <- list_rapeseed_matrix_modelkn3_hourly[[i]][[m]][temp1,]  
    
    # Calculate temperature effect
    temp2$effect <- as.matrix(temp2[,6:7]) %*% coef(list_farm_rapeseed_panelmodelkn3_modeloutputs [[i]][[m]])[1:2]
    
    temp_list[[m]] <- temp2  
    rm(temp1,temp2)  
  }
  list_farm_rapeseed_hourly_effect_kn3[[i]] <- temp_list
  rm(temp_list)
}

# -------------------------------------------
# Historical payouts & premium
# -------------------------------------------

# Different strike level temperatures
strike_temperature <- c(seq(13,36,by=1))

# Historical payouts with following structure: [[i]] == of farm i, [[m]] with model specification m, 
# [[m]] contains a data.frame with [strike level, year]
list_farm_rapeseed_payouts_kn3 <- list() 

for (i in 1:nrow(farm_yields[[2]])){
  
  # temp list for specification[[m]]
  temp_list <- list() 
  for (m in 1:ncol(list_knot_locations_kn3_rapeseed[[i]])){
    
    temp_model <- list_farm_rapeseed_hourly_effect_kn3 [[i]] [[m]]
    
    # DF for payouts
    temp_df <- as.data.frame(matrix(NA, nrow=length(strike_temperature), ncol=length(seq(1995,last_year,1))))
    row.names(temp_df) <- strike_temperature 
    colnames(temp_df) <- seq(1995,last_year,1)
    
    for (t in 1:ncol(temp_df)){
      # Subset of year
      sub_model_df <- subset(temp_model, temp_model$year == as.numeric(colnames(temp_df)[t]))
      # Subset only negative expected impacts trigger a payout 
      sub2_model_df <- subset(sub_model_df, sub_model_df$effect < 0)
      
      for (strike in 1:length(strike_temperature)){
        
        # Subset only temperature above the strike level trigger a payout
        sub3_model_df <- subset(sub2_model_df, sub2_model_df$Temperature >= strike_temperature[strike]) 
        
        # Calculate the cumulative payout.
        # -1 because negative yield response means yield reduction for which the payout formula compensates.
        temp_df[strike,t] <- sum(sub3_model_df$effect) * (-1)
        rm(sub3_model_df)
      }
      
      rm(sub2_model_df, sub_model_df)  
    }
    temp_list[[m]] <- temp_df
    rm(temp_df,temp_model)
  }
  list_farm_rapeseed_payouts_kn3[[i]] <- temp_list
  rm(temp_list)
} 

# Calculate the premium
# Structure: [farm i, model m, strike level]
array_premium_rapeseed_kn3 <- array(dim=c(nrow(farm_yields[[2]]), ncol(list_knot_locations_kn3_rapeseed[[1]]),length(strike_temperature)))
dimnames(array_premium_rapeseed_kn3)[[3]] <- strike_temperature

for (i in 1:nrow(farm_yields[[2]])){
  for (m in 1:ncol(list_knot_locations_kn3_rapeseed[[i]])){
    for (strike in 1:length(strike_temperature)){
      temp1 <- which(!is.na(detr_rapeseed_yields[i,]))
      array_premium_rapeseed_kn3[i,m,strike] <- mean(as.numeric(list_farm_rapeseed_payouts_kn3[[i]][[m]][strike, temp1]), na.rm=T)
      rm(temp1)
    }  
  }
}

# -------------------------------------------
# Insured revenue and expected utility
# -------------------------------------------

# Calculation of revenue [[strike]][m,i,t]
list_farm_rapeseed_wealth_modelskn3_strike <- list()

for (strike in 1:length(strike_temperature)){
  
  array_farm_rapeseed_wealth_modelskn3 <- array(dim=c(ncol(list_knot_locations_kn3_rapeseed[[1]]), nrow(farm_yields[[2]]), length(seq(1995,last_year,1))))
  dimnames(array_farm_rapeseed_wealth_modelskn3) [[3]] <-seq(1995,last_year,1) 
  
  for (m in 1:ncol(list_knot_locations_kn3_rapeseed[[1]])){
    for (i in 1:nrow(farm_yields[[2]])){
      for (t in 1: length(seq(1995,last_year,1))){
        
        array_farm_rapeseed_wealth_modelskn3 [m,i,t] <- as.numeric(detr_rapeseed_yields[i,t]) + as.numeric(list_farm_rapeseed_payouts_kn3 [[i]][[m]] [strike,t]) - as.numeric(array_premium_rapeseed_kn3[i,m,strike]) 
        
      }
    }
  }
  
  list_farm_rapeseed_wealth_modelskn3_strike [[strike]] <- array_farm_rapeseed_wealth_modelskn3 
  rm(array_farm_rapeseed_wealth_modelskn3)
}

# Utility of being insured with specification m.
# List has structure: [[strike]][mi,i,t]
list_farm_rapeseed_wealth_modelkn3_strike_utility <- list()

for (strike in 1:length(strike_temperature)){
  
  temp_array <- array(dim=c(ncol(list_knot_locations_kn3_rapeseed[[1]]), nrow(farm_yields[[2]]), length(seq(1995,last_year,1))))
  dimnames(temp_array) [[3]] <-seq(1995,last_year,1) 
  
  for (m in 1:ncol(list_knot_locations_kn3_rapeseed[[1]])){
    for (i in 1:nrow(farm_yields[[2]])){
      for (t in 1: length(seq(1995,last_year,1))){
        
        temp_array [m,i,t] <- utility_function(alpha = alpha, revenue = list_farm_rapeseed_wealth_modelskn3_strike [[strike]][m,i,t])
        
      }
    }
  }
  
  list_farm_rapeseed_wealth_modelkn3_strike_utility  [[strike]] <- temp_array 
  rm(temp_array)
}

# Expected utility statistics of being insured
# Structure of list: 
summary_EU_rapeseed_insured_kn3 <- list() # [[strike]] [m,i,c(EU,E,CE,RP)]

for (strike in 1:length(strike_temperature)){
  temp_array <- array(dim=c(ncol(list_knot_locations_kn3_rapeseed[[1]]), nrow(farm_yields[[2]]), 4))
  dimnames(temp_array) [[3]] <- c("EU", "Expected Revenue", "CE", "RP")
  
  for (m in 1:ncol(list_knot_locations_kn3_rapeseed[[1]])){
    for (i in 1:nrow(farm_yields[[2]])){
      
      # Expected utility
      temp_array  [m,i,1] <- mean(as.numeric(list_farm_rapeseed_wealth_modelkn3_strike_utility[[strike]][m,i,]), na.rm=T)
      
      # Expected Revenue
      temp_array  [m,i,2] <- mean(as.numeric(list_farm_rapeseed_wealth_modelskn3_strike[[strike]][m,i,]), na.rm=T)
      
      # Certainty Equivalent
      temp_array  [m,i,3] <- inverse_utility_function(alpha = alpha, eu = as.numeric(temp_array[m,i,1]))
      
      # Risk Premium
      temp_array  [m,i,4] <- temp_array  [m,i,2] -  temp_array  [m,i,3]
    }
  }
  summary_EU_rapeseed_insured_kn3 [[strike]] <- temp_array
  rm(temp_array)
}

# -------------------------------------------
# Testing
# -------------------------------------------

average_rel_change_rapeseed_knot3_uninsured <- as.data.frame(matrix(nrow=ncol(list_knot_locations_kn3_rapeseed[[1]]), ncol=length(strike_temperature)))
colnames(average_rel_change_rapeseed_knot3_uninsured) <- strike_temperature

for (m in 1:ncol(list_knot_locations_kn3_rapeseed[[1]])){
  for (strike in 1:length(strike_temperature)){
    
    average_rel_change_rapeseed_knot3_uninsured [m,strike] <- mean((summary_EU_rapeseed_insured_kn3[[strike]][m,,4] - summary_EU_rapeseed_uninsured[,4]) / summary_EU_rapeseed_uninsured[,4])
    
  }
}

average_rel_change_rapeseed_knot3_uninsured_test <- as.data.frame(matrix(nrow=ncol(list_knot_locations_kn3_rapeseed[[1]]), ncol=length(strike_temperature)))
colnames(average_rel_change_rapeseed_knot3_uninsured_test) <- strike_temperature

for (m in 1:ncol(list_knot_locations_kn3_rapeseed[[1]])){
  for (strike in 1:length(strike_temperature)){
    
    average_rel_change_rapeseed_knot3_uninsured_test [m,strike] <-wilcox.test(summary_EU_rapeseed_insured_kn3[[strike]][m,,4], summary_EU_rapeseed_uninsured[,4], paired=T, alternative="l")$p.value
    
  }
}

# ==============================================================================================================
# 9. Winter wheat: out-of-sample calibration and testing for five knots
# ==============================================================================================================

# -------------------------------------------
# Definition of knots
# -------------------------------------------

knot_locations_kn5_wheat            <- as.data.frame(matrix(NA, nrow=5, ncol=1))
row.names(knot_locations_kn5_wheat) <- c("knot1", "knot2", "knot3","knot4","knot5")
colnames(knot_locations_kn5_wheat)  <- c("Schlenker") 

# Equally spaced across temperature support
knot_locations_kn5_wheat [,1] <- c(5,10,15,20,25)

'
# 10%, 50% and 90% quantile
knot_locations_kn5_wheat [,2] <- c(quantile(temp_wheat_exposure2$Temperature,0.05, type=1),
                                   quantile(temp_wheat_exposure2$Temperature,0.275, type=1),
                                   quantile(temp_wheat_exposure2$Temperature,0.5, type=1),
                                   quantile(temp_wheat_exposure2$Temperature,0.725, type=1),
                                   quantile(temp_wheat_exposure2$Temperature,0.95, type=1))'

# -------------------------------------------
# Out-of-sample calibration: leave-out-farm i
# -------------------------------------------

# Yields together with new hourly temperature time series: [[m]] == new time series for knot specification m
list_wheat_matrix_modelkn5_aggregated <- list() # [[m]]

# New time series of hourly temperature exposures: [[m]] == new time series for knot specification m
list_wheat_matrix_modelkn5_hourly <- list() # [[m]]

for (m in 1:ncol(knot_locations_kn5_wheat)){
  
  # Get new time series for RCS   
  temp1 <- as.matrix(rcspline.eval(temp_wheat_exposure2$Temperature,knots= knot_locations_kn5_wheat[,m], inclx=T))
  temp2 <- cbind(temp_wheat_exposure2, temp1)
  
  # Aggregate to yearly values
  temp3<-temp2 %>%
    group_by(Farm,year)%>%
    summarise_at(vars(x,V2,V3,V4),sum)
  
  list_wheat_matrix_modelkn5_hourly[[m]]     <- temp2
  list_wheat_matrix_modelkn5_aggregated[[m]] <- join(wheat_yield_melted,temp3, by=c("Farm", "year"))
  
  rm(temp1,temp2,temp3)
}

# Now, we have all data together (yields and temperature time series)
# We continue with regression. See the paper for more details.

# List with model outputs. [[i]] == farm i; [[m]] == knot specification m
list_farm_wheat_panelmodelkn5_modeloutputs <- list()

for (i in 1:nrow(farm_yields[[1]])){
  
  # Contains the m models of farm i.
  list_farm_wheat_panelmodelkn5_models.of.farm <- list() 
  for (m in 1:ncol(knot_locations_kn5_wheat)){
    
    # Aggregated data without farm i (leave-out-farm i for out-of-sample training)
    temp1 <- subset(list_wheat_matrix_modelkn5_aggregated[[m]], list_wheat_matrix_modelkn5_aggregated[[m]] [,"Farm"] != i)
    
    # Run panel regression
    temp_reg <- lm(yield ~ x + V2+V3+V4 + year + I(year^2) + as.factor(Farm)-1, data = temp1)
    list_farm_wheat_panelmodelkn5_models.of.farm[[m]] <- temp_reg
    
    rm(temp1,temp_reg)
    
  }
  
  list_farm_wheat_panelmodelkn5_modeloutputs [[i]] <- list_farm_wheat_panelmodelkn5_models.of.farm   
  rm(list_farm_wheat_panelmodelkn5_models.of.farm)
}


# Calculation of hourly effects (marginal temperature effect)
list_farm_wheat_hourly_effect_kn5 <- list() # [[farm]] [[model]][]

for (i in 1:nrow(farm_yields[[1]])){
  
  temp_list <- list()
  for (m in 1:ncol(knot_locations_kn5_wheat)){
    
    # Subset of farm i
    temp1 <- which(list_wheat_matrix_modelkn5_hourly[[m]][,"Farm"] == i) 
    temp2 <- list_wheat_matrix_modelkn5_hourly[[m]][temp1,]  
    
    # Calculate temperature effect
    temp2$effect <- as.matrix(temp2[,6:9]) %*% coef(list_farm_wheat_panelmodelkn5_modeloutputs [[i]][[m]])[1:4]
    
    temp_list[[m]] <- temp2  
    rm(temp1,temp2)  
  }
  list_farm_wheat_hourly_effect_kn5[[i]] <- temp_list
  rm(temp_list)
}

# -------------------------------------------
# Historical payouts & premium
# -------------------------------------------

# Different strike level temperatures
strike_temperature <- c(seq(13,36,by=1))

# Historical payouts with following structure: [[i]] == of farm i, [[m]] with model specification m, 
# [[m]] contains a data.frame with [strike level, year]
list_farm_wheat_payouts_kn5 <- list() 
for (i in 1:nrow(farm_yields[[1]])){
  
  # temp list for specification[[m]]
  temp_list <- list() 
  for (m in 1:ncol(knot_locations_kn5_wheat)){
    
    temp_model <- list_farm_wheat_hourly_effect_kn5 [[i]] [[m]]
    
    # DF for payouts
    temp_df <- as.data.frame(matrix(NA, nrow=length(strike_temperature), ncol=length(seq(1995,last_year,1))))
    row.names(temp_df) <- strike_temperature 
    colnames(temp_df) <- seq(1995,last_year,1)
    
    for (t in 1:ncol(temp_df)){
      # Subset of year
      sub_model_df <- subset(temp_model, temp_model$year == as.numeric(colnames(temp_df)[t]))
      # Subset only negative expected impacts trigger a payout 
      sub2_model_df <- subset(sub_model_df, sub_model_df$effect < 0)
      
      for (strike in 1:length(strike_temperature)){
        
        # Subset only temperature above the strike level trigger a payout
        sub3_model_df <- subset(sub2_model_df, sub2_model_df$Temperature >= strike_temperature[strike]) 
        
        # Calculate the cumulative payout
        temp_df[strike,t] <- sum(sub3_model_df$effect) * (-1)
        rm(sub3_model_df)
      }
      
      rm(sub2_model_df, sub_model_df)  
    }
    temp_list[[m]] <- temp_df
    rm(temp_df, temp_model)
  }
  list_farm_wheat_payouts_kn5[[i]] <- temp_list
  rm(temp_list)
} 

# Calculate the premium
# Structure: [farm i, model m, strike level]
array_premium_wheat_kn5 <- array(dim=c(nrow(farm_yields[[1]]), ncol(knot_locations_kn5_wheat),length(strike_temperature)))
dimnames(array_premium_wheat_kn5)[[3]] <- strike_temperature

for (i in 1:nrow(farm_yields[[1]])){
  for (m in 1:ncol(knot_locations_kn5_wheat)){
    for (strike in 1:length(strike_temperature)){
      temp1 <- which(!is.na(detr_wheat_yields[i,]))
      array_premium_wheat_kn5[i,m,strike] <- mean(as.numeric(list_farm_wheat_payouts_kn5[[i]][[m]][strike, temp1]), na.rm=T)
      rm(temp1)
    }  
  }
}

# -------------------------------------------
# Insured revenue and expected utility
# -------------------------------------------

# Calculation of revenue [[strike]][m,i,t]
list_farm_wheat_wealth_modelskn5_strike <- list()

for (strike in 1:length(strike_temperature)){
  
  array_farm_wheat_wealth_modelskn5 <- array(dim=c(length(knot_locations_kn5_wheat), nrow(farm_yields[[1]]), length(seq(1995,last_year,1))))
  dimnames(array_farm_wheat_wealth_modelskn5) [[3]] <-seq(1995,last_year,1) 
  
  for (m in 1:ncol(knot_locations_kn5_wheat)){
    for (i in 1:nrow(farm_yields[[1]])){
      for (t in 1: length(seq(1995,last_year,1))){
        
        array_farm_wheat_wealth_modelskn5 [m,i,t] <- as.numeric(detr_wheat_yields[i,t]) + as.numeric(list_farm_wheat_payouts_kn5 [[i]][[m]] [strike,t]) - as.numeric(array_premium_wheat_kn5[i,m,strike]) 
        
      }
    }
  }
  
  list_farm_wheat_wealth_modelskn5_strike [[strike]] <- array_farm_wheat_wealth_modelskn5 
  rm(array_farm_wheat_wealth_modelskn5)
}

# Utility of being insured with specification m.
# List has structure: [[strike]][m,i,t]
list_farm_wheat_wealth_modelkn5_strike_utility <- list()

for (strike in 1:length(strike_temperature)){
  
  temp_array <- array(dim=c(length(knot_locations_kn5_wheat), nrow(farm_yields[[1]]), length(seq(1995,last_year,1))))
  dimnames(temp_array) [[3]] <-seq(1995,last_year,1) 
  
  for (m in 1:ncol(knot_locations_kn5_wheat)){
    for (i in 1:nrow(farm_yields[[1]])){
      for (t in 1: length(seq(1995,last_year,1))){
        
        temp_array [m,i,t] <- utility_function(alpha = alpha, revenue = list_farm_wheat_wealth_modelskn5_strike [[strike]][m,i,t])
        
      }
    }
  }
  
  list_farm_wheat_wealth_modelkn5_strike_utility  [[strike]] <- temp_array 
  rm(temp_array)
}

# Expected utility statitsics of being insured
# Structure of list: 
summary_EU_wheat_insured_kn5 <- list() # [[strike]] [m,i,c(EU,E,CE,RP)]

for (strike in 1:length(strike_temperature)){
  temp_array <- array(dim=c(length(knot_locations_kn5_wheat), nrow(farm_yields[[1]]), 4))
  dimnames(temp_array) [[3]] <- c("EU", "Expected Revenue", "CE", "RP")
  
  for (m in 1:ncol(knot_locations_kn5_wheat)){
    for (i in 1:nrow(farm_yields[[1]])){
      
      # Expected utility
      temp_array  [m,i,1] <- mean(as.numeric(list_farm_wheat_wealth_modelkn5_strike_utility[[strike]][m,i,]), na.rm=T)
      
      # Expected Revenue
      temp_array  [m,i,2] <- mean(as.numeric(list_farm_wheat_wealth_modelskn5_strike[[strike]][m,i,]), na.rm=T)
      
      # Certainty Equivalent
      temp_array  [m,i,3] <- inverse_utility_function(alpha = alpha, eu = as.numeric(temp_array[m,i,1]))
      
      # Risk Premium
      temp_array  [m,i,4] <- temp_array  [m,i,2] -  temp_array  [m,i,3]
    }
  }
  summary_EU_wheat_insured_kn5 [[strike]] <- temp_array
  rm(temp_array)
}

# -------------------------------------------
# Testing
# -------------------------------------------

average_rel_change_wheat_knot5_uninsured <- as.data.frame(matrix(nrow=ncol(knot_locations_kn5_wheat), ncol=length(strike_temperature)))
colnames(average_rel_change_wheat_knot5_uninsured) <- strike_temperature

for (m in 1:ncol(knot_locations_kn5_wheat)){
  for (strike in 1:length(strike_temperature)){
    
    average_rel_change_wheat_knot5_uninsured [m,strike] <- mean((summary_EU_wheat_insured_kn5[[strike]][m,,4] - summary_EU_wheat_uninsured[,4]) / summary_EU_wheat_uninsured[,4])
    
  }
}

average_rel_change_wheat_knot5_uninsured_test <- as.data.frame(matrix(nrow=ncol(knot_locations_kn5_wheat), ncol=length(strike_temperature)))
colnames(average_rel_change_wheat_knot5_uninsured_test) <- strike_temperature

for (m in 1:ncol(knot_locations_kn5_wheat)){
  for (strike in 1:length(strike_temperature)){
    
    average_rel_change_wheat_knot5_uninsured_test [m,strike] <-wilcox.test(summary_EU_wheat_insured_kn5[[strike]][m,,4], summary_EU_wheat_uninsured[,4], paired=T, alternative="l")$p.value
    
  }
}



# ==============================================================================================================
# 10. Winter rapeseed: out-of-sample calibration and testing for five knots
# ==============================================================================================================

# -------------------------------------------
# Definition of knots
# -------------------------------------------

# Here only 1 model with 5 knots
knot_locations_kn5_rapeseed            <- as.data.frame(matrix(NA, nrow=5, ncol=1))
row.names(knot_locations_kn5_rapeseed) <- c("knot1", "knot2", "knot3","knot4","knot5")
colnames(knot_locations_kn5_rapeseed)  <- c("Schlenker") 

# Equally spaced across temperature support
knot_locations_kn5_rapeseed [,1] <- c(5,10,15,20,25)

'
# 5%, 27.5%, 50%, 72.5% and 95% quantile
knot_locations_kn5_rapeseed [,2] <- NA
  
  
  c(quantile(temp_rapeseed_exposure2$Temperature,0.05, type=1),
                                    quantile(temp_rapeseed_exposure2$Temperature,0.275, type=1),
                                    quantile(temp_rapeseed_exposure2$Temperature,0.5, type=1),
                                    quantile(temp_rapeseed_exposure2$Temperature,0.725, type=1),
                                    quantile(temp_rapeseed_exposure2$Temperature,0.95, type=1))
                                    '

# -------------------------------------------
# Out-of-sample calibration: leave-out-farm i
# -------------------------------------------

# Yields together with new hourly temperature time series: [[m]] == new time series for knot specification m
list_rapeseed_matrix_modelkn5_aggregated <- list() 

# New time series of hourly temperature exposures: [[m]] == new time series for knot specification m
list_rapeseed_matrix_modelkn5_hourly <- list() 

for (m in 1:ncol(knot_locations_kn5_rapeseed)){
  
  # Get new time series for RCS   
  temp1 <- as.matrix(rcspline.eval(temp_rapeseed_exposure2$Temperature,knots= knot_locations_kn5_rapeseed[,m], inclx=T))
  temp2 <- cbind(temp_rapeseed_exposure2, temp1)
  
  # Aggregate to yearly values
  temp3<-temp2 %>%
    group_by(Farm,year)%>%
    summarise_at(vars(x,V2,V3,V4),sum)
  
  list_rapeseed_matrix_modelkn5_hourly[[m]]     <- temp2
  list_rapeseed_matrix_modelkn5_aggregated[[m]] <- join(rapeseed_yield_melted,temp3, by=c("Farm", "year"))
  
  rm(temp1,temp2,temp3)
}

# Now, we have all data together (yields and temperature time series)
# We continue with regression. See the paper for more details.

# List with model outputs. [[i]] == farm i; [[m]] == knot specification m
list_farm_rapeseed_panelmodelkn5_modeloutputs <- list()

for (i in 1:nrow(farm_yields[[2]])){
  
  # Contains the m models of farm i.
  list_farm_rapeseed_panelmodelkn5_models.of.farm <- list() 
  for (m in 1:ncol(knot_locations_kn5_rapeseed)){
    
    # Aggregated data without farm i (leave-out-farm i for out-of-sample training)
    temp1 <- subset(list_rapeseed_matrix_modelkn5_aggregated[[m]], list_rapeseed_matrix_modelkn5_aggregated[[m]] [,"Farm"] != i)
    
    # Run panel regression
    temp_reg <- lm(yield ~ x + V2+V3+V4 + year + I(year^2) + as.factor(Farm)-1, data = temp1)
    list_farm_rapeseed_panelmodelkn5_models.of.farm[[m]] <- temp_reg
    
    rm(temp1,temp_reg)
    
  }
  
  list_farm_rapeseed_panelmodelkn5_modeloutputs [[i]] <- list_farm_rapeseed_panelmodelkn5_models.of.farm   
  rm(list_farm_rapeseed_panelmodelkn5_models.of.farm)
}


# Calculation of hourly effects (marginal temperature effect)
list_farm_rapeseed_hourly_effect_kn5 <- list() # [[farm]] [[model]][]

for (i in 1:nrow(farm_yields[[2]])){
  
  temp_list <- list()
  for (m in 1:ncol(knot_locations_kn5_rapeseed)){
    
    # Subset of farm i
    temp1 <- which(list_rapeseed_matrix_modelkn5_hourly[[m]][,"Farm"] == i) 
    temp2 <- list_rapeseed_matrix_modelkn5_hourly[[m]][temp1,]  
    
    # Calculate temperature effect
    temp2$effect <- as.matrix(temp2[,6:9]) %*% coef(list_farm_rapeseed_panelmodelkn5_modeloutputs [[i]][[m]])[1:4]
    
    temp_list[[m]] <- temp2  
    rm(temp1,temp2)  
  }
  list_farm_rapeseed_hourly_effect_kn5[[i]] <- temp_list
  rm(temp_list)
}

# -------------------------------------------
# Historical payouts & premium
# -------------------------------------------

# Different strike level temperatures
strike_temperature <- c(seq(13,36,by=1))

# Historical payouts with following structure: [[i]] == of farm i, [[m]] with model specification m, 
# [[m]] contains a data.frame with [strike level, year]
list_farm_rapeseed_payouts_kn5 <- list() 
for (i in 1:nrow(farm_yields[[2]])){
  
  # temp list for specification[[m]]
  temp_list <- list() 
  for (m in 1:ncol(knot_locations_kn5_rapeseed)){
    
    temp_model <- list_farm_rapeseed_hourly_effect_kn5 [[i]] [[m]]
    
    # DF for payouts
    temp_df <- as.data.frame(matrix(NA, nrow=length(strike_temperature), ncol=length(seq(1995,last_year,1))))
    row.names(temp_df) <- strike_temperature 
    colnames(temp_df) <- seq(1995,last_year,1)
    
    for (t in 1:ncol(temp_df)){
      # Subset of year
      sub_model_df <- subset(temp_model, temp_model$year == as.numeric(colnames(temp_df)[t]))
      # Subset only negative expected impacts trigger a payout 
      sub2_model_df <- subset(sub_model_df, sub_model_df$effect < 0)
      
      for (strike in 1:length(strike_temperature)){
        
        # Subset only temperature above the strike level trigger a payout
        sub3_model_df <- subset(sub2_model_df, sub2_model_df$Temperature >= strike_temperature[strike]) 
        
        # Calculate the cumulative payout
        temp_df[strike,t] <- sum(sub3_model_df$effect) * (-1)
        rm(sub3_model_df)
      }
      
      rm(sub2_model_df, sub_model_df)  
    }
    temp_list[[m]] <- temp_df
    rm(temp_df, temp_model)
  }
  list_farm_rapeseed_payouts_kn5[[i]] <- temp_list
  rm(temp_list)
} 

# Calculate the premium
# Structure: [farm i, model m, strike level]
array_premium_rapeseed_kn5 <- array(dim=c(nrow(farm_yields[[1]]), ncol(knot_locations_kn5_rapeseed),length(strike_temperature)))
dimnames(array_premium_rapeseed_kn5)[[3]] <- strike_temperature

for (i in 1:nrow(farm_yields[[2]])){
  for (m in 1:ncol(knot_locations_kn5_rapeseed)){
    for (strike in 1:length(strike_temperature)){
      temp1 <- which(!is.na(detr_rapeseed_yields[i,]))
      array_premium_rapeseed_kn5[i,m,strike] <- mean(as.numeric(list_farm_rapeseed_payouts_kn5[[i]][[m]][strike,temp1 ]), na.rm=T)
      rm(temp1)
    }  
  }
}

# -------------------------------------------
# Insured revenue and expected utility
# -------------------------------------------

# Calculation of revenue [[strike]][m,i,t]
list_farm_rapeseed_wealth_modelskn5_strike <- list()

for (strike in 1:length(strike_temperature)){
  
  array_farm_rapeseed_wealth_modelskn5 <- array(dim=c(length(knot_locations_kn5_rapeseed), nrow(farm_yields[[2]]), length(seq(1995,last_year,1))))
  dimnames(array_farm_rapeseed_wealth_modelskn5) [[3]] <-seq(1995,last_year,1) 
  
  for (m in 1:ncol(knot_locations_kn5_rapeseed)){
    for (i in 1:nrow(farm_yields[[2]])){
      for (t in 1: length(seq(1995,last_year,1))){
        
        array_farm_rapeseed_wealth_modelskn5 [m,i,t] <- as.numeric(detr_rapeseed_yields[i,t]) + as.numeric(list_farm_rapeseed_payouts_kn5 [[i]][[m]] [strike,t]) - as.numeric(array_premium_rapeseed_kn5[i,m,strike]) 
        
      }
    }
  }
  
  list_farm_rapeseed_wealth_modelskn5_strike [[strike]] <- array_farm_rapeseed_wealth_modelskn5 
  rm(array_farm_rapeseed_wealth_modelskn5)
}

# Utility of being insured with specification m.
# List has structure: [[strike]][mi,i,t]
list_farm_rapeseed_wealth_modelkn5_strike_utility <- list()

for (strike in 1:length(strike_temperature)){
  
  temp_array <- array(dim=c(length(knot_locations_kn5_rapeseed), nrow(farm_yields[[2]]), length(seq(1995,last_year,1))))
  dimnames(temp_array) [[3]] <-seq(1995,last_year,1) 
  
  for (m in 1:ncol(knot_locations_kn5_rapeseed)){
    for (i in 1:nrow(farm_yields[[2]])){
      for (t in 1: length(seq(1995,last_year,1))){
        
        temp_array [m,i,t] <- utility_function(alpha = alpha, revenue = list_farm_rapeseed_wealth_modelskn5_strike [[strike]][m,i,t])
        
      }
    }
  }
  
  list_farm_rapeseed_wealth_modelkn5_strike_utility  [[strike]] <- temp_array 
  rm(temp_array)
}

# Expected utility statitsics of being insured
# Structure of list: 
summary_EU_rapeseed_insured_kn5 <- list() # [[strike]] [m,i,c(EU,E,CE,RP)]

for (strike in 1:length(strike_temperature)){
  temp_array <- array(dim=c(length(knot_locations_kn5_rapeseed), nrow(farm_yields[[2]]), 4))
  dimnames(temp_array) [[3]] <- c("EU", "Expected Revenue", "CE", "RP")
  
  for (m in 1:ncol(knot_locations_kn5_rapeseed)){
    for (i in 1:nrow(farm_yields[[2]])){
      
      # Expected utility
      temp_array  [m,i,1] <- mean(as.numeric(list_farm_rapeseed_wealth_modelkn5_strike_utility[[strike]][m,i,]), na.rm=T)
      
      # Expected Revenue
      temp_array  [m,i,2] <- mean(as.numeric(list_farm_rapeseed_wealth_modelskn5_strike[[strike]][m,i,]), na.rm=T)
      
      # Certainty Equivalent
      temp_array  [m,i,3] <- inverse_utility_function(alpha = alpha, eu = as.numeric(temp_array[m,i,1]))
      
      # Risk Premium
      temp_array  [m,i,4] <- temp_array  [m,i,2] -  temp_array  [m,i,3]
    }
  }
  summary_EU_rapeseed_insured_kn5 [[strike]] <- temp_array
  rm(temp_array)
}

# -------------------------------------------
# Testing
# -------------------------------------------

average_rel_change_rapeseed_kn5_uninsured <- as.data.frame(matrix(nrow=ncol(knot_locations_kn5_rapeseed), ncol=length(strike_temperature)))
colnames(average_rel_change_rapeseed_kn5_uninsured) <- strike_temperature

for (m in 1:ncol(knot_locations_kn5_rapeseed)){
  for (strike in 1:length(strike_temperature)){
    
    average_rel_change_rapeseed_kn5_uninsured [m,strike] <- mean((summary_EU_rapeseed_insured_kn5[[strike]][m,,4] - summary_EU_rapeseed_uninsured[,4]) / summary_EU_rapeseed_uninsured[,4])
    
  }
}

average_rel_change_rapeseed_kn5_uninsured_test <- as.data.frame(matrix(nrow=ncol(knot_locations_kn5_rapeseed), ncol=length(strike_temperature)))
colnames(average_rel_change_rapeseed_kn5_uninsured_test) <- strike_temperature

for (m in 1:ncol(knot_locations_kn5_rapeseed)){
  for (strike in 1:length(strike_temperature)){
    
    average_rel_change_rapeseed_kn5_uninsured_test [m,strike] <- wilcox.test(summary_EU_rapeseed_insured_kn5[[strike]][m,,4], summary_EU_rapeseed_uninsured[,4], paired=T, alternative="l")$p.value
    
  }
}
