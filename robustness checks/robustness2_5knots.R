# =========================================================================================================
# ---------------------------------------------------------------------------------------------------------
#
# Temperature Effects on Crop Yields in Heat Index Insurance
#
# Supplementary R code for robustness check with 5 knots (out-of-sample calibration & risk evaluation)
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

# Before you run the robustness check, run main_codes_masterfile.R

# ==============================================================================================================
# Winter wheat: out-of-sample calibration and testing for five knots
# ==============================================================================================================

# -------------------------------------------
# Definition of knots (out-of-sample)
# -------------------------------------------

# List with knots for each farm and model: [[i]] == farm; [knot,model]
list_knot_locations_kn5_wheat <- list()

for (i in 1:nrow(farm_yields[[1]])){
  
  # Leave-out data from farm i
  temp_sample <- temp_wheat_exposure2[which(temp_wheat_exposure2$Farm != i),]
  temp_knot_df <- as.data.frame(matrix(NA, nrow=5, ncol=4))
  row.names(temp_knot_df) <- c("knot1", "knot2", "knot3","knot4","knot5")
  colnames(temp_knot_df)  <- c("Best fit","Equal", "Quantile","Schlenker") 
  
  # Model 1: Best fit
  
  # lower and upper bound define lowest and highest knot location
  lower_bound <- ceiling(quantile(temp_sample$Temperature,0.05, type=1))
  upper_bound <- floor(quantile(temp_sample$Temperature,0.95, type=1))
  temp_range <- seq(lower_bound, upper_bound,1) 
  kn5_combi <- combn(temp_range,5, simplify = T)
  
  # Minimum space between 2 knots
  required_space <- 4
  
  # Differences between knots
  differences <- matrix(NA, nrow=10, ncol=ncol(kn5_combi))
  differences[1,] <- abs(kn5_combi[1,] - kn5_combi[2,])
  differences[2,] <- abs(kn5_combi[1,] - kn5_combi[3,])
  differences[3,] <- abs(kn5_combi[1,] - kn5_combi[4,])
  differences[4,] <- abs(kn5_combi[1,] - kn5_combi[5,])
  
  differences[5,] <- abs(kn5_combi[2,] - kn5_combi[3,])
  differences[6,] <- abs(kn5_combi[2,] - kn5_combi[4,])
  differences[7,] <- abs(kn5_combi[2,] - kn5_combi[5,])
  
  differences[8,] <- abs(kn5_combi[3,] - kn5_combi[4,])
  differences[9,] <- abs(kn5_combi[3,] - kn5_combi[5,])
  
  differences[10,] <- abs(kn5_combi[4,] - kn5_combi[5,])
  
  knots_5_combi_good <- kn5_combi[,which(differences[1,] >= required_space & differences[2,] >= required_space & differences[3,] >= required_space & differences[4,] >= required_space & differences[5,] >= required_space & differences[6,] >= required_space & differences[7,] >= required_space & differences[8,] >= required_space & differences[9,] >= required_space & differences[10,] >= required_space  )]
  
  # vector containing residual sum of squares (RSS)
  RSS_5knots <- vector(length=ncol(knots_5_combi_good ))
  
  # Get the RSS for each model
  for (c in 1:ncol(knots_5_combi_good)){
    
    # Get new time series and aggregate hourly values  
    temp1 <- as.matrix(rcspline.eval(temp_sample$Temperature,knots= knots_5_combi_good[,c], inclx=T))
    temp2 <- cbind(temp_sample,temp1)
    
    temp3<-temp2 %>%
      group_by(Farm,year)%>%
      summarise_at(vars(x,V2,V3,V4),sum)
    
    # Match yearly values with yearly yields
    sub_wheat_yield_melted <- wheat_yield_melted[which(wheat_yield_melted$Farm != i),]
    temp4 <- join(sub_wheat_yield_melted,temp3, by=c("Farm", "year"))
    
    # Run regression
    temp_reg <- lm(yield ~ x + V2 + V3 + V4 + year + I(year^2) + as.factor(Farm)-1, data = temp4)
    RSS_5knots[c] <- sum(resid(temp_reg)^2)
    
    rm(temp_reg,temp4,sub_wheat_yield_melted,temp3,temp2,temp1)
    
    print(paste(c/ncol(knots_5_combi_good),"for farm i", i/nrow(farm_yields[[1]])))
    
  }
  
  # Best number of knots
  temp_knot_df[,1] <- knots_5_combi_good[,which.min(RSS_5knots)]
  rm(RSS_5knots, knots_5_combi_good, differences, kn5_combi,temp_range,upper_bound, lower_bound)
  
  
  # Model 2: equally spaced
  temp_knot_df[,2] <- c(min(temp_sample$Temperature) + ((max(temp_sample$Temperature)-min(temp_sample$Temperature)) / 6),
                        min(temp_sample$Temperature) + 2* ((max(temp_sample$Temperature)-min(temp_sample$Temperature)) / 6),
                        min(temp_sample$Temperature) + 3* ((max(temp_sample$Temperature)-min(temp_sample$Temperature)) / 6),
                        min(temp_sample$Temperature) + 4* ((max(temp_sample$Temperature)-min(temp_sample$Temperature)) / 6),
                        min(temp_sample$Temperature) + 5* ((max(temp_sample$Temperature)-min(temp_sample$Temperature)) / 6))
                   
  
  
  # Model 3: 10%, 50%, 90% quantile
  temp_knot_df[,3] <- c(quantile(temp_sample$Temperature,0.05, type=1),
                        quantile(temp_sample$Temperature,0.275, type=1),
                        quantile(temp_sample$Temperature,0.5, type=1),
                        quantile(temp_sample$Temperature,0.725, type=1),
                        quantile(temp_sample$Temperature,0.95, type=1))
  
  # Model 4: Schlenker
  temp_knot_df[,4] <- c(5,10,15,20,25)
  
  
  list_knot_locations_kn5_wheat[[i]] <-  temp_knot_df
  rm(temp_sample, temp_knot_df)
  print(i / nrow(farm_yields[[1]]))
}

# Testing whether model with best fit is better than linear model
# Leave-out data from farm i

# 1 = cubic spline is superior to linear model
comparison_vector <- vector()

for (i in 1:nrow(farm_yields[[1]])){
  temp_sample <- temp_wheat_exposure2[which(temp_wheat_exposure2$Farm != i),]
  
  # The RCS model with 5 knots
  # Get new time series and aggregate hourly values  
  temp1 <- as.matrix(rcspline.eval(temp_sample$Temperature,knots= list_knot_locations_kn5_wheat[[i]][,1], inclx=T))
  temp2 <- cbind(temp_sample,temp1)
  
  temp3<-temp2 %>%
    group_by(Farm,year)%>%
    summarise_at(vars(x,V2,V3,V4),sum)
  
  # Match yearly values with yearly yields
  sub_wheat_yield_melted <- wheat_yield_melted[which(wheat_yield_melted$Farm != i),]
  temp4 <- join(sub_wheat_yield_melted,temp3, by=c("Farm", "year"))
  
  # Run cubic model (restricted cubic spline model)
  temp_cub_reg <- lm(yield ~ x + V2 + V3 +V4 + year + I(year^2) + as.factor(Farm)-1, data = temp4)
  
  
  # The RCS model with 3 knots
  # Get new time series and aggregate hourly values 
  
  temp_sample_3kn <- as.matrix(rcspline.eval(temp_sample$Temperature,knots= list_knot_locations_kn3_wheat[[i]][,1], inclx=T))
  temp2_kn3 <- cbind(temp_sample,temp_sample_3kn)
  
  temp3_kn3<-temp2_kn3 %>%
    group_by(Farm,year)%>%
    summarise_at(vars(x,V2),sum)
  
  temp4_kn3 <- join(sub_wheat_yield_melted,temp3_kn3, by=c("Farm", "year"))
  
  temp_cub_reg_kn3 <- lm(yield ~ x + V2  + year + I(year^2) + as.factor(Farm)-1, data = temp4_kn3)
  
  # Compare models
  
  AIC_cub_kn3 <- AIC(temp_cub_reg_kn3)
  AIC_cub_kn5 <- AIC(temp_cub_reg)
  
  if (AIC_cub_kn5 < AIC_cub_kn3){
    comparison_vector[i] <- 1
  } else {comparison_vector[i] <- 0}
  
  rm(temp_sample, temp1, temp2, temp3,sub_wheat_yield_melted,temp_cub_reg, temp_lin_reg,AIC_cub_kn3,  AIC_cub_kn5)
  rm(temp_sample_3kn,  temp2_kn3, temp3_kn3, temp4_kn3)
  print(i / nrow(farm_yields[[1]]))
}

# Check 
comparison_vector
rm(comparison_vector)


# -------------------------------------------
# Out-of-sample calibration: leave-out-farm i
# -------------------------------------------

# Yields together with new hourly temperature time series: [[i]] == farm; [[m]] == new time series for knot specification m
list_wheat_matrix_modelkn5_aggregated <- list() # [[i]][[m]]

# New time series of hourly temperature exposures: [[i]] == farm; [[m]] == new time series for knot specification m
list_wheat_matrix_modelkn5_hourly <- list() # [[i]][[m]]

for (i in 1:nrow(farm_yields[[1]])){
  temp_hourly <- list() #[[m]]
  temp_aggregated <- list() # [[m]]
  
  for (m in 1:ncol(list_knot_locations_kn5_wheat[[i]])){
    
    # Get new time series for RCS   
    temp1 <- as.matrix(rcspline.eval(temp_wheat_exposure2$Temperature,knots= list_knot_locations_kn5_wheat[[i]][,m], inclx=T))
    temp2 <- cbind(temp_wheat_exposure2, temp1)
    
    # Aggregate to yearly values
    temp3<-temp2 %>%
      group_by(Farm,year)%>%
      summarise_at(vars(x,V2,V3,V4),sum)
    
    temp_hourly[[m]]     <- temp2
    temp_aggregated[[m]] <- join(wheat_yield_melted,temp3, by=c("Farm", "year"))
    
    rm(temp1,temp2,temp3)
  }
  list_wheat_matrix_modelkn5_aggregated[[i]] <- temp_aggregated
  list_wheat_matrix_modelkn5_hourly[[i]] <- temp_hourly
  
  rm(temp_aggregated, temp_hourly)
}

# Now, we have all data together (yields and temperature time series)
# We continue with out-of-sample regression. See the paper for more details.

# List with model outputs. [[i]] == farm i; [[m]] == knot specification m
list_farm_wheat_panelmodelkn5_modeloutputs <- list()

for (i in 1:nrow(farm_yields[[1]])){
  
  # Contains the m models of farm i.
  list_farm_wheat_panelmodelkn5_models.of.farm <- list() 
  for (m in 1:ncol(list_knot_locations_kn5_wheat[[i]])){
    
    # Aggregated data without farm i (leave-out-farm i for out-of-sample training)
    temp1 <- subset(list_wheat_matrix_modelkn5_aggregated[[i]][[m]], list_wheat_matrix_modelkn5_aggregated[[i]][[m]] [,"Farm"] != i)
    
    # Run panel regression
    temp_reg <- lm(yield ~ x + V2 +V3 + V4 + year + I(year^2) + as.factor(Farm)-1, data = temp1)
    list_farm_wheat_panelmodelkn5_models.of.farm[[m]] <- temp_reg
    
    rm(temp1,temp_reg)
    
  }
  
  list_farm_wheat_panelmodelkn5_modeloutputs [[i]] <- list_farm_wheat_panelmodelkn5_models.of.farm   
  rm(list_farm_wheat_panelmodelkn5_models.of.farm)
}


# Calculation of hourly effects (marginal temperature effect estimated at farm i)
list_farm_wheat_hourly_effect_kn5 <- list() # [[farm]] [[model]][]

for (i in 1:nrow(farm_yields[[1]])){
  
  temp_list <- list()
  for (m in 1:ncol(list_knot_locations_kn5_wheat[[i]])){
    
    # Subset of farm i
    temp1 <- which(list_wheat_matrix_modelkn5_hourly[[i]][[m]][,"Farm"] == i) 
    temp2 <- list_wheat_matrix_modelkn5_hourly[[i]][[m]][temp1,]  
    
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
  for (m in 1:ncol(list_knot_locations_kn5_wheat[[i]])){
    
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
    rm(temp_df,temp_model)
  }
  list_farm_wheat_payouts_kn5[[i]] <- temp_list
  rm(temp_list)
} 

# Calculate the premium
# Structure: [farm i, model m, strike level]
array_premium_wheat_kn5 <- array(dim=c(nrow(farm_yields[[1]]), ncol(list_knot_locations_kn5_wheat[[i]]),length(strike_temperature)))
dimnames(array_premium_wheat_kn5)[[3]] <- strike_temperature

for (i in 1:nrow(farm_yields[[1]])){
  for (m in 1:ncol(list_knot_locations_kn5_wheat[[i]])){
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
  
  array_farm_wheat_wealth_modelskn5 <- array(dim=c(ncol(list_knot_locations_kn5_wheat[[1]]), nrow(farm_yields[[1]]), length(seq(1995,last_year,1))))
  dimnames(array_farm_wheat_wealth_modelskn5) [[3]] <-seq(1995,last_year,1) 
  
  for (m in 1:ncol(list_knot_locations_kn5_wheat[[1]])){
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
# List has structure: [[strike]][mi,i,t]
list_farm_wheat_wealth_modelkn5_strike_utility <- list()

for (strike in 1:length(strike_temperature)){
  
  temp_array <- array(dim=c(ncol(list_knot_locations_kn5_wheat[[i]]), nrow(farm_yields[[1]]), length(seq(1995,last_year,1))))
  dimnames(temp_array) [[3]] <-seq(1995,last_year,1) 
  
  for (m in 1:ncol(list_knot_locations_kn5_wheat[[i]])){
    for (i in 1:nrow(farm_yields[[1]])){
      for (t in 1: length(seq(1995,last_year,1))){
        
        temp_array [m,i,t] <- utility_function(alpha = alpha, revenue = list_farm_wheat_wealth_modelskn5_strike [[strike]][m,i,t])
        
      }
    }
  }
  
  list_farm_wheat_wealth_modelkn5_strike_utility  [[strike]] <- temp_array 
  rm(temp_array)
}

# Expected utility statistics of being insured
# Structure of list: 
summary_EU_wheat_insured_kn5 <- list() # [[strike]] [m,i,c(EU,E,CE,RP)]

for (strike in 1:length(strike_temperature)){
  temp_array <- array(dim=c(ncol(list_knot_locations_kn5_wheat[[1]]), nrow(farm_yields[[1]]), 4))
  dimnames(temp_array) [[3]] <- c("EU", "Expected Revenue", "CE", "RP")
  
  for (m in 1:ncol(list_knot_locations_kn5_wheat[[1]])){
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

average_rel_change_wheat_knot5_uninsured <- as.data.frame(matrix(nrow=ncol(list_knot_locations_kn5_wheat[[1]]), ncol=length(strike_temperature)))
colnames(average_rel_change_wheat_knot5_uninsured) <- strike_temperature

for (m in 1:ncol(list_knot_locations_kn5_wheat[[1]])){
  for (strike in 1:length(strike_temperature)){
    
    average_rel_change_wheat_knot5_uninsured [m,strike] <- median((summary_EU_wheat_insured_kn5[[strike]][m,,4] - summary_EU_wheat_uninsured[,4]) / summary_EU_wheat_uninsured[,4])
    
  }
}

average_rel_change_wheat_knot5_uninsured_test <- as.data.frame(matrix(nrow=ncol(list_knot_locations_kn5_wheat[[1]]), ncol=length(strike_temperature)))
colnames(average_rel_change_wheat_knot5_uninsured_test) <- strike_temperature

for (m in 1:ncol(list_knot_locations_kn5_wheat[[1]])){
  for (strike in 1:length(strike_temperature)){
    
    average_rel_change_wheat_knot5_uninsured_test [m,strike] <- wilcox.test(summary_EU_wheat_insured_kn5[[strike]][m,,4], summary_EU_wheat_uninsured[,4], paired=T, alternative = "l")$p.value
    
  }
}

# ==============================================================================================================
# Winter rapeseed: out-of-sample calibration and testing for five knots
# ==============================================================================================================

# -------------------------------------------
# Definition of knots (out-of-sample)
# -------------------------------------------

# List with knots for each farm and model: [[i]] == farm; [knot,model]
list_knot_locations_kn5_rapeseed <- list()

for (i in 1:nrow(farm_yields[[2]])){
  
  # Leave-out data from farm i
  temp_sample <- temp_rapeseed_exposure2[which(temp_rapeseed_exposure2$Farm != i),]
  temp_knot_df <- as.data.frame(matrix(NA, nrow=5, ncol=4))
  row.names(temp_knot_df) <- c("knot1", "knot2", "knot3","knot4","knot5")
  colnames(temp_knot_df)  <- c("Best fit","Equal", "Quantile","Schlenker") 
  
  # Model 1: Best fit
  
  # lower and upper bound define lowest and highest knot location
  lower_bound <- ceiling(quantile(temp_sample$Temperature,0.05, type=1))
  upper_bound <- floor(quantile(temp_sample$Temperature,0.95, type=1))
  temp_range <- seq(lower_bound, upper_bound,1) 
  kn5_combi <- combn(temp_range,5, simplify = T)
  
  # Minimum space between 2 knots
  required_space <- 4
  
  # Differences between knots
  differences <- matrix(NA, nrow=10, ncol=ncol(kn5_combi))
  differences[1,] <- abs(kn5_combi[1,] - kn5_combi[2,])
  differences[2,] <- abs(kn5_combi[1,] - kn5_combi[3,])
  differences[3,] <- abs(kn5_combi[1,] - kn5_combi[4,])
  differences[4,] <- abs(kn5_combi[1,] - kn5_combi[5,])
  
  differences[5,] <- abs(kn5_combi[2,] - kn5_combi[3,])
  differences[6,] <- abs(kn5_combi[2,] - kn5_combi[4,])
  differences[7,] <- abs(kn5_combi[2,] - kn5_combi[5,])
  
  differences[8,] <- abs(kn5_combi[3,] - kn5_combi[4,])
  differences[9,] <- abs(kn5_combi[3,] - kn5_combi[5,])
  
  differences[10,] <- abs(kn5_combi[4,] - kn5_combi[5,])
  
  knots_5_combi_good <- kn5_combi[,which(differences[1,] >= required_space & differences[2,] >= required_space & differences[3,] >= required_space & differences[4,] >= required_space & differences[5,] >= required_space & differences[6,] >= required_space & differences[7,] >= required_space & differences[8,] >= required_space & differences[9,] >= required_space & differences[10,] >= required_space  )]
  
  # vector containing residual sum of squares (RSS)
  RSS_5knots <- vector(length=ncol(knots_5_combi_good ))
  
  # Get the RSS for each model
  for (c in 1:ncol(knots_5_combi_good)){
    
    # Get new time series and aggregate hourly values  
    temp1 <- as.matrix(rcspline.eval(temp_sample$Temperature,knots= knots_5_combi_good[,c], inclx=T))
    temp2 <- cbind(temp_sample,temp1)
    
    temp3<-temp2 %>%
      group_by(Farm,year)%>%
      summarise_at(vars(x,V2,V3,V4),sum)
    
    # Match yearly values with yearly yields
    sub_rapeseed_yield_melted <- rapeseed_yield_melted[which(rapeseed_yield_melted$Farm != i),]
    temp4 <- join(sub_rapeseed_yield_melted,temp3, by=c("Farm", "year"))
    
    # Run regression
    temp_reg <- lm(yield ~ x + V2 + V3 + V4 + year + I(year^2) + as.factor(Farm)-1, data = temp4)
    RSS_5knots[c] <- sum(resid(temp_reg)^2)
    
    rm(temp_reg,temp4,sub_rapeseed_yield_melted,temp3,temp2,temp1)
    
    print(paste(c/ncol(knots_5_combi_good),"for farm i", i/nrow(farm_yields[[2]])))
    
  }
  
  # Best number of knots
  temp_knot_df[,1] <- knots_5_combi_good[,which.min(RSS_5knots)]
  rm(RSS_5knots, knots_5_combi_good, differences, kn5_combi,temp_range,upper_bound, lower_bound)
  
  
  # Model 2: equally spaced
  temp_knot_df[,2] <- c(min(temp_sample$Temperature) + ((max(temp_sample$Temperature)-min(temp_sample$Temperature)) / 6),
                        min(temp_sample$Temperature) + 2* ((max(temp_sample$Temperature)-min(temp_sample$Temperature)) / 6),
                        min(temp_sample$Temperature) + 3* ((max(temp_sample$Temperature)-min(temp_sample$Temperature)) / 6),
                        min(temp_sample$Temperature) + 4* ((max(temp_sample$Temperature)-min(temp_sample$Temperature)) / 6),
                        min(temp_sample$Temperature) + 5* ((max(temp_sample$Temperature)-min(temp_sample$Temperature)) / 6))
  
  
  
  # Model 3: 10%, 50%, 90% quantile
  temp_knot_df[,3] <- c(quantile(temp_sample$Temperature,0.05, type=1),
                        quantile(temp_sample$Temperature,0.275, type=1),
                        quantile(temp_sample$Temperature,0.5, type=1),
                        quantile(temp_sample$Temperature,0.725, type=1),
                        quantile(temp_sample$Temperature,0.95, type=1))
  
  # Model 4: Schlenker
  temp_knot_df[,4] <- c(5,10,15,20,25)
  
  
  list_knot_locations_kn5_rapeseed[[i]] <-  temp_knot_df
  rm(temp_sample, temp_knot_df)
  print(i / nrow(farm_yields[[2]]))
}

# Testing whether model with best fit is better than linear model
# Leave-out data from farm i

# 1 = cubic spline is superior to linear model
comparison_vector <- vector()

for (i in 1:nrow(farm_yields[[2]])){
  temp_sample <- temp_rapeseed_exposure2[which(temp_rapeseed_exposure2$Farm != i),]
  
  # The RCS model with 5 knots
  # Get new time series and aggregate hourly values  
  temp1 <- as.matrix(rcspline.eval(temp_sample$Temperature,knots= list_knot_locations_kn5_rapeseed[[i]][,1], inclx=T))
  temp2 <- cbind(temp_sample,temp1)
  
  temp3<-temp2 %>%
    group_by(Farm,year)%>%
    summarise_at(vars(x,V2,V3,V4),sum)
  
  # Match yearly values with yearly yields
  sub_rapeseed_yield_melted <- rapeseed_yield_melted[which(rapeseed_yield_melted$Farm != i),]
  temp4 <- join(sub_rapeseed_yield_melted,temp3, by=c("Farm", "year"))
  
  # Run cubic model (restricted cubic spline model)
  temp_cub_reg <- lm(yield ~ x + V2 + V3 +V4 + year + I(year^2) + as.factor(Farm)-1, data = temp4)
  
  
  # The RCS model with 3 knots
  # Get new time series and aggregate hourly values 
  
  temp_sample_3kn <- as.matrix(rcspline.eval(temp_sample$Temperature,knots= list_knot_locations_kn3_rapeseed[[i]][,1], inclx=T))
  temp2_kn3 <- cbind(temp_sample,temp_sample_3kn)
  
  temp3_kn3<-temp2_kn3 %>%
    group_by(Farm,year)%>%
    summarise_at(vars(x,V2),sum)
  
  temp4_kn3 <- join(sub_rapeseed_yield_melted,temp3_kn3, by=c("Farm", "year"))
  
  temp_cub_reg_kn3 <- lm(yield ~ x + V2  + year + I(year^2) + as.factor(Farm)-1, data = temp4_kn3)
  
  # Compare models
  
  AIC_cub_kn3 <- AIC(temp_cub_reg_kn3)
  AIC_cub_kn5 <- AIC(temp_cub_reg)
  
  if (AIC_cub_kn5 < AIC_cub_kn3){
    comparison_vector[i] <- 1
  } else {comparison_vector[i] <- 0}
  
  rm(temp_sample, temp1, temp2, temp3,sub_rapeseed_yield_melted,temp_cub_reg, temp_lin_reg,AIC_cub_kn3,  AIC_cub_kn5)
  rm(temp_sample_3kn,  temp2_kn3, temp3_kn3, temp4_kn3)
  print(i / nrow(farm_yields[[2]]))
}

# Check 
comparison_vector
rm(comparison_vector)


# -------------------------------------------
# Out-of-sample calibration: leave-out-farm i
# -------------------------------------------

# Yields together with new hourly temperature time series: [[i]] == farm; [[m]] == new time series for knot specification m
list_rapeseed_matrix_modelkn5_aggregated <- list() # [[i]][[m]]

# New time series of hourly temperature exposures: [[i]] == farm; [[m]] == new time series for knot specification m
list_rapeseed_matrix_modelkn5_hourly <- list() # [[i]][[m]]

for (i in 1:nrow(farm_yields[[2]])){
  temp_hourly <- list() #[[m]]
  temp_aggregated <- list() # [[m]]
  
  for (m in 1:ncol(list_knot_locations_kn5_rapeseed[[i]])){
    
    # Get new time series for RCS   
    temp1 <- as.matrix(rcspline.eval(temp_rapeseed_exposure2$Temperature,knots= list_knot_locations_kn5_rapeseed[[i]][,m], inclx=T))
    temp2 <- cbind(temp_rapeseed_exposure2, temp1)
    
    # Aggregate to yearly values
    temp3<-temp2 %>%
      group_by(Farm,year)%>%
      summarise_at(vars(x,V2,V3,V4),sum)
    
    temp_hourly[[m]]     <- temp2
    temp_aggregated[[m]] <- join(rapeseed_yield_melted,temp3, by=c("Farm", "year"))
    
    rm(temp1,temp2,temp3)
  }
  list_rapeseed_matrix_modelkn5_aggregated[[i]] <- temp_aggregated
  list_rapeseed_matrix_modelkn5_hourly[[i]] <- temp_hourly
  
  rm(temp_aggregated, temp_hourly)
}

# Now, we have all data together (yields and temperature time series)
# We continue with out-of-sample regression. See the paper for more details.

# List with model outputs. [[i]] == farm i; [[m]] == knot specification m
list_farm_rapeseed_panelmodelkn5_modeloutputs <- list()

for (i in 1:nrow(farm_yields[[2]])){
  
  # Contains the m models of farm i.
  list_farm_rapeseed_panelmodelkn5_models.of.farm <- list() 
  for (m in 1:ncol(list_knot_locations_kn5_rapeseed[[i]])){
    
    # Aggregated data without farm i (leave-out-farm i for out-of-sample training)
    temp1 <- subset(list_rapeseed_matrix_modelkn5_aggregated[[i]][[m]], list_rapeseed_matrix_modelkn5_aggregated[[i]][[m]] [,"Farm"] != i)
    
    # Run panel regression
    temp_reg <- lm(yield ~ x + V2 +V3 + V4 + year + I(year^2) + as.factor(Farm)-1, data = temp1)
    list_farm_rapeseed_panelmodelkn5_models.of.farm[[m]] <- temp_reg
    
    rm(temp1,temp_reg)
    
  }
  
  list_farm_rapeseed_panelmodelkn5_modeloutputs [[i]] <- list_farm_rapeseed_panelmodelkn5_models.of.farm   
  rm(list_farm_rapeseed_panelmodelkn5_models.of.farm)
}


# Calculation of hourly effects (marginal temperature effect estimated at farm i)
list_farm_rapeseed_hourly_effect_kn5 <- list() # [[farm]] [[model]][]

for (i in 1:nrow(farm_yields[[2]])){
  
  temp_list <- list()
  for (m in 1:ncol(list_knot_locations_kn5_rapeseed[[i]])){
    
    # Subset of farm i
    temp1 <- which(list_rapeseed_matrix_modelkn5_hourly[[i]][[m]][,"Farm"] == i) 
    temp2 <- list_rapeseed_matrix_modelkn5_hourly[[i]][[m]][temp1,]  
    
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
  for (m in 1:ncol(list_knot_locations_kn5_rapeseed[[i]])){
    
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
    rm(temp_df,temp_model)
  }
  list_farm_rapeseed_payouts_kn5[[i]] <- temp_list
  rm(temp_list)
} 

# Calculate the premium
# Structure: [farm i, model m, strike level]
array_premium_rapeseed_kn5 <- array(dim=c(nrow(farm_yields[[2]]), ncol(list_knot_locations_kn5_rapeseed[[i]]),length(strike_temperature)))
dimnames(array_premium_rapeseed_kn5)[[3]] <- strike_temperature

for (i in 1:nrow(farm_yields[[2]])){
  for (m in 1:ncol(list_knot_locations_kn5_rapeseed[[i]])){
    for (strike in 1:length(strike_temperature)){
      temp1 <- which(!is.na(detr_rapeseed_yields[i,]))
      array_premium_rapeseed_kn5[i,m,strike] <- mean(as.numeric(list_farm_rapeseed_payouts_kn5[[i]][[m]][strike, temp1]), na.rm=T)
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
  
  array_farm_rapeseed_wealth_modelskn5 <- array(dim=c(ncol(list_knot_locations_kn5_rapeseed[[1]]), nrow(farm_yields[[2]]), length(seq(1995,last_year,1))))
  dimnames(array_farm_rapeseed_wealth_modelskn5) [[3]] <-seq(1995,last_year,1) 
  
  for (m in 1:ncol(list_knot_locations_kn5_rapeseed[[1]])){
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
  
  temp_array <- array(dim=c(ncol(list_knot_locations_kn5_rapeseed[[i]]), nrow(farm_yields[[2]]), length(seq(1995,last_year,1))))
  dimnames(temp_array) [[3]] <-seq(1995,last_year,1) 
  
  for (m in 1:ncol(list_knot_locations_kn5_rapeseed[[i]])){
    for (i in 1:nrow(farm_yields[[2]])){
      for (t in 1: length(seq(1995,last_year,1))){
        
        temp_array [m,i,t] <- utility_function(alpha = alpha, revenue = list_farm_rapeseed_wealth_modelskn5_strike [[strike]][m,i,t])
        
      }
    }
  }
  
  list_farm_rapeseed_wealth_modelkn5_strike_utility  [[strike]] <- temp_array 
  rm(temp_array)
}

# Expected utility statistics of being insured
# Structure of list: 
summary_EU_rapeseed_insured_kn5 <- list() # [[strike]] [m,i,c(EU,E,CE,RP)]

for (strike in 1:length(strike_temperature)){
  temp_array <- array(dim=c(ncol(list_knot_locations_kn5_rapeseed[[1]]), nrow(farm_yields[[2]]), 4))
  dimnames(temp_array) [[3]] <- c("EU", "Expected Revenue", "CE", "RP")
  
  for (m in 1:ncol(list_knot_locations_kn5_rapeseed[[1]])){
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

average_rel_change_rapeseed_knot5_uninsured <- as.data.frame(matrix(nrow=ncol(list_knot_locations_kn5_rapeseed[[1]]), ncol=length(strike_temperature)))
colnames(average_rel_change_rapeseed_knot5_uninsured) <- strike_temperature

for (m in 1:ncol(list_knot_locations_kn5_rapeseed[[1]])){
  for (strike in 1:length(strike_temperature)){
    
    average_rel_change_rapeseed_knot5_uninsured [m,strike] <- median((summary_EU_rapeseed_insured_kn5[[strike]][m,,4] - summary_EU_rapeseed_uninsured[,4]) / summary_EU_rapeseed_uninsured[,4])
    
  }
}

average_rel_change_rapeseed_knot5_uninsured_test <- as.data.frame(matrix(nrow=ncol(list_knot_locations_kn5_rapeseed[[1]]), ncol=length(strike_temperature)))
colnames(average_rel_change_rapeseed_knot5_uninsured_test) <- strike_temperature

for (m in 1:ncol(list_knot_locations_kn5_rapeseed[[1]])){
  for (strike in 1:length(strike_temperature)){
    
    average_rel_change_rapeseed_knot5_uninsured_test [m,strike] <- wilcox.test(summary_EU_rapeseed_insured_kn5[[strike]][m,,4], summary_EU_rapeseed_uninsured[,4], paired=T, alternative = "l")$p.value
    
  }
}

