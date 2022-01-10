# =========================================================================================================
# ---------------------------------------------------------------------------------------------------------
#
# Temperature Effects on Crop Yields in Heat Index Insurance
#
# Supplementary R code for plots in the main paper
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

# This file shows the codes to create figures 2-4 of the main paper-
# Run the scrip main_codes_masterfile.R before you run this script.

# Content:
# 1) Figure 2: hourly temperature effects and payouts for winter wheat
# 2) Figure 3: hourly temperature effects and payouts for winter rapeseed
# 3) Figure 4: Out-of-sample risk reductions

# Load the packages
library(fixest)
library(ggplot2)
library(gridExtra)
library(robustbase)
library(egg)
library(matrixStats)
library(tidyverse)

# ====================================================================
#
# 1) Figure 2: hourly temperature effects and payouts for winter wheat
#
# ====================================================================

# The plots are representative for a farm drawn at random.
# Note that plots slightly vary between farms because of the out-of-sample calibration
set.seed(123)
farm_ID <- round(runif(1,1,nrow(farm_yields[[1]])))

# -------------------------------------------
# Histogram of temperature exposure
# -------------------------------------------

sub_temperature_exposure_wheat <- temp_wheat_exposure2[which(temp_wheat_exposure2$Farm != farm_ID),]

temperature_wheat <- ggplot(sub_temperature_exposure_wheat, aes(x=Temperature)) + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "grey25", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "grey25"),
        axis.text = element_text(family="Times New Roman", size=10),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  ylab("Total exposure \n (1'000 hours)")+ xlab("Temperature (°C)")+
  geom_histogram(binwidth=1, color="black", fill="brown1")+
  scale_x_continuous(breaks = seq(-5,39, by=5), limits = c(-5,40))+
  scale_y_continuous(breaks = c(0,100000,200000), labels = c(0,100,200))

# -------------------------------------------
# Model 1: 
# -------------------------------------------

# Temperature time series
temperatures_ts_M1 <- rcspline.eval(sub_temperature_exposure_wheat$Temperature,knots= list_knot_locations_kn3_wheat[[farm_ID]][,1], inclx=T)  

# Add new ts to dates
temp1 <- cbind(sub_temperature_exposure_wheat, temperatures_ts_M1)

# Aggregation to yearly values
final_wheat_temperature_aggregated_M1 <-temp1 %>%
  group_by(Farm,year)%>%
  summarise_at(vars(x,V2),sum)
rm(temp1)

# Combine yearly values with yields
temp1 <- wheat_yield_melted[which(wheat_yield_melted$Farm != farm_ID),]
data_wheat_final_M1 <- join(temp1, final_wheat_temperature_aggregated_M1, by=c("Farm","year"))
rm(temperatures_ts_M1,temp1,final_wheat_temperature_aggregated_M1)

# New TS for temperature support of interest
marginal_temperatures_M1 <- rcspline.eval(seq(-10,39,0.225),knots=list_knot_locations_kn3_wheat[[farm_ID]][,1], inclx=T)

# Regression output with reliable R-squared
temp_reg1 <- feols(yield ~ x + V2 + year + I(year^2) | Farm, data=data_wheat_final_M1)
summary(temp_reg1)
rm(temp_reg, temp_reg1)

# 
# Statistical Uncertainty
# 

# Run 1000 models to derive uncertainty in estimates
spline_fitted_M1 <- matrix(NA,nrow=1000, ncol=nrow(marginal_temperatures_M1))

for (l in 1: 1000){
  # Sample some years
  data_temp <- data_wheat_final_M1[which(data_wheat_final_M1$year %in% sample(unique(data_wheat_final_M1$year),replace=T)),]
  
  # Run regression with random subsample (data_temp)
  fit_boot  <- lm( yield ~ x + V2 +year + I(year^2) + as.factor(Farm)-1, data=data_temp)
  
  # Save marginal effect of temperatures 
  spline_fitted_M1[l,]<-as.vector((marginal_temperatures_M1 %*% coef(fit_boot)[1:2]))
  
  # Show progress in %
  print(l/1000 *100)
  rm(fit_boot,data_temp)
}

# Change structure
spline_sub_melt_M1 <- melt(spline_fitted_M1)

# Get median response for each column (= a temperature from seq(3,39,0.225))
colmedians_spline_M1 <- colMedians(spline_fitted_M1)

# Get distance from median response for each value
dist_colmedian_M1 <- matrix(NA, nrow=nrow(spline_fitted_M1),ncol=ncol(spline_fitted_M1))

for(l in 1:ncol(spline_fitted_M1)){
  dist_colmedian_M1[,l] <- spline_fitted_M1[,l] - colmedians_spline_M1[l]
}

# Normalize distance to median
dist_colmedian_M1_normalized <- matrix(NA, nrow=nrow(spline_fitted_M1),ncol=ncol(spline_fitted_M1))

for(l in 1:ncol(spline_fitted_M1)){
  dist_colmedian_M1_normalized[,l]<-abs(dist_colmedian_M1[,l]/max(abs(dist_colmedian_M1[,l])))
}

# Change structure
dist_colmedian_M1_normalized_melted <- melt(dist_colmedian_M1_normalized)

# Last preparation before plotting

spline_merged_M1 <- merge(spline_sub_melt_M1,dist_colmedian_M1_normalized_melted,by=c("Var1","Var2"))
spline_merged_M1$value.y_minus1 <- 1-spline_merged_M1$value.y

temp1 <- aggregate(spline_merged_M1$value.y_minus1, list(spline_merged_M1$Var1), mean, na.rm=T) 
temp2 <- aggregate(spline_merged_M1$value.x, list(spline_merged_M1$Var2), mean, na.rm=T)

temp3 <- aggregate(spline_merged_M1$value.x, list(spline_merged_M1$Var2), quantile, probs=0.025) 
temp4 <- aggregate(spline_merged_M1$value.x, list(spline_merged_M1$Var2), quantile, probs=0.975)

spline_merged_M1 <- merge(spline_merged_M1, temp1, by.x=c("Var1"), by.y=c("Group.1"))
spline_merged_M1 <- merge(spline_merged_M1, temp2, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M1)[7]<-"average_effect"

spline_merged_M1<-merge(spline_merged_M1, temp3, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M1)[8]<-"lower_bound"

spline_merged_M1<-merge(spline_merged_M1, temp4, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M1)[9]<-"upper_bound"

temp5 <- cbind(seq(-10,39,0.225),seq(1:length(seq(-10,39,0.225))))
colnames(temp5) <- c("Temp_bin", "Var2")
spline_merged_M1 <-merge(spline_merged_M1, temp5, by=c("Var2"))

# 
# Plot: hourly temperature effect
# 

model_1 <- ggplot()+ ggtitle("Hourly temperature effects: Model 1") + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "grey25", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "grey25"),
        axis.text = element_text(family="Times New Roman", size=10),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  xlab("Temperature (°C)") + ylab("Hourly yield response \n (dt/ha)")+
  scale_x_continuous(breaks = seq(-5,39, by=5), limits = c(-5,40))+
  scale_y_continuous(breaks = c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1), limits = c(-0.4,0.1))+
  geom_hline(yintercept=c(seq(-0.5,0.1,by=0.1)), colour="gray")+
  geom_vline(xintercept=list_knot_locations_kn3_wheat[[farm_ID]][,1], linetype="dotted", colour="gray")+
  geom_line(data=spline_merged_M1, aes(y=value.x, x=Temp_bin, group=Var1, alpha=x.x^200/200), size=.01, colour="brown1")+
  geom_line(data=spline_merged_M1, aes(y=average_effect, x=Temp_bin), alpha=0.85, size=.75, colour="brown4")+
  geom_line(data=spline_merged_M1, aes(y=lower_bound, x=Temp_bin)   , alpha=0.85, size=.75, colour="brown4")+       
  geom_line(data=spline_merged_M1, aes(y=upper_bound, x=Temp_bin)   , alpha=0.85, size=.75, colour="brown4")

rm(temp1,temp2,temp3,temp4,temp5)
rm(dist_colmedian_M1_normalized_melted, colmedians_spline_M1, spline_fitted_M1,dist_colmedian_M1_normalized,dist_colmedian_M1,spline_sub_melt_M1)

# 
# Preparation of hourly payouts
# 

# We change the sign of the effect (payouts are a positive number)

spline_merged_M1$average_effect <- spline_merged_M1$average_effect * (-1) 
spline_merged_M1$lower_bound <- spline_merged_M1$lower_bound * (-1) 
spline_merged_M1$upper_bound <- spline_merged_M1$upper_bound * (-1) 
spline_merged_M1$value.x <- spline_merged_M1$value.x * (-1) 

# Prepare lower bound (no values below 0)
temp <- subset(spline_merged_M1[,c(2,3,6,7,8,9,10)])
colnames(temp)[5] <- "upper_bound"
colnames(temp)[6] <- "lower_bound"

# Remove negative values(required for geom_ribbon())

temp[which(temp$lower_bound<0),6] <- 0

payout_wheat_model_1 <- ggplot(temp) +
                        ggtitle("Hourly payouts: Model 1")+ theme(panel.grid.major = element_blank(),
                        plot.title = element_text(family="Times New Roman", size=12, colour = "grey25", hjust=0.5),
                        axis.title = element_text(family="Times New Roman", size=10, colour = "grey25"),
                        axis.text = element_text(family="Times New Roman", size=10),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        legend.position="none")+
                        xlab("Temperature (°C)") + ylab("Payout (dt/ha)")+
                        scale_x_continuous(breaks = seq(10,39, by=5), limits = c(10,40))+
                        scale_y_continuous(breaks = c(0,0.1, 0.2,0.3,0.4), limits = c(0,0.46))+
                        geom_hline(yintercept=c(seq(0,0.4,by=0.1)), colour="gray")+
                        geom_hline(yintercept=0, colour="darkgray")+
                        geom_vline(xintercept= seq(5,38, by=5), colour="gray", linetype="dashed")+
  
  # Add lines
  
  geom_line(data=temp, aes(y=value.x, x=Temp_bin, group=Var1, alpha=x.x^200/200), size=.01, colour="royalblue1")+
  geom_line(data=temp, aes(y=average_effect, x=Temp_bin), size=0.75, color="royalblue4")+
  geom_line(data=temp, aes(y=lower_bound, x=Temp_bin), size=0.75, color="royalblue4")+
  geom_line(data=temp, aes(y=upper_bound, x=Temp_bin), size=0.75, color="royalblue4")

rm(temp,spline_merged_M1)

# -------------------------------------------
# Model 2: 
# -------------------------------------------

# Temperature time series
temperatures_ts_M2 <- rcspline.eval(sub_temperature_exposure_wheat$Temperature,knots= list_knot_locations_kn3_wheat[[farm_ID]][,2], inclx=T)  

# Add new ts to dates
temp1 <- cbind(sub_temperature_exposure_wheat, temperatures_ts_M2)

# Aggregation to yearly values
final_wheat_temperature_aggregated_M2 <-temp1 %>%
  group_by(Farm,year)%>%
  summarise_at(vars(x,V2),sum)
rm(temp1)

# Combine yearly values with yields
temp1 <- wheat_yield_melted[which(wheat_yield_melted$Farm != farm_ID),]
data_wheat_final_M2 <- join(temp1, final_wheat_temperature_aggregated_M2, by=c("Farm","year"))
rm(temperatures_ts_M2,temp1,final_wheat_temperature_aggregated_M2)

# New TS for temperature support of interest
marginal_temperatures_M2 <- rcspline.eval(seq(-10,39,0.225),knots=list_knot_locations_kn3_wheat[[farm_ID]][,2], inclx=T)

# Regression output with reliable R-squared
temp_reg1 <- feols(yield ~ x + V2 + year + I(year^2) | Farm, data=data_wheat_final_M2)
summary(temp_reg1)
rm(temp_reg1)

# 
# Statistical uncertainty 
# 

# Run 1000 models to derive uncertainty in estimates
spline_fitted_M2 <- matrix(NA,nrow=1000, ncol=nrow(marginal_temperatures_M2))

for (l in 1: 1000){
  # Sample some years
  data_temp <- data_wheat_final_M2[which(data_wheat_final_M2$year %in% sample(unique(data_wheat_final_M2$year),replace=T)),]
  
  # Run regression with random subsample (data_temp)
  fit_boot  <- lm( yield ~ x + V2 + year + I(year^2)  + as.factor(Farm)-1, data=data_temp)
  
  # Save marginal effect of temperatures 
  spline_fitted_M2[l,]<-as.vector((marginal_temperatures_M2 %*% coef(fit_boot)[1:2]))
  
  # Show progress in %
  print(l/1000 *100)
  rm(fit_boot,data_temp)
}

# Change structure
spline_sub_melt_M2 <- melt(spline_fitted_M2)

# Get median response for each column (= a temperature from seq(3,39,0.225))
colmedians_spline_M2 <- colMedians(spline_fitted_M2)

# Get distance from median response for each value
dist_colmedian_M2 <- matrix(NA, nrow=nrow(spline_fitted_M2),ncol=ncol(spline_fitted_M2))

for(l in 1:ncol(spline_fitted_M2)){
  dist_colmedian_M2[,l] <- spline_fitted_M2[,l] - colmedians_spline_M2[l]
}

# Normalize distance to median
dist_colmedian_M2_normalized <- matrix(NA, nrow=nrow(spline_fitted_M2),ncol=ncol(spline_fitted_M2))

for(l in 1:ncol(spline_fitted_M2)){
  dist_colmedian_M2_normalized[,l]<-abs(dist_colmedian_M2[,l]/max(abs(dist_colmedian_M2[,l])))
}

# Change structure
dist_colmedian_M2_normalized_melted <- melt(dist_colmedian_M2_normalized)

# Last preparation before plotting

spline_merged_M2 <- merge(spline_sub_melt_M2,dist_colmedian_M2_normalized_melted,by=c("Var1","Var2"))
spline_merged_M2$value.y_minus1 <- 1-spline_merged_M2$value.y

temp1 <- aggregate(spline_merged_M2$value.y_minus1, list(spline_merged_M2$Var1), mean, na.rm=T) 
temp2 <- aggregate(spline_merged_M2$value.x, list(spline_merged_M2$Var2), mean, na.rm=T)

temp3 <- aggregate(spline_merged_M2$value.x, list(spline_merged_M2$Var2), quantile, probs=0.025) 
temp4 <- aggregate(spline_merged_M2$value.x, list(spline_merged_M2$Var2), quantile, probs=0.975)

spline_merged_M2 <- merge(spline_merged_M2, temp1, by.x=c("Var1"), by.y=c("Group.1"))
spline_merged_M2 <- merge(spline_merged_M2, temp2, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M2)[7]<-"average_effect"

spline_merged_M2<-merge(spline_merged_M2, temp3, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M2)[8]<-"lower_bound"

spline_merged_M2<-merge(spline_merged_M2, temp4, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M2)[9]<-"upper_bound"

temp5 <- cbind(seq(-10,39,0.225),seq(1:length(seq(-10,39,0.225))))
colnames(temp5) <- c("Temp_bin", "Var2")
spline_merged_M2 <-merge(spline_merged_M2, temp5, by=c("Var2"))

# 
# Plot: hourly temperature effect
# 

model_2 <- ggplot()+ ggtitle("Hourly temperature effects: Model 2") + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "grey25", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "grey25"),
        axis.text = element_text(family="Times New Roman", size=10),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  xlab("Temperature (°C)") + ylab("Hourly yield response \n (dt/ha)")+
  scale_x_continuous(breaks = seq(-5,39, by=5), limits = c(-5,40))+
  scale_y_continuous(breaks = c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1), limits = c(-0.4,0.1))+
  geom_hline(yintercept=c(seq(-0.5,0.1,by=0.1)), colour="gray")+
  geom_vline(xintercept=list_knot_locations_kn3_wheat[[farm_ID]][,2], linetype="dotted", colour="gray")+
  geom_line(data=spline_merged_M2, aes(y=value.x, x=Temp_bin, group=Var1, alpha=x.x^200/200), size=.01, colour="brown1")+
  geom_line(data=spline_merged_M2, aes(y=average_effect, x=Temp_bin), alpha=0.85, size=.75, colour="brown4")+
  geom_line(data=spline_merged_M2, aes(y=lower_bound, x=Temp_bin)   , alpha=0.85, size=.75, colour="brown4")+       
  geom_line(data=spline_merged_M2, aes(y=upper_bound, x=Temp_bin)   , alpha=0.85, size=.75, colour="brown4")  

rm(temp1,temp2,temp3,temp4,temp5)
rm(dist_colmedian_M2_normalized_melted, colmedians_spline_M2, spline_fitted_M2,dist_colmedian_M2_normalized,dist_colmedian_M2,spline_sub_melt_M2)

# 
# Preparation of hourly payouts
# 

# We change the sign of the effect (payouts are a positive number)

spline_merged_M2$average_effect <- spline_merged_M2$average_effect * (-1) 
spline_merged_M2$lower_bound <- spline_merged_M2$lower_bound * (-1) 
spline_merged_M2$upper_bound <- spline_merged_M2$upper_bound * (-1) 
spline_merged_M2$value.x <- spline_merged_M2$value.x * (-1) 

# Prepare lower bound (no values below 0)
temp <- subset(spline_merged_M2[,c(2,3,6,7,8,9,10)])
colnames(temp)[5] <- "upper_bound"
colnames(temp)[6] <- "lower_bound"

# Remove negative values(required for geom_ribbon())

temp[which(temp$lower_bound<0),6] <- 0

# 
# Plot: hourly payout function
# 

payout_wheat_model_2 <- ggplot(temp) +
  ggtitle("Hourly payouts: Model 2")+ theme(panel.grid.major = element_blank(),
                                                    plot.title = element_text(family="Times New Roman", size=12, colour = "grey25", hjust=0.5),
                                                    axis.title = element_text(family="Times New Roman", size=10, colour = "grey25"),
                                                    axis.text = element_text(family="Times New Roman", size=10),
                                                    panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(),
                                                    axis.line = element_line(colour = "black"),
                                                    legend.position="none")+
  xlab("Temperature (°C)") + ylab("Payout (dt/ha)")+
  scale_x_continuous(breaks = seq(10,39, by=5), limits = c(10,40))+
  scale_y_continuous(breaks = c(0,0.1, 0.2,0.3,0.4), limits = c(0,0.46))+
  geom_hline(yintercept=c(seq(0,0.4,by=0.1)), colour="gray")+
  geom_hline(yintercept=0, colour="darkgray")+
  geom_vline(xintercept= seq(5,38, by=5), colour="gray", linetype="dashed")+
  
  # Add lines
  
  geom_line(data=temp, aes(y=value.x, x=Temp_bin, group=Var1, alpha=x.x^200/200), size=.01, colour="royalblue1")+
  geom_line(data=temp, aes(y=average_effect, x=Temp_bin), size=0.75, color="royalblue4")+
  geom_line(data=temp, aes(y=lower_bound, x=Temp_bin), size=0.75, color="royalblue4")+
  geom_line(data=temp, aes(y=upper_bound, x=Temp_bin), size=0.75, color="royalblue4")

rm(temp,spline_merged_M2)

# -------------------------------------------
# Model 3: 
# -------------------------------------------

# Temperature time series
temperatures_ts_M3 <- rcspline.eval(sub_temperature_exposure_wheat$Temperature,knots= list_knot_locations_kn3_wheat[[farm_ID]][,3], inclx=T)  

# Add new ts to dates
temp1 <- cbind(sub_temperature_exposure_wheat, temperatures_ts_M3)

# Aggregation to yearly values
final_wheat_temperature_aggregated_M3 <-temp1 %>%
  group_by(Farm,year)%>%
  summarise_at(vars(x,V2),sum)
rm(temp1)

# Combine yearly values with yields
temp1 <- wheat_yield_melted[which(wheat_yield_melted$Farm != farm_ID),]
data_wheat_final_M3 <- join(temp1, final_wheat_temperature_aggregated_M3, by=c("Farm","year"))
rm(temperatures_ts_M3,temp1,final_wheat_temperature_aggregated_M3)

# New TS for temperature support of interest
marginal_temperatures_M3 <- rcspline.eval(seq(-10,39,0.225),knots=list_knot_locations_kn3_wheat[[farm_ID]][,3], inclx=T)

temp_reg1 <- feols(yield ~ x + V2 + year + I(year^2) | Farm, data=data_wheat_final_M3)
summary(temp_reg1)
rm(temp_reg1)

# 
# Statistical uncertainty
# 

# Run 1000 models to derive uncertainty in estimates
spline_fitted_M3 <- matrix(NA,nrow=1000, ncol=nrow(marginal_temperatures_M3))

for (l in 1: 1000){
  # Sample some years
  data_temp <- data_wheat_final_M3[which(data_wheat_final_M3$year %in% sample(unique(data_wheat_final_M3$year),replace=T)),]
  data_temp$year2 <- (data_temp$year)^2 
  
  # Run regression with random subsample (data_temp)
  fit_boot  <- lm( yield ~ x + V2  + year +I(year^2)+ as.factor(Farm)-1, data=data_temp)
  
  # Save marginal effect of temperatures 
  spline_fitted_M3[l,]<-as.vector((marginal_temperatures_M3 %*% coef(fit_boot)[1:2]))
  
  # Show progress in %
  print(l/1000 *100)
  rm(fit_boot,data_temp)
}

# Change structure
spline_sub_melt_M3 <- melt(spline_fitted_M3)

# Get median response for each column (= a temperature from seq(3,39,0.225))
colmedians_spline_M3 <- colMedians(spline_fitted_M3)

# Get distance from median response for each value
dist_colmedian_M3 <- matrix(NA, nrow=nrow(spline_fitted_M3),ncol=ncol(spline_fitted_M3))

for(l in 1:ncol(spline_fitted_M3)){
  dist_colmedian_M3[,l] <- spline_fitted_M3[,l] - colmedians_spline_M3[l]
}

# Normalize distance to median
dist_colmedian_M3_normalized <- matrix(NA, nrow=nrow(spline_fitted_M3),ncol=ncol(spline_fitted_M3))

for(l in 1:ncol(spline_fitted_M3)){
  dist_colmedian_M3_normalized[,l]<-abs(dist_colmedian_M3[,l]/max(abs(dist_colmedian_M3[,l])))
}

# Change structure
dist_colmedian_M3_normalized_melted <- melt(dist_colmedian_M3_normalized)

# Last preparation before plotting

spline_merged_M3 <- merge(spline_sub_melt_M3,dist_colmedian_M3_normalized_melted,by=c("Var1","Var2"))
spline_merged_M3$value.y_minus1 <- 1-spline_merged_M3$value.y

temp1 <- aggregate(spline_merged_M3$value.y_minus1, list(spline_merged_M3$Var1), mean, na.rm=T) 
temp2 <- aggregate(spline_merged_M3$value.x, list(spline_merged_M3$Var2), mean, na.rm=T)

temp3 <- aggregate(spline_merged_M3$value.x, list(spline_merged_M3$Var2), quantile, probs=0.025) 
temp4 <- aggregate(spline_merged_M3$value.x, list(spline_merged_M3$Var2), quantile, probs=0.975)

spline_merged_M3 <- merge(spline_merged_M3, temp1, by.x=c("Var1"), by.y=c("Group.1"))
spline_merged_M3 <- merge(spline_merged_M3, temp2, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M3)[7]<-"average_effect"

spline_merged_M3<-merge(spline_merged_M3, temp3, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M3)[8]<-"lower_bound"

spline_merged_M3<-merge(spline_merged_M3, temp4, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M3)[9]<-"upper_bound"

temp5 <- cbind(seq(-10,39,0.225),seq(1:length(seq(-10,39,0.225))))
colnames(temp5) <- c("Temp_bin", "Var2")
spline_merged_M3 <-merge(spline_merged_M3, temp5, by=c("Var2"))

# 
# Plot: hourly temperature effect
# 

model_3 <- ggplot()+ ggtitle("Hourly temperature effects: Model 3") + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "grey25", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "grey25"),
        axis.text = element_text(family="Times New Roman", size=10),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  xlab("Temperature (°C)") + ylab("Hourly yield response \n (dt/ha)")+
  scale_x_continuous(breaks = seq(-5,39, by=5), limits = c(-5,40))+
  scale_y_continuous(breaks = c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1), limits = c(-0.4,0.1))+
  geom_hline(yintercept=c(seq(-0.5,0.1,by=0.1)), colour="gray")+
  geom_vline(xintercept=list_knot_locations_kn3_wheat[[farm_ID]][,3], linetype="dotted", colour="gray")+
  geom_line(data=spline_merged_M3, aes(y=value.x, x=Temp_bin, group=Var1, alpha=x.x^200/200), size=.01, colour="brown1")+
  geom_line(data=spline_merged_M3, aes(y=average_effect, x=Temp_bin), alpha=0.85, size=.75, colour="brown4")+
  geom_line(data=spline_merged_M3, aes(y=lower_bound, x=Temp_bin)   , alpha=0.85, size=.75, colour="brown4")+       
  geom_line(data=spline_merged_M3, aes(y=upper_bound, x=Temp_bin)   , alpha=0.85, size=.75, colour="brown4")  

rm(temp1,temp2,temp3,temp4,temp5)
rm(dist_colmedian_M3_normalized_melted, colmedians_spline_M3, spline_fitted_M3,dist_colmedian_M3_normalized,dist_colmedian_M3,spline_sub_melt_M3)

# 
# Preparation of hourly payouts
# 

# We change the sign of the effect (payouts are a positive number)

spline_merged_M3$average_effect <- spline_merged_M3$average_effect * (-1) 
spline_merged_M3$lower_bound <- spline_merged_M3$lower_bound * (-1) 
spline_merged_M3$upper_bound <- spline_merged_M3$upper_bound * (-1) 
spline_merged_M3$value.x <- spline_merged_M3$value.x * (-1) 

# Prepare lower bound (no values below 0)
temp <- subset(spline_merged_M3[,c(2,3,6,7,8,9,10)])
colnames(temp)[5] <- "upper_bound"
colnames(temp)[6] <- "lower_bound"

# Remove negative values(required for geom_ribbon())

temp[which(temp$lower_bound<0),6] <- 0

# 
# Plot: hourly payout function
# 

payout_wheat_model_3 <- ggplot(temp) +
  ggtitle("Hourly payouts: Model 3")+ theme(panel.grid.major = element_blank(),
                                                    plot.title = element_text(family="Times New Roman", size=12, colour = "grey25", hjust=0.5),
                                                    axis.title = element_text(family="Times New Roman", size=10, colour = "grey25"),
                                                    axis.text = element_text(family="Times New Roman", size=10),
                                                    panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(),
                                                    axis.line = element_line(colour = "black"),
                                                    legend.position="none")+
  xlab("Temperature (°C)") + ylab("Payout (dt/ha)")+
  scale_x_continuous(breaks = seq(10,39, by=5), limits = c(10,40))+
  scale_y_continuous(breaks = c(0,0.1, 0.2,0.3,0.4), limits = c(0,0.46))+
  geom_hline(yintercept=c(seq(0,0.4,by=0.1)), colour="gray")+
  geom_hline(yintercept=0, colour="darkgray")+
  geom_vline(xintercept= seq(5,38, by=5), colour="gray", linetype="dashed")+
  
  # Add lines
  
  geom_line(data=temp, aes(y=value.x, x=Temp_bin, group=Var1, alpha=x.x^200/200), size=.01, colour="royalblue1")+
  geom_line(data=temp, aes(y=average_effect, x=Temp_bin), size=0.75, color="royalblue4")+
  geom_line(data=temp, aes(y=lower_bound, x=Temp_bin), size=0.75, color="royalblue4")+
  geom_line(data=temp, aes(y=upper_bound, x=Temp_bin), size=0.75, color="royalblue4")

rm(temp,spline_merged_M3)


# -------------------------------------------
# Model 4: 
# -------------------------------------------

# Temperature time series
temperatures_ts_M4 <- rcspline.eval(sub_temperature_exposure_wheat$Temperature,knots= c(5,10,15,20,25), inclx=T)  

# Add new ts to dates
temp1 <- cbind(sub_temperature_exposure_wheat, temperatures_ts_M4)

# Aggregation to yearly values
final_wheat_temperature_aggregated_M4 <-temp1 %>%
  group_by(Farm,year)%>%
  summarise_at(vars(x,V2,V3,V4),sum)
rm(temp1)

# Combine yearly values with yields
temp1 <- wheat_yield_melted[which(wheat_yield_melted$Farm != farm_ID),]
data_wheat_final_M4 <- join(temp1, final_wheat_temperature_aggregated_M4, by=c("Farm","year"))
rm(temperatures_ts_M4,temp1,final_wheat_temperature_aggregated_M4)

# New TS for temperature support of interest
marginal_temperatures_M4 <- rcspline.eval(seq(-10,39,0.225),knots=knot_locations_kn5_wheat[,1], inclx=T)

temp_reg1 <- feols(yield ~ x + V2 +V3 +V4 + year + I(year^2) | Farm, data=data_wheat_final_M4)
summary(temp_reg1)
rm(temp_reg1)

# -------------------------------------------
# Statistical uncertainty
# -------------------------------------------

# Run 1000 models to derive uncertainty in estimates
spline_fitted_M4 <- matrix(NA,nrow=1000, ncol=nrow(marginal_temperatures_M4))

for (l in 1: 1000){
  # Sample some years
  data_temp <- data_wheat_final_M4[which(data_wheat_final_M4$year %in% sample(unique(data_wheat_final_M4$year),replace=T)),]
  
  # Run regression with random subsample (data_temp)
  fit_boot  <- lm( yield ~ x + V2 +V3 +V4 + year + I(year^2)  + as.factor(Farm) -1, data=data_temp)
  
  # Save marginal effect of temperatures 
  spline_fitted_M4[l,]<-as.vector((marginal_temperatures_M4 %*% coef(fit_boot)[1:4]))
  
  # Show progress in %
  print(l/1000 *100)
  rm(fit_boot,data_temp)
}

# Change structure
spline_sub_melt_M4 <- melt(spline_fitted_M4)

# Get median response for each column (= a temperature from seq(-10,39,0.225))
colmedians_spline_M4 <- colMedians(spline_fitted_M4)

# Get distance from median response for each value
dist_colmedian_M4 <- matrix(NA, nrow=nrow(spline_fitted_M4),ncol=ncol(spline_fitted_M4))

for(l in 1:ncol(spline_fitted_M4)){
  dist_colmedian_M4[,l] <- spline_fitted_M4[,l] - colmedians_spline_M4[l]
}

# Normalize distance to median
# 1= largest value; 0= no distance
dist_colmedian_M4_normalized <- matrix(NA, nrow=nrow(spline_fitted_M4),ncol=ncol(spline_fitted_M4))

for(l in 1:ncol(spline_fitted_M4)){
  dist_colmedian_M4_normalized[,l]<-abs(dist_colmedian_M4[,l]/max(abs(dist_colmedian_M4[,l])))
}

# Change structure
dist_colmedian_M4_normalized_melted <- melt(dist_colmedian_M4_normalized)

# Last preparation before plotting
# Var1 = 1:1000; Var2= temperature
# value.y ist normalized distance to median; value.x is estimated effect
spline_merged_M4 <- merge(spline_sub_melt_M4,dist_colmedian_M4_normalized_melted,by=c("Var1","Var2"))
spline_merged_M4$value.y_minus1 <- 1-spline_merged_M4$value.y

temp1 <- aggregate(spline_merged_M4$value.y_minus1, list(spline_merged_M4$Var1), mean, na.rm=T) 
temp2 <- aggregate(spline_merged_M4$value.x, list(spline_merged_M4$Var2), mean, na.rm=T)

temp3 <- aggregate(spline_merged_M4$value.x, list(spline_merged_M4$Var2), quantile, probs=0.025) 
temp4 <- aggregate(spline_merged_M4$value.x, list(spline_merged_M4$Var2), quantile, probs=0.975)

spline_merged_M4 <- merge(spline_merged_M4, temp1, by.x=c("Var1"), by.y=c("Group.1"))
spline_merged_M4 <- merge(spline_merged_M4, temp2, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M4)[7]<-"average_effect"

spline_merged_M4<-merge(spline_merged_M4, temp3, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M4)[8]<-"lower_bound"

spline_merged_M4<-merge(spline_merged_M4, temp4, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M4)[9]<-"upper_bound"

temp5 <- cbind(seq(-10,39,0.225),seq(1:length(seq(-10,39,0.225))))
colnames(temp5) <- c("Temp_bin", "Var2")
spline_merged_M4 <-merge(spline_merged_M4, temp5, by=c("Var2"))

# 
# Plot: hourly temperature effect
# 

model_4 <- ggplot()+ ggtitle("Hourly temperature effects: Model 4") + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "grey25", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "grey25"),
        axis.text = element_text(family="Times New Roman", size=10),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  xlab("Temperature (°C)") + ylab("Hourly yield response \n (dt/ha)")+
  scale_x_continuous(breaks = seq(-5,39, by=5), limits = c(-5,40))+
  scale_y_continuous(breaks = c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1), limits = c(-0.4,0.1))+
  geom_hline(yintercept=c(seq(-0.5,0.1,by=0.1)), colour="gray")+
  geom_vline(xintercept=knot_locations_kn5_wheat[,1], linetype="dotted", colour="gray")+
  geom_line(data=spline_merged_M4, aes(y=value.x, x=Temp_bin, group=Var1, alpha=x.x^200/200), size=.01, colour="brown1")+
  geom_line(data=spline_merged_M4, aes(y=average_effect, x=Temp_bin), alpha=0.85, size=.75, colour="brown4")+
  geom_line(data=spline_merged_M4, aes(y=lower_bound, x=Temp_bin)   , alpha=0.85, size=.75, colour="brown4")+       
  geom_line(data=spline_merged_M4, aes(y=upper_bound, x=Temp_bin)   , alpha=0.85, size=.75, colour="brown4")  

rm(temp1,temp2,temp3,temp4,temp5)
rm(dist_colmedian_M4_normalized_melted, colmedians_spline_M4, spline_fitted_M4,dist_colmedian_M4_normalized,dist_colmedian_M4,spline_sub_melt_M4)

# 
# Preparation of hourly payouts
# 

# We change the sign of the effect (payouts are a positive number)

spline_merged_M4$average_effect <- spline_merged_M4$average_effect * (-1) 
spline_merged_M4$lower_bound <- spline_merged_M4$lower_bound * (-1) 
spline_merged_M4$upper_bound <- spline_merged_M4$upper_bound * (-1) 
spline_merged_M4$value.x <- spline_merged_M4$value.x * (-1) 

# Prepare lower bound (no values below 0)
temp <- subset(spline_merged_M4[,c(2,3,6,7,8,9,10)])
colnames(temp)[5] <- "upper_bound"
colnames(temp)[6] <- "lower_bound"

# Remove negative values(required for geom_ribbon())

temp[which(temp$lower_bound<0),6] <- 0

# 
# Plot: hourly payout function
# 

payout_wheat_model_4 <- ggplot(temp) +
  ggtitle("Hourly payouts: Model 4")+ theme(panel.grid.major = element_blank(),
                                                    plot.title = element_text(family="Times New Roman", size=12, colour = "grey25", hjust=0.5),
                                                    axis.title = element_text(family="Times New Roman", size=10, colour = "grey25"),
                                                    axis.text = element_text(family="Times New Roman", size=10),
                                                    panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(),
                                                    axis.line = element_line(colour = "black"),
                                                    legend.position="none")+
  xlab("Temperature (°C)") + ylab("Payout (dt/ha)")+
  scale_x_continuous(breaks = seq(10,39, by=5), limits = c(10,40))+
  scale_y_continuous(breaks = c(0,0.1, 0.2,0.3,0.4), limits = c(0,0.46))+
  geom_hline(yintercept=c(seq(0,0.4,by=0.1)), colour="gray")+
  geom_hline(yintercept=0, colour="darkgray")+
  geom_vline(xintercept= seq(5,38, by=5), colour="gray", linetype="dashed")+
  
  # Add lines
  
  geom_line(data=temp, aes(y=value.x, x=Temp_bin, group=Var1, alpha=x.x^200/200), size=.01, colour="royalblue1")+
  geom_line(data=temp, aes(y=average_effect, x=Temp_bin), size=0.75, color="royalblue4")+
  geom_line(data=temp, aes(y=lower_bound, x=Temp_bin), size=0.75, color="royalblue4")+
  geom_line(data=temp, aes(y=upper_bound, x=Temp_bin), size=0.75, color="royalblue4")

rm(temp,spline_merged_M4)


models_payouts_wheat <- ggarrange(model_1,payout_wheat_model_1,model_2,payout_wheat_model_2,model_3,payout_wheat_model_3,model_4,payout_wheat_model_4,temperature_wheat, ncol=2, byrow=T, heights = c(1,1,1,1,0.5))
#ggsave(plot=models_payouts_wheat, file="Plots/T Effect/models_payouts_wheat.png", width = 16, height = 20, unit="cm")
rm(model_1,payout_wheat_model_1,model_2,payout_wheat_model_2,model_3,payout_wheat_model_3,model_4,payout_wheat_model_4,temperature_wheat)

# =======================================================================
#
# 2) Figure 3: hourly temperature effects and payouts for winter rapeseed
#
# =======================================================================

# The plots are representative for a farm drawn at random.
# Note that plots slightly vary between farms because of the out-of-sample calibration
set.seed(123)
farm_ID <- round(runif(1,1,nrow(farm_yields[[2]])))

# -------------------------------------------
# Histogram of temperature exposure
# -------------------------------------------

sub_temperature_exposure_rapeseed <- temp_rapeseed_exposure2[which(temp_rapeseed_exposure2$Farm != farm_ID),]

temperature_rapeseed <- ggplot(sub_temperature_exposure_rapeseed, aes(x=Temperature)) + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "grey25", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "grey25"),
        axis.text = element_text(family="Times New Roman", size=10),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  ylab("Total exposure \n (1'000 hours)")+ xlab("Temperature (°C)")+
  geom_histogram(binwidth=1, color="black", fill="brown1")+
  scale_x_continuous(breaks = seq(-5,39, by=5), limits = c(-5,40))+
  scale_y_continuous(breaks = c(0,100000,200000), labels = c(0,100,200))

# -------------------------------------------
# Model 1: 
# -------------------------------------------

# Temperature time series
temperatures_ts_M1 <- rcspline.eval(sub_temperature_exposure_rapeseed$Temperature,knots= list_knot_locations_kn3_rapeseed[[farm_ID]][,1], inclx=T)  

# Add new ts to dates
temp1 <- cbind(sub_temperature_exposure_rapeseed, temperatures_ts_M1)

# Aggregation to yearly values
final_rapeseed_temperature_aggregated_M1 <-temp1 %>%
  group_by(Farm,year)%>%
  summarise_at(vars(x,V2),sum)
rm(temp1)

# Combine yearly values with yields
temp1 <- rapeseed_yield_melted[which(rapeseed_yield_melted$Farm != farm_ID),]
data_rapeseed_final_M1 <- join(temp1, final_rapeseed_temperature_aggregated_M1, by=c("Farm","year"))
rm(temperatures_ts_M1,temp1,final_rapeseed_temperature_aggregated_M1)

# New TS for temperature support of interest
marginal_temperatures_M1 <- rcspline.eval(seq(-10,39,0.225),knots=list_knot_locations_kn3_rapeseed[[farm_ID]][,1], inclx=T)

# Report regression output
temp_reg1 <- feols(yield ~ x + V2  + year + I(year^2) | Farm, data=data_rapeseed_final_M1)
summary(temp_reg1)
rm(temp_reg1)

# -------------------------------------------
# Statistical uncertainty
# -------------------------------------------

# Run 1000 models to derive uncertainty in estimates
spline_fitted_M1 <- matrix(NA,nrow=1000, ncol=nrow(marginal_temperatures_M1))

for (l in 1: 1000){
  # Sample some years
  data_temp <- data_rapeseed_final_M1[which(data_rapeseed_final_M1$year %in% sample(unique(data_rapeseed_final_M1$year),replace=T)),]
  
  # Run regression with random subsample (data_temp)
  fit_boot  <- lm( yield ~ x + V2 +year + I(year^2) + as.factor(Farm)-1, data=data_temp)
  
  # Save marginal effect of temperatures 
  spline_fitted_M1[l,]<-as.vector((marginal_temperatures_M1 %*% coef(fit_boot)[1:2]))
  
  # Show progress in %
  print(l/1000 *100)
  rm(fit_boot,data_temp)
}

# Change structure
spline_sub_melt_M1 <- melt(spline_fitted_M1)

# Get median response for each column (= a temperature from seq(3,39,0.225))
colmedians_spline_M1 <- colMedians(spline_fitted_M1)

# Get distance from median response for each value
dist_colmedian_M1 <- matrix(NA, nrow=nrow(spline_fitted_M1),ncol=ncol(spline_fitted_M1))

for(l in 1:ncol(spline_fitted_M1)){
  dist_colmedian_M1[,l] <- spline_fitted_M1[,l] - colmedians_spline_M1[l]
}

# Normalize distance to median
dist_colmedian_M1_normalized <- matrix(NA, nrow=nrow(spline_fitted_M1),ncol=ncol(spline_fitted_M1))

for(l in 1:ncol(spline_fitted_M1)){
  dist_colmedian_M1_normalized[,l]<-abs(dist_colmedian_M1[,l]/max(abs(dist_colmedian_M1[,l])))
}

# Change structure
dist_colmedian_M1_normalized_melted <- melt(dist_colmedian_M1_normalized)

# Last preparation before plotting

spline_merged_M1 <- merge(spline_sub_melt_M1,dist_colmedian_M1_normalized_melted,by=c("Var1","Var2"))
spline_merged_M1$value.y_minus1 <- 1-spline_merged_M1$value.y

temp1 <- aggregate(spline_merged_M1$value.y_minus1, list(spline_merged_M1$Var1), mean, na.rm=T) 
temp2 <- aggregate(spline_merged_M1$value.x, list(spline_merged_M1$Var2), mean, na.rm=T)

temp3 <- aggregate(spline_merged_M1$value.x, list(spline_merged_M1$Var2), quantile, probs=0.025) 
temp4 <- aggregate(spline_merged_M1$value.x, list(spline_merged_M1$Var2), quantile, probs=0.975)

spline_merged_M1 <- merge(spline_merged_M1, temp1, by.x=c("Var1"), by.y=c("Group.1"))
spline_merged_M1 <- merge(spline_merged_M1, temp2, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M1)[7]<-"average_effect"

spline_merged_M1<-merge(spline_merged_M1, temp3, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M1)[8]<-"lower_bound"

spline_merged_M1<-merge(spline_merged_M1, temp4, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M1)[9]<-"upper_bound"

temp5 <- cbind(seq(-10,39,0.225),seq(1:length(seq(-10,39,0.225))))
colnames(temp5) <- c("Temp_bin", "Var2")
spline_merged_M1 <-merge(spline_merged_M1, temp5, by=c("Var2"))

# 
# Plot: hourly temperature effect
# 

model_1 <- ggplot()+ ggtitle("Hourly temperature effects: Model 1") + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "grey25", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "grey25"),
        axis.text = element_text(family="Times New Roman", size=10),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  xlab("Temperature (°C)") + ylab("Hourly yield response \n (dt/ha)")+
  scale_x_continuous(breaks = seq(-5,39, by=5), limits = c(-5,40))+
  scale_y_continuous(breaks = c(-0.2,-0.1,0,0.1), limits = c(-0.2,0.1))+
  geom_hline(yintercept=c(seq(-0.5,0.1,by=0.1)), colour="gray")+
  geom_vline(xintercept=list_knot_locations_kn3_rapeseed[[farm_ID]][,1], linetype="dotted", colour="gray")+
  geom_line(data=spline_merged_M1, aes(y=value.x, x=Temp_bin, group=Var1, alpha=x.x^200/200), size=.01, colour="brown1")+
  geom_line(data=spline_merged_M1, aes(y=average_effect, x=Temp_bin), alpha=0.85, size=.75, colour="brown4")+
  geom_line(data=spline_merged_M1, aes(y=lower_bound, x=Temp_bin)   , alpha=0.85, size=.75, colour="brown4")+       
  geom_line(data=spline_merged_M1, aes(y=upper_bound, x=Temp_bin)   , alpha=0.85, size=.75, colour="brown4")

rm(temp1,temp2,temp3,temp4,temp5)
rm(dist_colmedian_M1_normalized_melted, colmedians_spline_M1, spline_fitted_M1,dist_colmedian_M1_normalized,dist_colmedian_M1,spline_sub_melt_M1)

# 
# Preparation of hourly payouts
# 

# We change the sign of the effect (payouts are a positive number)

spline_merged_M1$average_effect <- spline_merged_M1$average_effect * (-1) 
spline_merged_M1$lower_bound <- spline_merged_M1$lower_bound * (-1) 
spline_merged_M1$upper_bound <- spline_merged_M1$upper_bound * (-1) 
spline_merged_M1$value.x <- spline_merged_M1$value.x * (-1) 

# Prepare lower bound (no values below 0)
temp <- subset(spline_merged_M1[,c(2,3,6,7,8,9,10)])
colnames(temp)[5] <- "upper_bound"
colnames(temp)[6] <- "lower_bound"

# Remove negative values(required for geom_ribbon())

temp[which(temp$lower_bound<0),6] <- 0

# 
# Plot: hourly payout function
# 

payout_rapeseed_model_1 <- ggplot(temp) +
  ggtitle("Hourly payouts: Model 1")+ theme(panel.grid.major = element_blank(),
                                                    plot.title = element_text(family="Times New Roman", size=12, colour = "grey25", hjust=0.5),
                                                    axis.title = element_text(family="Times New Roman", size=10, colour = "grey25"),
                                                    axis.text = element_text(family="Times New Roman", size=10),
                                                    panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(),
                                                    axis.line = element_line(colour = "black"),
                                                    legend.position="none")+
  xlab("Temperature (°C)") + ylab("Payout (dt/ha)")+
  scale_x_continuous(breaks = seq(10,39, by=5), limits = c(10,40))+
  scale_y_continuous(breaks = c(0,0.1, 0.2,0.3,0.4), limits = c(0,0.2))+
  geom_hline(yintercept=c(seq(0,0.4,by=0.1)), colour="gray")+
  geom_hline(yintercept=0, colour="darkgray")+
  geom_vline(xintercept= seq(5,38, by=5), colour="gray", linetype="dashed")+
  
  # Add lines
  
  geom_line(data=temp, aes(y=value.x, x=Temp_bin, group=Var1, alpha=x.x^200/200), size=.01, colour="royalblue1")+
  geom_line(data=temp, aes(y=average_effect, x=Temp_bin), size=0.75, color="royalblue4")+
  geom_line(data=temp, aes(y=lower_bound, x=Temp_bin), size=0.75, color="royalblue4")+
  geom_line(data=temp, aes(y=upper_bound, x=Temp_bin), size=0.75, color="royalblue4")

rm(temp,spline_merged_M1)

# -------------------------------------------
# Model 2: 
# -------------------------------------------

# Temperature time series
temperatures_ts_M2 <- rcspline.eval(sub_temperature_exposure_rapeseed$Temperature,knots= list_knot_locations_kn3_rapeseed[[farm_ID]][,2], inclx=T)  

# Add new ts to dates
temp1 <- cbind(sub_temperature_exposure_rapeseed, temperatures_ts_M2)

# Aggregation to yearly values
final_rapeseed_temperature_aggregated_M2 <-temp1 %>%
  group_by(Farm,year)%>%
  summarise_at(vars(x,V2),sum)
rm(temp1)

# Combine yearly values with yields
temp1 <- rapeseed_yield_melted[which(rapeseed_yield_melted$Farm != farm_ID),]
data_rapeseed_final_M2 <- join(temp1, final_rapeseed_temperature_aggregated_M2, by=c("Farm","year"))
rm(temperatures_ts_M2,temp1,final_rapeseed_temperature_aggregated_M2)

# New TS for temperature support of interest
marginal_temperatures_M2 <- rcspline.eval(seq(-10,39,0.225),knots=list_knot_locations_kn3_rapeseed[[farm_ID]][,2], inclx=T)

# Report regression output
temp_reg1 <- feols(yield ~ x + V2  + year + I(year^2) | Farm, data=data_rapeseed_final_M2)
summary(temp_reg1)
rm(temp_reg1)

# -------------------------------------------
# Statistical uncertainty
# -------------------------------------------

# Run 1000 models to derive uncertainty in estimates
spline_fitted_M2 <- matrix(NA,nrow=1000, ncol=nrow(marginal_temperatures_M2))

for (l in 1: 1000){
  # Sample some years
  data_temp <- data_rapeseed_final_M2[which(data_rapeseed_final_M2$year %in% sample(unique(data_rapeseed_final_M2$year),replace=T)),]
  
  # Run regression with random subsample (data_temp)
  fit_boot  <- lm( yield ~ x + V2 +year + I(year^2) + as.factor(Farm)-1, data=data_temp)
  
  # Save marginal effect of temperatures 
  spline_fitted_M2[l,]<-as.vector((marginal_temperatures_M2 %*% coef(fit_boot)[1:2]))
  
  # Show progress in %
  print(l/1000 *100)
  rm(fit_boot,data_temp)
}

# Change structure
spline_sub_melt_M2 <- melt(spline_fitted_M2)

# Get median response for each column (= a temperature from seq(3,39,0.225))
colmedians_spline_M2 <- colMedians(spline_fitted_M2)

# Get distance from median response for each value
dist_colmedian_M2 <- matrix(NA, nrow=nrow(spline_fitted_M2),ncol=ncol(spline_fitted_M2))

for(l in 1:ncol(spline_fitted_M2)){
  dist_colmedian_M2[,l] <- spline_fitted_M2[,l] - colmedians_spline_M2[l]
}

# Normalize distance to median
dist_colmedian_M2_normalized <- matrix(NA, nrow=nrow(spline_fitted_M2),ncol=ncol(spline_fitted_M2))

for(l in 1:ncol(spline_fitted_M2)){
  dist_colmedian_M2_normalized[,l]<-abs(dist_colmedian_M2[,l]/max(abs(dist_colmedian_M2[,l])))
}

# Change structure
dist_colmedian_M2_normalized_melted <- melt(dist_colmedian_M2_normalized)

# Last preparation before plotting

spline_merged_M2 <- merge(spline_sub_melt_M2,dist_colmedian_M2_normalized_melted,by=c("Var1","Var2"))
spline_merged_M2$value.y_minus1 <- 1-spline_merged_M2$value.y

temp1 <- aggregate(spline_merged_M2$value.y_minus1, list(spline_merged_M2$Var1), mean, na.rm=T) 
temp2 <- aggregate(spline_merged_M2$value.x, list(spline_merged_M2$Var2), mean, na.rm=T)

temp3 <- aggregate(spline_merged_M2$value.x, list(spline_merged_M2$Var2), quantile, probs=0.025) 
temp4 <- aggregate(spline_merged_M2$value.x, list(spline_merged_M2$Var2), quantile, probs=0.975)

spline_merged_M2 <- merge(spline_merged_M2, temp1, by.x=c("Var1"), by.y=c("Group.1"))
spline_merged_M2 <- merge(spline_merged_M2, temp2, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M2)[7]<-"average_effect"

spline_merged_M2<-merge(spline_merged_M2, temp3, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M2)[8]<-"lower_bound"

spline_merged_M2<-merge(spline_merged_M2, temp4, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M2)[9]<-"upper_bound"

temp5 <- cbind(seq(-10,39,0.225),seq(1:length(seq(-10,39,0.225))))
colnames(temp5) <- c("Temp_bin", "Var2")
spline_merged_M2 <-merge(spline_merged_M2, temp5, by=c("Var2"))

# 
# Plot: hourly temperature effect
# 

model_2 <- ggplot()+ ggtitle("Hourly temperature effects: Model 2") + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "grey25", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "grey25"),
        axis.text = element_text(family="Times New Roman", size=10),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  xlab("Temperature (°C)") + ylab("Hourly yield response \n (dt/ha)")+
  scale_x_continuous(breaks = seq(-5,39, by=5), limits = c(-5,40))+
  scale_y_continuous(breaks = c(-0.2,-0.1,0,0.1), limits = c(-0.2,0.1))+
  geom_hline(yintercept=c(seq(-0.5,0.1,by=0.1)), colour="gray")+
  geom_vline(xintercept=list_knot_locations_kn3_rapeseed[[farm_ID]][,2], linetype="dotted", colour="gray")+
  geom_line(data=spline_merged_M2, aes(y=value.x, x=Temp_bin, group=Var1, alpha=x.x^200/200), size=.01, colour="brown1")+
  geom_line(data=spline_merged_M2, aes(y=average_effect, x=Temp_bin), alpha=0.85, size=.75, colour="brown4")+
  geom_line(data=spline_merged_M2, aes(y=lower_bound, x=Temp_bin)   , alpha=0.85, size=.75, colour="brown4")+       
  geom_line(data=spline_merged_M2, aes(y=upper_bound, x=Temp_bin)   , alpha=0.85, size=.75, colour="brown4")

rm(temp1,temp2,temp3,temp4,temp5)
rm(dist_colmedian_M2_normalized_melted, colmedians_spline_M2, spline_fitted_M2,dist_colmedian_M2_normalized,dist_colmedian_M2,spline_sub_melt_M2)

# 
# Preparation of hourly payouts
# 

# We change the sign of the effect (payouts are a positive number)

spline_merged_M2$average_effect <- spline_merged_M2$average_effect * (-1) 
spline_merged_M2$lower_bound <- spline_merged_M2$lower_bound * (-1) 
spline_merged_M2$upper_bound <- spline_merged_M2$upper_bound * (-1) 
spline_merged_M2$value.x <- spline_merged_M2$value.x * (-1) 

# Prepare lower bound (no values below 0)
temp <- subset(spline_merged_M2[,c(2,3,6,7,8,9,10)])
colnames(temp)[5] <- "upper_bound"
colnames(temp)[6] <- "lower_bound"

# Remove negative values(required for geom_ribbon())

temp[which(temp$lower_bound<0),6] <- 0

# 
# Plot: hourly payout function
# 

payout_rapeseed_model_2 <- ggplot(temp) +
  ggtitle("Hourly payouts: Model 2")+ theme(panel.grid.major = element_blank(),
                                                    plot.title = element_text(family="Times New Roman", size=12, colour = "grey25", hjust=0.5),
                                                    axis.title = element_text(family="Times New Roman", size=10, colour = "grey25"),
                                                    axis.text = element_text(family="Times New Roman", size=10),
                                                    panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(),
                                                    axis.line = element_line(colour = "black"),
                                                    legend.position="none")+
  xlab("Temperature (°C)") + ylab("Payout (dt/ha)")+
  scale_x_continuous(breaks = seq(10,39, by=5), limits = c(10,40))+
  scale_y_continuous(breaks = c(0,0.1, 0.2), limits = c(0,0.2))+
  geom_hline(yintercept=c(seq(0,0.4,by=0.1)), colour="gray")+
  geom_hline(yintercept=0, colour="darkgray")+
  geom_vline(xintercept= seq(5,38, by=5), colour="gray", linetype="dashed")+
  
  # Add lines
  
  geom_line(data=temp, aes(y=value.x, x=Temp_bin, group=Var1, alpha=x.x^200/200), size=.01, colour="royalblue1")+
  geom_line(data=temp, aes(y=average_effect, x=Temp_bin), size=0.75, color="royalblue4")+
  geom_line(data=temp, aes(y=lower_bound, x=Temp_bin), size=0.75, color="royalblue4")+
  geom_line(data=temp, aes(y=upper_bound, x=Temp_bin), size=0.75, color="royalblue4")

rm(temp,spline_merged_M2)

# -------------------------------------------
# Model 3: 
# -------------------------------------------

# Temperature time series
temperatures_ts_M3 <- rcspline.eval(sub_temperature_exposure_rapeseed$Temperature,knots= list_knot_locations_kn3_rapeseed[[farm_ID]][,3], inclx=T)  

# Add new ts to dates
temp1 <- cbind(sub_temperature_exposure_rapeseed, temperatures_ts_M3)

# Aggregation to yearly values
final_rapeseed_temperature_aggregated_M3 <-temp1 %>%
  group_by(Farm,year)%>%
  summarise_at(vars(x,V2),sum)
rm(temp1)

# Combine yearly values with yields
temp1 <- rapeseed_yield_melted[which(rapeseed_yield_melted$Farm != farm_ID),]
data_rapeseed_final_M3 <- join(temp1, final_rapeseed_temperature_aggregated_M3, by=c("Farm","year"))
rm(temperatures_ts_M3,temp1,final_rapeseed_temperature_aggregated_M3)

# New TS for temperature support of interest
marginal_temperatures_M3 <- rcspline.eval(seq(-10,39,0.225),knots=list_knot_locations_kn3_rapeseed[[farm_ID]][,3], inclx=T)

temp_reg1 <- feols(yield ~ x + V2  + year + I(year^2) | Farm, data=data_rapeseed_final_M3)
summary(temp_reg1)
rm(temp_reg1)

# -------------------------------------------
# Statistical uncertainty 
# -------------------------------------------

# Run 1000 models to derive uncertainty in estimates
spline_fitted_M3 <- matrix(NA,nrow=1000, ncol=nrow(marginal_temperatures_M3))

for (l in 1: 1000){
  # Sample some years
  data_temp <- data_rapeseed_final_M3[which(data_rapeseed_final_M3$year %in% sample(unique(data_rapeseed_final_M3$year),replace=T)),]
  
  # Run regression with random subsample (data_temp)
  fit_boot  <- lm( yield ~ x + V2 +year + I(year^2) + as.factor(Farm)-1, data=data_temp)
  
  # Save marginal effect of temperatures 
  spline_fitted_M3[l,]<-as.vector((marginal_temperatures_M3 %*% coef(fit_boot)[1:2]))
  
  # Show progress in %
  print(l/1000 *100)
  rm(fit_boot,data_temp)
}

# Change structure
spline_sub_melt_M3 <- melt(spline_fitted_M3)


# Get median response for each column (= a temperature from seq(3,39,0.225))
colmedians_spline_M3 <- colMedians(spline_fitted_M3)

# Get distance from median response for each value
dist_colmedian_M3 <- matrix(NA, nrow=nrow(spline_fitted_M3),ncol=ncol(spline_fitted_M3))

for(l in 1:ncol(spline_fitted_M3)){
  dist_colmedian_M3[,l] <- spline_fitted_M3[,l] - colmedians_spline_M3[l]
}

# Normalize distance to median
dist_colmedian_M3_normalized <- matrix(NA, nrow=nrow(spline_fitted_M3),ncol=ncol(spline_fitted_M3))

for(l in 1:ncol(spline_fitted_M3)){
  dist_colmedian_M3_normalized[,l]<-abs(dist_colmedian_M3[,l]/max(abs(dist_colmedian_M3[,l])))
}

# Change structure
dist_colmedian_M3_normalized_melted <- melt(dist_colmedian_M3_normalized)

# Last preparation before plotting

spline_merged_M3 <- merge(spline_sub_melt_M3,dist_colmedian_M3_normalized_melted,by=c("Var1","Var2"))
spline_merged_M3$value.y_minus1 <- 1-spline_merged_M3$value.y

temp1 <- aggregate(spline_merged_M3$value.y_minus1, list(spline_merged_M3$Var1), mean, na.rm=T) 
temp2 <- aggregate(spline_merged_M3$value.x, list(spline_merged_M3$Var2), mean, na.rm=T)

temp3 <- aggregate(spline_merged_M3$value.x, list(spline_merged_M3$Var2), quantile, probs=0.025) 
temp4 <- aggregate(spline_merged_M3$value.x, list(spline_merged_M3$Var2), quantile, probs=0.975)

spline_merged_M3 <- merge(spline_merged_M3, temp1, by.x=c("Var1"), by.y=c("Group.1"))
spline_merged_M3 <- merge(spline_merged_M3, temp2, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M3)[7]<-"average_effect"

spline_merged_M3<-merge(spline_merged_M3, temp3, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M3)[8]<-"lower_bound"

spline_merged_M3<-merge(spline_merged_M3, temp4, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M3)[9]<-"upper_bound"

temp5 <- cbind(seq(-10,39,0.225),seq(1:length(seq(-10,39,0.225))))
colnames(temp5) <- c("Temp_bin", "Var2")
spline_merged_M3 <-merge(spline_merged_M3, temp5, by=c("Var2"))

# 
# Plot: hourly temperature effect
#

model_3 <- ggplot()+ ggtitle("Hourly temperature effects: Model 3") + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "grey25", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "grey25"),
        axis.text = element_text(family="Times New Roman", size=10),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  xlab("Temperature (°C)") + ylab("Hourly yield response \n (dt/ha)")+
  scale_x_continuous(breaks = seq(-5,39, by=5), limits = c(-5,40))+
  scale_y_continuous(breaks = c(-0.2,-0.1,0,0.1), limits = c(-0.2,0.1))+
  geom_hline(yintercept=c(seq(-0.5,0.1,by=0.1)), colour="gray")+
  geom_vline(xintercept=list_knot_locations_kn3_rapeseed[[farm_ID]][,3], linetype="dotted", colour="gray")+
  geom_line(data=spline_merged_M3, aes(y=value.x, x=Temp_bin, group=Var1, alpha=x.x^200/200), size=.01, colour="brown1")+
  geom_line(data=spline_merged_M3, aes(y=average_effect, x=Temp_bin), alpha=0.85, size=.75, colour="brown4")+
  geom_line(data=spline_merged_M3, aes(y=lower_bound, x=Temp_bin)   , alpha=0.85, size=.75, colour="brown4")+       
  geom_line(data=spline_merged_M3, aes(y=upper_bound, x=Temp_bin)   , alpha=0.85, size=.75, colour="brown4")

rm(temp1,temp2,temp3,temp4,temp5)
rm(dist_colmedian_M3_normalized_melted, colmedians_spline_M3, spline_fitted_M3,dist_colmedian_M3_normalized,dist_colmedian_M3,spline_sub_melt_M3)

# 
# Preparation of hourly payouts
#

# We change the sign of the effect (payouts are a positive number)

spline_merged_M3$average_effect <- spline_merged_M3$average_effect * (-1) 
spline_merged_M3$lower_bound <- spline_merged_M3$lower_bound * (-1) 
spline_merged_M3$upper_bound <- spline_merged_M3$upper_bound * (-1) 
spline_merged_M3$value.x <- spline_merged_M3$value.x * (-1) 

# Prepare lower bound (no values below 0)
temp <- subset(spline_merged_M3[,c(2,3,6,7,8,9,10)])
colnames(temp)[5] <- "upper_bound"
colnames(temp)[6] <- "lower_bound"

# Remove negative values(required for geom_ribbon())

temp[which(temp$lower_bound<0),6] <- 0

#
# Plot: hourly payout function
#

payout_rapeseed_model_3 <- ggplot(temp) +
  ggtitle("Hourly payouts: Model 3")+ theme(panel.grid.major = element_blank(),
                                                    plot.title = element_text(family="Times New Roman", size=12, colour = "grey25", hjust=0.5),
                                                    axis.title = element_text(family="Times New Roman", size=10, colour = "grey25"),
                                                    axis.text = element_text(family="Times New Roman", size=10),
                                                    panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(),
                                                    axis.line = element_line(colour = "black"),
                                                    legend.position="none")+
  xlab("Temperature (°C)") + ylab("Payout (dt/ha)")+
  scale_x_continuous(breaks = seq(10,39, by=5), limits = c(10,40))+
  scale_y_continuous(breaks = c(0,0.1, 0.2), limits = c(0,0.2))+
  geom_hline(yintercept=c(seq(0,0.4,by=0.1)), colour="gray")+
  geom_hline(yintercept=0, colour="darkgray")+
  geom_vline(xintercept= seq(5,38, by=5), colour="gray", linetype="dashed")+
  
  # Add lines
  
  geom_line(data=temp, aes(y=value.x, x=Temp_bin, group=Var1, alpha=x.x^200/200), size=.01, colour="royalblue1")+
  geom_line(data=temp, aes(y=average_effect, x=Temp_bin), size=0.75, color="royalblue4")+
  geom_line(data=temp, aes(y=lower_bound, x=Temp_bin), size=0.75, color="royalblue4")+
  geom_line(data=temp, aes(y=upper_bound, x=Temp_bin), size=0.75, color="royalblue4")

rm(temp,spline_merged_M3)

# -------------------------------------------
# Model 4: 
# -------------------------------------------

# Temperature time series
temperatures_ts_M4 <- rcspline.eval(sub_temperature_exposure_rapeseed$Temperature,knots= c(5,10,15,20,25), inclx=T)  

# Add new ts to dates
temp1 <- cbind(sub_temperature_exposure_rapeseed, temperatures_ts_M4)

# Aggregation to yearly values
final_rapeseed_temperature_aggregated_M4 <-temp1 %>%
  group_by(Farm,year)%>%
  summarise_at(vars(x,V2,V3,V4),sum)
rm(temp1)

# Combine yearly values with yields
temp1 <- rapeseed_yield_melted[which(rapeseed_yield_melted$Farm != farm_ID),]
data_rapeseed_final_M4 <- join(temp1, final_rapeseed_temperature_aggregated_M4, by=c("Farm","year"))
rm(temperatures_ts_M4,temp1,final_rapeseed_temperature_aggregated_M4)

# New TS for temperature support of interest
marginal_temperatures_M4 <- rcspline.eval(seq(-10,39,0.225),knots=c(5,10,15,20,25), inclx=T)

temp_reg1 <- feols(yield ~ x + V2 + V3 + V4 + year + I(year^2) | Farm, data=data_rapeseed_final_M4)
summary(temp_reg1)
rm(temp_reg1)

# -------------------------------------------
# Statistical uncertainty
# -------------------------------------------

# Run 1000 models to derive uncertainty in estimates
spline_fitted_M4 <- matrix(NA,nrow=1000, ncol=nrow(marginal_temperatures_M4))

for (l in 1: 1000){
  # Sample some years
  data_temp <- data_rapeseed_final_M4[which(data_rapeseed_final_M4$year %in% sample(unique(data_rapeseed_final_M4$year),replace=T)),]
  
  # Run regression with random subsample (data_temp)
  fit_boot  <- lm( yield ~ x + V2+V3+V4 +year + I(year^2) + as.factor(Farm)-1, data=data_temp)
  
  # Save marginal effect of temperatures 
  spline_fitted_M4[l,]<-as.vector((marginal_temperatures_M4 %*% coef(fit_boot)[1:4]))
  
  # Show progress in %
  print(l/1000 *100)
  rm(fit_boot,data_temp)
}

# Change structure
spline_sub_melt_M4 <- melt(spline_fitted_M4)

# Get median response for each column (= a temperature from seq(3,39,0.225))
colmedians_spline_M4 <- colMedians(spline_fitted_M4)

# Get distance from median response for each value
dist_colmedian_M4 <- matrix(NA, nrow=nrow(spline_fitted_M4),ncol=ncol(spline_fitted_M4))

for(l in 1:ncol(spline_fitted_M4)){
  dist_colmedian_M4[,l] <- spline_fitted_M4[,l] - colmedians_spline_M4[l]
}

# Normalize distance to median
dist_colmedian_M4_normalized <- matrix(NA, nrow=nrow(spline_fitted_M4),ncol=ncol(spline_fitted_M4))

for(l in 1:ncol(spline_fitted_M4)){
  dist_colmedian_M4_normalized[,l]<-abs(dist_colmedian_M4[,l]/max(abs(dist_colmedian_M4[,l])))
}

# Change structure
dist_colmedian_M4_normalized_melted <- melt(dist_colmedian_M4_normalized)

# Last preparation before plotting

spline_merged_M4 <- merge(spline_sub_melt_M4,dist_colmedian_M4_normalized_melted,by=c("Var1","Var2"))
spline_merged_M4$value.y_minus1 <- 1-spline_merged_M4$value.y

temp1 <- aggregate(spline_merged_M4$value.y_minus1, list(spline_merged_M4$Var1), mean, na.rm=T) 
temp2 <- aggregate(spline_merged_M4$value.x, list(spline_merged_M4$Var2), mean, na.rm=T)

temp3 <- aggregate(spline_merged_M4$value.x, list(spline_merged_M4$Var2), quantile, probs=0.025) 
temp4 <- aggregate(spline_merged_M4$value.x, list(spline_merged_M4$Var2), quantile, probs=0.975)

spline_merged_M4 <- merge(spline_merged_M4, temp1, by.x=c("Var1"), by.y=c("Group.1"))
spline_merged_M4 <- merge(spline_merged_M4, temp2, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M4)[7]<-"average_effect"

spline_merged_M4<-merge(spline_merged_M4, temp3, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M4)[8]<-"lower_bound"

spline_merged_M4<-merge(spline_merged_M4, temp4, by.x=c("Var2"), by.y=c("Group.1"))
colnames(spline_merged_M4)[9]<-"upper_bound"

temp5 <- cbind(seq(-10,39,0.225),seq(1:length(seq(-10,39,0.225))))
colnames(temp5) <- c("Temp_bin", "Var2")
spline_merged_M4 <-merge(spline_merged_M4, temp5, by=c("Var2"))

# 
# Plot: hourly temperature effect
# 

model_4 <- ggplot()+ ggtitle("Hourly temperature effects: Model 4") + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "grey25", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "grey25"),
        axis.text = element_text(family="Times New Roman", size=10),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  xlab("Temperature (°C)") + ylab("Hourly yield response \n (dt/ha)")+
  scale_x_continuous(breaks = seq(-5,39, by=5), limits = c(-5,40))+
  scale_y_continuous(breaks = c(-0.2,-0.1,0,0.1), limits = c(-0.2,0.1))+
  geom_hline(yintercept=c(seq(-0.5,0.1,by=0.1)), colour="gray")+
  geom_vline(xintercept=knot_locations_kn5_rapeseed[,1], linetype="dotted", colour="gray")+
  geom_line(data=spline_merged_M4, aes(y=value.x, x=Temp_bin, group=Var1, alpha=x.x^200/200), size=.01, colour="brown1")+
  geom_line(data=spline_merged_M4, aes(y=average_effect, x=Temp_bin), alpha=0.85, size=.75, colour="brown4")+
  geom_line(data=spline_merged_M4, aes(y=lower_bound, x=Temp_bin)   , alpha=0.85, size=.75, colour="brown4")+       
  geom_line(data=spline_merged_M4, aes(y=upper_bound, x=Temp_bin)   , alpha=0.85, size=.75, colour="brown4")

rm(temp1,temp2,temp3,temp4,temp5)
rm(dist_colmedian_M4_normalized_melted, colmedians_spline_M4, spline_fitted_M4,dist_colmedian_M4_normalized,dist_colmedian_M4,spline_sub_melt_M4)

# 
# Preparation of hourly payouts
# 

# We change the sign of the effect (payouts are a positive number)

spline_merged_M4$average_effect <- spline_merged_M4$average_effect * (-1) 
spline_merged_M4$lower_bound <- spline_merged_M4$lower_bound * (-1) 
spline_merged_M4$upper_bound <- spline_merged_M4$upper_bound * (-1) 
spline_merged_M4$value.x <- spline_merged_M4$value.x * (-1) 

# Prepare lower bound (no values below 0)
temp <- subset(spline_merged_M4[,c(2,3,6,7,8,9,10)])
colnames(temp)[5] <- "upper_bound"
colnames(temp)[6] <- "lower_bound"

# Remove negative values(required for geom_ribbon())

temp[which(temp$lower_bound<0),6] <- 0

# 
# Plot: hourly payout function
# 

payout_rapeseed_model_4 <- ggplot(temp) +
  ggtitle("Hourly payouts: Model 4")+ theme(panel.grid.major = element_blank(),
                                                    plot.title = element_text(family="Times New Roman", size=12, colour = "grey25", hjust=0.5),
                                                    axis.title = element_text(family="Times New Roman", size=10, colour = "grey25"),
                                                    axis.text = element_text(family="Times New Roman", size=10),
                                                    panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(),
                                                    axis.line = element_line(colour = "black"),
                                                    legend.position="none")+
  xlab("Temperature (°C)") + ylab("Payout (dt/ha)")+
  scale_x_continuous(breaks = seq(10,39, by=5), limits = c(10,40))+
  scale_y_continuous(breaks = c(0,0.1, 0.2), limits = c(0,0.2))+
  geom_hline(yintercept=c(seq(0,0.4,by=0.1)), colour="gray")+
  geom_hline(yintercept=0, colour="darkgray")+
  geom_vline(xintercept= seq(5,38, by=5), colour="gray", linetype="dashed")+
  
  # Add lines
  
  geom_line(data=temp, aes(y=value.x, x=Temp_bin, group=Var1, alpha=x.x^200/200), size=.01, colour="royalblue1")+
  geom_line(data=temp, aes(y=average_effect, x=Temp_bin), size=0.75, color="royalblue4")+
  geom_line(data=temp, aes(y=lower_bound, x=Temp_bin), size=0.75, color="royalblue4")+
  geom_line(data=temp, aes(y=upper_bound, x=Temp_bin), size=0.75, color="royalblue4")

rm(temp,spline_merged_M4)

models_payouts_rapeseed <- ggarrange(model_1,payout_rapeseed_model_1,model_2,payout_rapeseed_model_2,model_3,payout_rapeseed_model_3,model_4,payout_rapeseed_model_4,temperature_rapeseed, ncol=2, byrow=T, heights = c(1,1,1,1,0.5))
#ggsave(plot=models_payouts_rapeseed, file="Plots/T Effect/models_payouts_rapeseed.png", width = 16, height = 20, unit="cm")
rm(model_1,payout_rapeseed_model_1,model_2,payout_rapeseed_model_2,model_3,payout_rapeseed_model_3,model_4,payout_rapeseed_model_4,temperature_rapeseed)

# ===========================================================================================================

# 3) Figure 4: Out-of-sample risk reductions

# ===========================================================================================================


# 15, 20 and 25 degrees have the following position in "strike_temperature"
reported_strikes <- c(15,20,25)
strike_temperature_levels <- which( strike_temperature %in% reported_strikes)

# -------------------------------------------
# For winter wheat
# -------------------------------------------

# Relative change in risk premium (for plots)
# Positive means risk-decreasing

# Calculation of changes in risk premium per farm [strike, farm, model]
rel_RP_wheat_models <- array(dim=c(3,nrow(farm_yields[[1]]), ncol=4))
dimnames(rel_RP_wheat_models) [[1]] <- c("15", "20", "25")
dimnames(rel_RP_wheat_models) [[3]] <- c("Model 1", "Model 2", "Model 3", "Model 4")

for (l in 1:length(strike_temperature_levels)){
  for (i in 1:nrow(farm_yields[[1]])){
    
    # Model 1: Best fit
    rel_RP_wheat_models[l,i,1] <- (-100)*(summary_EU_wheat_insured_kn3[[strike_temperature_levels[l]]] [1,i,4] - summary_EU_wheat_uninsured[i,4]) / summary_EU_wheat_uninsured[i,4]
    
    # Model 2: Equal spaces
    rel_RP_wheat_models[l,i,2] <- (-100)*(summary_EU_wheat_insured_kn3[[strike_temperature_levels[l]]] [2,i,4] - summary_EU_wheat_uninsured[i,4]) / summary_EU_wheat_uninsured[i,4]
    
    # Model 3: Quantiles
    rel_RP_wheat_models[l,i,3] <- (-100)*(summary_EU_wheat_insured_kn3[[strike_temperature_levels[l]]] [3,i,4] - summary_EU_wheat_uninsured[i,4]) / summary_EU_wheat_uninsured[i,4]
    
    # Model 4: Schlenker (5 space)
    rel_RP_wheat_models[l,i,4] <- (-100)*(summary_EU_wheat_insured_kn5[[strike_temperature_levels[l]]] [1,i,4] - summary_EU_wheat_uninsured[i,4]) / summary_EU_wheat_uninsured[i,4]
  }
}

# Plot for strike level 20
# rel_RP_wheat_models[2,,] contains RP changes for strike level temperature 20

# data.frame for Asteriks
label.df <- data.frame(Var2=c("Model 1","Model 2","Model 3","Model 4"), Value = c(45,45,45,45))

rel_RP_wheat_models_melted <- melt(as.matrix(rel_RP_wheat_models[2,,]))

wheat_20 <- ggplot(rel_RP_wheat_models_melted, aes(x=Var2, y=value)) + ggtitle("a) winter wheat \n strike level: 20°C") +
  geom_boxplot(notch = F, fill="gray", outlier.size = 0.1, coef=0) + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "black", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "black"),
        axis.text = element_text(family="Times New Roman", size=10, color="black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  ylab("Reduction in risk premium (%)") + xlab("")+
  scale_y_continuous(breaks = c(-20,-10,0,10,20,30,40), limits = c(-23, 50))+
  geom_hline(yintercept=c(-20,-10,0,10,20,30,40), linetype="dashed", colour="gray", alpha = 0.5)+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20", alpha=0.6)+
  geom_text(data=label.df, aes(x=Var2, y=c(45,45,45,45)), label="***", size=4)


# Plot for strike level 25
rel_RP_wheat_models_melted_25 <- melt(as.matrix(rel_RP_wheat_models[3,,]))

wheat_25 <- ggplot(rel_RP_wheat_models_melted_25, aes(x=Var2, y=value)) + ggtitle("b) winter wheat \n strike level: 25°C") +
  geom_boxplot(notch = F, fill="gray", outlier.size = 0.1, coef=0) + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "black", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "black"),
        axis.text = element_text(family="Times New Roman", size=10, color="black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  ylab("Reduction in risk premium (%)") + xlab("")+
  scale_y_continuous(breaks = c(-20,-10,0,10,20,30,40), limits = c(-23, 50))+
  geom_hline(yintercept=c(-20,-10,0,10,20,30,40), linetype="dashed", colour="gray", alpha=0.5)+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20", alpha=0.6)+
  geom_text(data=label.df, aes(x=Var2, y=c(45,45,45,45)), label="***", size=4)

# -------------------------------------------
# For rapeseed
# -------------------------------------------

# Relative change in risk premium (for plots)
# Positive means risk-decreasing

rel_RP_rapeseed_models <- array(dim=c(3,nrow(farm_yields[[2]]), ncol=4))
dimnames(rel_RP_rapeseed_models) [[1]] <- c("15", "20", "25")
dimnames(rel_RP_rapeseed_models) [[3]] <- c("Model 1", "Model 2", "Model 3", "Model 4")

for (l in 1:length(strike_temperature_levels)){
  for (i in 1:nrow(farm_yields[[2]])){
    
    # Model 1: Best fit
    rel_RP_rapeseed_models[l,i,1] <- (-100)*(summary_EU_rapeseed_insured_kn3[[strike_temperature_levels[l]]] [1,i,4] - summary_EU_rapeseed_uninsured[i,4]) / summary_EU_rapeseed_uninsured[i,4]
    # Model 2: Equal spaces
    rel_RP_rapeseed_models[l,i,2] <- (-100)*(summary_EU_rapeseed_insured_kn3[[strike_temperature_levels[l]]] [2,i,4] - summary_EU_rapeseed_uninsured[i,4]) / summary_EU_rapeseed_uninsured[i,4]
    # Model 3: Quantiles
    rel_RP_rapeseed_models[l,i,3] <- (-100)*(summary_EU_rapeseed_insured_kn3[[strike_temperature_levels[l]]] [3,i,4] - summary_EU_rapeseed_uninsured[i,4]) / summary_EU_rapeseed_uninsured[i,4]
    # Model 4: (5 spaces)
    rel_RP_rapeseed_models[l,i,4] <- (-100)*(summary_EU_rapeseed_insured_kn5[[strike_temperature_levels[l]]] [1,i,4] - summary_EU_rapeseed_uninsured[i,4]) / summary_EU_rapeseed_uninsured[i,4]
  }
}

# Plot for strike level 15
rel_RP_rapeseed_models_melted <- melt(as.matrix(rel_RP_rapeseed_models[1,,]))

rapeseed_15 <- ggplot(rel_RP_rapeseed_models_melted, aes(x=Var2, y=value)) + ggtitle("c) winter rapeseed \n strike level: 15°C") +
  geom_boxplot(notch = F, fill="gray", outlier.size = 0.1, coef=0) + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "black", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "black"),
        axis.text = element_text(family="Times New Roman", size=10, color="black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  ylab("Reduction in risk premium (%)") + xlab("")+
  scale_y_continuous(breaks = c(-40,-20,0,20,40,60),limits = c(-45, 70))+ 
  geom_hline(yintercept=c(-40,-20,0,20,40,60), linetype="dashed", colour="gray", alpha=0.5)+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20", alpha=0.6)+
  geom_text(data=label.df, aes(x=Var2, y=c(65,65,65,65)), label="***", size=4)

# Plot for strike level 20
rel_RP_rapeseed_models_melted_20 <- melt(as.matrix(rel_RP_rapeseed_models[2,,]))

rapeseed_20 <- ggplot(rel_RP_rapeseed_models_melted_20, aes(x=Var2, y=value)) + ggtitle("d) winter rapeseed \n strike level: 20°C") +
  geom_boxplot(notch = F, fill="gray", outlier.size = 0.1, coef=0) + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "black", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "black"),
        axis.text = element_text(family="Times New Roman", size=10, color="black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  ylab("Reduction in risk premium (%)") + xlab("")+ 
  scale_y_continuous(breaks = c(-40,-20,0,20,40,60), limits = c(-45, 70))+  
  geom_hline(yintercept=c(-40,-20,0,20,40,60), linetype="dashed", colour="gray", alpha=0.6)+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20", alpha=0.6)+
  geom_text(data=label.df, aes(x=Var2, y=c(65,65,65,65)), label="***", size=4)

plot_risk <- ggarrange(wheat_20, wheat_25,rapeseed_15,rapeseed_20, nrow=2, ncol=2)
#ggsave(plot=plot_risk, file="Plots/risk-reduction/risk-reductions_round2_extremeaversion.png", width = 16, height = 13, unit="cm")

