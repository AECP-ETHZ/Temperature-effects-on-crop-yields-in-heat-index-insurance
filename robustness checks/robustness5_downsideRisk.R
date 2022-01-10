# =========================================================================================================
# ---------------------------------------------------------------------------------------------------------
#
# Temperature Effects on Crop Yields in Heat Index Insurance
#
# Supplementary R code for robustness check with lower partial moments
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

# Before you run this robustness check, run main_codes_masterfile.R

# =========================================================================================================
# 1. Loading packages 
# =========================================================================================================

library(PerformanceAnalytics)
library(reshape2)
library(ggplot2)

# =========================================================================================================
# 2. Lower partial moments (from a revennue/yield quantile)
# =========================================================================================================

# The paper shows results for the 50% quantile (=median)

quantile_interest <- 0.5

# Get farm-specific wheat yield quantiles (threshold yield)
quantile_wheat <- vector()
for(i in 1:nrow(detr_wheat_yields)){
  quantile_wheat[i] <- quantile(detr_wheat_yields[i,], quantile_interest, type=1, na.rm=T)
}

# ---------------------
# Wheat: uninsured
# ---------------------

LPM1_uninsured <- vector()
LPM2_uninsured <- vector()

for (i in 1:nrow(detr_wheat_yields)){
  
  temp1 <- detr_wheat_yields[i,!is.na(detr_wheat_yields[i,])]
  
  # LPM1
  shortfall_uninsured <- quantile_wheat[i] - temp1
  shortfall_uninsured <- shortfall_uninsured[shortfall_uninsured>0]
  LPM1_uninsured[i] <- mean(shortfall_uninsured)
  
  #LPM2
  variance_uninsured  <- (shortfall_uninsured)^2
  LPM2_uninsured [i]          <- mean(variance_uninsured)
  
  rm(temp1,shortfall_uninsured,variance_uninsured)
}

# ---------------------------------------
# Wheat: Model 1 with strike level 20 °C
# ---------------------------------------

LPM1_M1_20 <- vector()
LPM2_M1_20 <- vector()

for (i in 1:nrow(detr_wheat_yields)){
  # Yields (=wealth) for 20 °C strike level (element number 8 in list; 1st dimension in array = model 1)
  temp1 <- list_farm_wheat_wealth_modelskn3_strike[[8]][1,i,!is.na(detr_wheat_yields[i,])]
  
  # LPM1
  shortfall_uninsured <- quantile_wheat[i] - temp1
  shortfall_uninsured <- shortfall_uninsured[shortfall_uninsured>0]
  LPM1_M1_20[i] <- mean(shortfall_uninsured)
  
  #LPM2
  variance_uninsured  <- (shortfall_uninsured)^2
  LPM2_M1_20[i]          <- mean(variance_uninsured)
  
  rm(temp1,shortfall_uninsured,variance_uninsured)
}

# ----------------------------------------
# Wheat: Model 2 with strike level 20 °C
# ----------------------------------------

LPM1_M2_20 <- vector()
LPM2_M2_20 <- vector()

for (i in 1:nrow(detr_wheat_yields)){
  temp1 <- list_farm_wheat_wealth_modelskn3_strike[[8]][2,i,!is.na(detr_wheat_yields[i,])]
  
  # LPM1
  shortfall_uninsured <- quantile_wheat[i] - temp1
  shortfall_uninsured <- shortfall_uninsured[shortfall_uninsured>0]
  LPM1_M2_20[i] <- mean(shortfall_uninsured)
  
  #LPM2
  variance_uninsured  <- (shortfall_uninsured)^2
  LPM2_M2_20[i]          <- mean(variance_uninsured)
  
  rm(temp1,shortfall_uninsured,variance_uninsured)
}

# ----------------------------------------
# Wheat: Model 3 with strike level 20 °C
# ----------------------------------------

LPM1_M3_20 <- vector()
LPM2_M3_20 <- vector()

for (i in 1:nrow(detr_wheat_yields)){
  temp1 <- list_farm_wheat_wealth_modelskn3_strike[[8]][3,i,!is.na(detr_wheat_yields[i,])]
  
  # LPM1
  shortfall_uninsured <- quantile_wheat[i] - temp1
  shortfall_uninsured <- shortfall_uninsured[shortfall_uninsured>0]
  LPM1_M3_20[i] <- mean(shortfall_uninsured)
  
  #LPM2
  variance_uninsured  <- (shortfall_uninsured)^2
  LPM2_M3_20[i]          <- mean(variance_uninsured)
  
  rm(temp1,shortfall_uninsured,variance_uninsured)
}

# ----------------------------------------
# Wheat: Model 4 with strike level 20 °C
# ----------------------------------------

LPM1_M4_20 <- vector()
LPM2_M4_20 <- vector()

for (i in 1:nrow(detr_wheat_yields)){
  temp1 <- list_farm_wheat_wealth_modelskn5_strike[[8]][1,i,!is.na(detr_wheat_yields[i,])]
  
  # LPM1
  shortfall_uninsured <- quantile_wheat[i] - temp1
  shortfall_uninsured <- shortfall_uninsured[shortfall_uninsured>0]
  LPM1_M4_20[i] <- mean(shortfall_uninsured)
  
  #LPM2
  variance_uninsured  <- (shortfall_uninsured)^2
  LPM2_M4_20[i]          <- mean(variance_uninsured)
  
  rm(temp1,shortfall_uninsured,variance_uninsured)
}

LPM1_wheat_20 <- cbind(LPM1_uninsured, LPM1_M1_20, LPM1_M2_20, LPM1_M3_20, LPM1_M4_20)
colnames(LPM1_wheat_20) <- c("uninsured","M1","M2","M3","M4")

LPM2_wheat_20 <- cbind(LPM2_uninsured, LPM2_M1_20, LPM2_M2_20, LPM2_M3_20, LPM2_M4_20)
colnames(LPM1_wheat_20) <- c("uninsured","M1","M2","M3","M4")
rm(LPM1_M1_20, LPM1_M2_20, LPM1_M3_20, LPM1_M4_20, LPM2_M1_20, LPM2_M2_20, LPM2_M3_20, LPM2_M4_20)

# ----------------------------------------
# Wheat: Model 1 with strike level 25 °C
# ----------------------------------------

LPM1_M1_25 <- vector()
LPM2_M1_25 <- vector()

for (i in 1:nrow(detr_wheat_yields)){
  temp1 <- list_farm_wheat_wealth_modelskn3_strike[[13]][1,i,!is.na(detr_wheat_yields[i,])]
  
  # LPM1
  shortfall_uninsured <- quantile_wheat[i] - temp1
  shortfall_uninsured <- shortfall_uninsured[shortfall_uninsured>0]
  LPM1_M1_25[i] <- mean(shortfall_uninsured)
  
  #LPM2
  variance_uninsured  <- (shortfall_uninsured)^2
  LPM2_M1_25[i]          <- mean(variance_uninsured)
  
  rm(temp1,shortfall_uninsured,variance_uninsured)
}

# ----------------------------------------
# Wheat: Model 2 with strike level 25 °C
# ----------------------------------------

LPM1_M2_25 <- vector()
LPM2_M2_25 <- vector()

for (i in 1:nrow(detr_wheat_yields)){
  temp1 <- list_farm_wheat_wealth_modelskn3_strike[[13]][2,i,!is.na(detr_wheat_yields[i,])]
  
  # LPM1
  shortfall_uninsured <- quantile_wheat[i] - temp1
  shortfall_uninsured <- shortfall_uninsured[shortfall_uninsured>0]
  LPM1_M2_25[i] <- mean(shortfall_uninsured)
  
  #LPM2
  variance_uninsured  <- (shortfall_uninsured)^2
  LPM2_M2_25[i]          <- mean(variance_uninsured)
  
  rm(temp1,shortfall_uninsured,variance_uninsured)
}

# ----------------------------------------
# Wheat: Model 3 with strike level 25 °C
# ----------------------------------------

LPM1_M3_25 <- vector()
LPM2_M3_25 <- vector()

for (i in 1:nrow(detr_wheat_yields)){
  temp1 <- list_farm_wheat_wealth_modelskn3_strike[[13]][3,i,!is.na(detr_wheat_yields[i,])]
  
  # LPM1
  shortfall_uninsured <- quantile_wheat[i] - temp1
  shortfall_uninsured <- shortfall_uninsured[shortfall_uninsured>0]
  LPM1_M3_25[i] <- mean(shortfall_uninsured)
  
  #LPM2
  variance_uninsured  <- (shortfall_uninsured)^2
  LPM2_M3_25[i]          <- mean(variance_uninsured)
  
  rm(temp1,shortfall_uninsured,variance_uninsured)
}

# ----------------------------------------
# Wheat: Model 4 with strike level 25 °C
# ----------------------------------------

LPM1_M4_25 <- vector()
LPM2_M4_25 <- vector()

for (i in 1:nrow(detr_wheat_yields)){
  temp1 <- list_farm_wheat_wealth_modelskn5_strike[[13]][1,i,!is.na(detr_wheat_yields[i,])]
  
  # LPM1
  shortfall_uninsured <- quantile_wheat[i] - temp1
  shortfall_uninsured <- shortfall_uninsured[shortfall_uninsured>0]
  LPM1_M4_25[i] <- mean(shortfall_uninsured)
  
  #LPM2
  variance_uninsured  <- (shortfall_uninsured)^2
  LPM2_M4_25[i]          <- mean(variance_uninsured)
  
  rm(temp1,shortfall_uninsured,variance_uninsured)
}

LPM1_wheat_25 <- cbind(LPM1_uninsured, LPM1_M1_25, LPM1_M2_25, LPM1_M3_25, LPM1_M4_25)
colnames(LPM1_wheat_25) <- c("uninsured","M1","M2","M3","M4")

LPM2_wheat_25 <- cbind(LPM2_uninsured, LPM2_M1_25, LPM2_M2_25, LPM2_M3_25, LPM2_M4_25)
colnames(LPM2_wheat_25) <- c("uninsured","M1","M2","M3","M4")
rm(LPM1_M1_25, LPM1_M2_25, LPM1_M3_25, LPM1_M4_25,  LPM2_M1_25, LPM2_M2_25, LPM2_M3_25, LPM2_M4_25)


# ---------------------
# rapeseed: uninsured
# ---------------------

# Get farm-specific rapeseed yield quantiles
quantile_rapeseed <- vector()

for(i in 1:nrow(detr_rapeseed_yields)){
  quantile_rapeseed[i] <- quantile(detr_rapeseed_yields[i,], quantile_interest, type=1, na.rm=T)
}

LPM1_uninsured <- vector()
LPM2_uninsured <- vector()

for (i in 1:nrow(detr_rapeseed_yields)){
  
  temp1 <- detr_rapeseed_yields[i,!is.na(detr_rapeseed_yields[i,])]
  
  # LPM1
  shortfall_uninsured <- quantile_rapeseed[i] - temp1
  shortfall_uninsured <- shortfall_uninsured[shortfall_uninsured>0]
  LPM1_uninsured[i] <- mean(shortfall_uninsured)
  
  #LPM2
  variance_uninsured  <- (shortfall_uninsured)^2
  LPM2_uninsured [i]          <- mean(variance_uninsured)
  
  rm(temp1,shortfall_uninsured,variance_uninsured)
}

# -----------------------------------------
# Rapeseed: Model 1 with strike level 20 °C
# -----------------------------------------

LPM1_M1_20 <- vector()
LPM2_M1_20 <- vector()

for (i in 1:nrow(detr_rapeseed_yields)){
  temp1 <- list_farm_rapeseed_wealth_modelskn3_strike[[8]][1,i,!is.na(detr_rapeseed_yields[i,])]
  
  # LPM1
  shortfall_uninsured <- quantile_rapeseed[i] - temp1
  shortfall_uninsured <- shortfall_uninsured[shortfall_uninsured>0]
  LPM1_M1_20[i] <- mean(shortfall_uninsured)
  
  #LPM2
  variance_uninsured  <- (shortfall_uninsured)^2
  LPM2_M1_20[i]          <- mean(variance_uninsured)
  
  rm(temp1,shortfall_uninsured,variance_uninsured)
}

# -----------------------------------------
# Rapeseed: Model 2 with strike level 20 °C
# -----------------------------------------

LPM1_M2_20 <- vector()
LPM2_M2_20 <- vector()

for (i in 1:nrow(detr_rapeseed_yields)){
  temp1 <- list_farm_rapeseed_wealth_modelskn3_strike[[8]][2,i,!is.na(detr_rapeseed_yields[i,])]
  
  # LPM1
  shortfall_uninsured <- quantile_rapeseed[i] - temp1
  shortfall_uninsured <- shortfall_uninsured[shortfall_uninsured>0]
  LPM1_M2_20[i] <- mean(shortfall_uninsured)
  
  #LPM2
  variance_uninsured  <- (shortfall_uninsured)^2
  LPM2_M2_20[i]          <- mean(variance_uninsured)
  
  rm(temp1,shortfall_uninsured,variance_uninsured)
}

# -----------------------------------------
# Rapeseed: Model 3 with strike level 20 °C
# -----------------------------------------

LPM1_M3_20 <- vector()
LPM2_M3_20 <- vector()

for (i in 1:nrow(detr_rapeseed_yields)){
  temp1 <- list_farm_rapeseed_wealth_modelskn3_strike[[8]][3,i,!is.na(detr_rapeseed_yields[i,])]
  
  # LPM1
  shortfall_uninsured <- quantile_rapeseed[i] - temp1
  shortfall_uninsured <- shortfall_uninsured[shortfall_uninsured>0]
  LPM1_M3_20[i] <- mean(shortfall_uninsured)
  
  #LPM2
  variance_uninsured  <- (shortfall_uninsured)^2
  LPM2_M3_20[i]          <- mean(variance_uninsured)
  
  rm(temp1,shortfall_uninsured,variance_uninsured)
}

# -----------------------------------------
# Rapeseed: Model 4 with strike level 20 °C
# -----------------------------------------

LPM1_M4_20 <- vector()
LPM2_M4_20 <- vector()

for (i in 1:nrow(detr_rapeseed_yields)){
  temp1 <- list_farm_rapeseed_wealth_modelskn5_strike[[8]][1,i,!is.na(detr_rapeseed_yields[i,])]
  
  # LPM1
  shortfall_uninsured <- quantile_rapeseed[i] - temp1
  shortfall_uninsured <- shortfall_uninsured[shortfall_uninsured>0]
  LPM1_M4_20[i] <- mean(shortfall_uninsured)
  
  #LPM2
  variance_uninsured  <- (shortfall_uninsured)^2
  LPM2_M4_20[i]          <- mean(variance_uninsured)
  
  rm(temp1,shortfall_uninsured,variance_uninsured)
}

LPM1_rapeseed_20 <- cbind(LPM1_uninsured, LPM1_M1_20, LPM1_M2_20, LPM1_M3_20, LPM1_M4_20)
colnames(LPM1_rapeseed_20) <- c("uninsured","M1","M2","M3","M4")

LPM2_rapeseed_20 <- cbind(LPM2_uninsured, LPM2_M1_20, LPM2_M2_20, LPM2_M3_20, LPM2_M4_20)
colnames(LPM1_rapeseed_20) <- c("uninsured","M1","M2","M3","M4")
rm( LPM1_M1_20, LPM1_M2_20, LPM1_M3_20, LPM1_M4_20,LPM2_M1_20, LPM2_M2_20, LPM2_M3_20, LPM2_M4_20)

# -----------------------------------------
# Rapeseed: Model 1 with strike level 15 °C
# -----------------------------------------

LPM1_M1_15 <- vector()
LPM2_M1_15 <- vector()

for (i in 1:nrow(detr_rapeseed_yields)){
  temp1 <- list_farm_rapeseed_wealth_modelskn3_strike[[3]][1,i,!is.na(detr_rapeseed_yields[i,])]
  
  # LPM1
  shortfall_uninsured <- quantile_rapeseed[i] - temp1
  shortfall_uninsured <- shortfall_uninsured[shortfall_uninsured>0]
  LPM1_M1_15[i] <- mean(shortfall_uninsured)
  
  #LPM2
  variance_uninsured  <- (shortfall_uninsured)^2
  LPM2_M1_15[i]          <- mean(variance_uninsured)
  
  rm(temp1,shortfall_uninsured,variance_uninsured)
}

# -----------------------------------------
# Rapeseed: Model 2 with strike level 15 °C
# -----------------------------------------

LPM1_M2_15 <- vector()
LPM2_M2_15 <- vector()

for (i in 1:nrow(detr_rapeseed_yields)){
  temp1 <- list_farm_rapeseed_wealth_modelskn3_strike[[3]][2,i,!is.na(detr_rapeseed_yields[i,])]
  
  # LPM1
  shortfall_uninsured <- quantile_rapeseed[i] - temp1
  shortfall_uninsured <- shortfall_uninsured[shortfall_uninsured>0]
  LPM1_M2_15[i] <- mean(shortfall_uninsured)
  
  #LPM2
  variance_uninsured  <- (shortfall_uninsured)^2
  LPM2_M2_15[i]          <- mean(variance_uninsured)
  
  rm(temp1,shortfall_uninsured,variance_uninsured)
}

# -----------------------------------------
# Rapeseed: Model 3 with strike level 15 °C
# -----------------------------------------

LPM1_M3_15 <- vector()
LPM2_M3_15 <- vector()

for (i in 1:nrow(detr_rapeseed_yields)){
  temp1 <- list_farm_rapeseed_wealth_modelskn3_strike[[3]][3,i,!is.na(detr_rapeseed_yields[i,])]
  
  # LPM1
  shortfall_uninsured <- quantile_rapeseed[i] - temp1
  shortfall_uninsured <- shortfall_uninsured[shortfall_uninsured>0]
  LPM1_M3_15[i] <- mean(shortfall_uninsured)
  
  #LPM2
  variance_uninsured  <- (shortfall_uninsured)^2
  LPM2_M3_15[i]          <- mean(variance_uninsured)
  
  rm(temp1,shortfall_uninsured,variance_uninsured)
}

# -----------------------------------------
# Rapeseed: Model 4 with strike level 15 °C
# -----------------------------------------

LPM1_M4_15 <- vector()
LPM2_M4_15 <- vector()

for (i in 1:nrow(detr_rapeseed_yields)){
  temp1 <- list_farm_rapeseed_wealth_modelskn5_strike[[3]][1,i,!is.na(detr_rapeseed_yields[i,])]
  
  # LPM1
  shortfall_uninsured <- quantile_rapeseed[i] - temp1
  shortfall_uninsured <- shortfall_uninsured[shortfall_uninsured>0]
  LPM1_M4_15[i] <- mean(shortfall_uninsured)
  
  #LPM2
  variance_uninsured  <- (shortfall_uninsured)^2
  LPM2_M4_15[i]          <- mean(variance_uninsured)
  
  rm(temp1,shortfall_uninsured,variance_uninsured)
}

LPM1_rapeseed_15 <- cbind(LPM1_uninsured, LPM1_M1_15, LPM1_M2_15, LPM1_M3_15, LPM1_M4_15)
colnames(LPM1_rapeseed_15) <- c("uninsured","M1","M2","M3","M4")

LPM2_rapeseed_15 <- cbind(LPM2_uninsured, LPM2_M1_15, LPM2_M2_15, LPM2_M3_15, LPM2_M4_15)
colnames(LPM2_rapeseed_15) <- c("uninsured","M1","M2","M3","M4")


# ---------------------
# Plots LPM1 
# ---------------------

# 
# Plot LPM1 wheat 20
# 

label.df <- data.frame(Var2=c("Model 1","Model 2","Model 3","Model 4"), Value = c(45,45,45,45))

# Check levels of significance manually: [,1] is uninsured. Then adjust significance levels manually in plots
wilcox.test(LPM1_wheat_20[,5], LPM1_wheat_20[,1], paired=T, alternative="l")$p.value

rel_red_LPM1_wheat_20 <- as.data.frame(matrix(NA, nrow=nrow(detr_wheat_yields), ncol=4))
colnames(rel_red_LPM1_wheat_20) <- c("Model 1","Model 2", "Model 3", "Model 4")
rel_red_LPM1_wheat_20[,1] <- ((LPM1_wheat_20[,2] - LPM1_wheat_20[,1]) / LPM1_wheat_20[,1]) * (-100)
rel_red_LPM1_wheat_20[,2] <- ((LPM1_wheat_20[,3] - LPM1_wheat_20[,1]) / LPM1_wheat_20[,1]) * (-100)
rel_red_LPM1_wheat_20[,3] <- ((LPM1_wheat_20[,4] - LPM1_wheat_20[,1]) / LPM1_wheat_20[,1]) * (-100)
rel_red_LPM1_wheat_20[,4] <- ((LPM1_wheat_20[,5] - LPM1_wheat_20[,1]) / LPM1_wheat_20[,1]) * (-100)
rel_red_LPM1_wheat_20_melted <- melt(rel_red_LPM1_wheat_20)

wheat_20 <- ggplot(rel_red_LPM1_wheat_20_melted, aes(x=variable, y=value)) + ggtitle("a) winter wheat \n strike level: 20°C") +
  geom_boxplot(notch = F, fill="gray", outlier.size = 0.1, coef=0) + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "black", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "black"),
        axis.text = element_text(family="Times New Roman", size=10, color="black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  ylab("Reduction in expected shortfall (%)") + xlab("")+
  scale_y_continuous(breaks = c(-40,-20,0,20,40), limits = c(-40, 55))+
  geom_hline(yintercept=c(-40,-20,0,20,40), linetype="dashed", colour="gray", alpha = 0.5)+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20", alpha=0.6)+
  geom_text(data=label.df, aes(x=Var2, y=c(55,55,55,55)), label="***", size=4)

# 
# Plot LPM1 wheat 25
# 

label.df <- data.frame(Var2=c("Model 1","Model 2","Model 3","Model 4"), Value = c(45,45,45,45))
wilcox.test(LPM1_wheat_25[,1], LPM1_wheat_25[,1], paired=T, alternative="l")$p.value

rel_red_LPM1_wheat_25 <- as.data.frame(matrix(NA, nrow=nrow(detr_wheat_yields), ncol=4))
colnames(rel_red_LPM1_wheat_25) <- c("Model 1","Model 2", "Model 3", "Model 4")
rel_red_LPM1_wheat_25[,1] <- ((LPM1_wheat_25[,2] - LPM1_wheat_25[,1]) / LPM1_wheat_25[,1]) * (-100)
rel_red_LPM1_wheat_25[,2] <- ((LPM1_wheat_25[,3] - LPM1_wheat_25[,1]) / LPM1_wheat_25[,1]) * (-100)
rel_red_LPM1_wheat_25[,3] <- ((LPM1_wheat_25[,4] - LPM1_wheat_25[,1]) / LPM1_wheat_25[,1]) * (-100)
rel_red_LPM1_wheat_25[,4] <- ((LPM1_wheat_25[,5] - LPM1_wheat_25[,1]) / LPM1_wheat_25[,1]) * (-100)
rel_red_LPM1_wheat_25_melted <- melt(rel_red_LPM1_wheat_25)

wheat_25 <- ggplot(rel_red_LPM1_wheat_25_melted, aes(x=variable, y=value)) + ggtitle("b) winter wheat \n strike level: 25°C") +
  geom_boxplot(notch = F, fill="gray", outlier.size = 0.1, coef=0) + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "black", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "black"),
        axis.text = element_text(family="Times New Roman", size=10, color="black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  ylab("Reduction in expected shortfall (%)") + xlab("")+
  scale_y_continuous(breaks = c(-40,-20,0,20,40), limits = c(-40, 55))+
  geom_hline(yintercept=c(-40,-20,0,20,40), linetype="dashed", colour="gray", alpha = 0.5)+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20", alpha=0.6)+
  geom_text(data=label.df, aes(x=Var2, y=c(55,55,55,55)), label="***", size=4)

# 
# Plot LPM1 rapeseed 15
# 

label.df <- data.frame(Var2=c("Model 1","Model 2","Model 3","Model 4"), Value = c(45,45,45,45))
wilcox.test(LPM1_rapeseed_15[,5], LPM1_rapeseed_15[,1], paired=T, alternative="l")$p.value

rel_red_LPM1_rapeseed_15 <- as.data.frame(matrix(NA, nrow=nrow(detr_rapeseed_yields), ncol=4))
colnames(rel_red_LPM1_rapeseed_15) <- c("Model 1","Model 2", "Model 3", "Model 4")
rel_red_LPM1_rapeseed_15[,1] <- ((LPM1_rapeseed_15[,2] - LPM1_rapeseed_15[,1]) / LPM1_rapeseed_15[,1]) * (-100)
rel_red_LPM1_rapeseed_15[,2] <- ((LPM1_rapeseed_15[,3] - LPM1_rapeseed_15[,1]) / LPM1_rapeseed_15[,1]) * (-100)
rel_red_LPM1_rapeseed_15[,3] <- ((LPM1_rapeseed_15[,4] - LPM1_rapeseed_15[,1]) / LPM1_rapeseed_15[,1]) * (-100)
rel_red_LPM1_rapeseed_15[,4] <- ((LPM1_rapeseed_15[,5] - LPM1_rapeseed_15[,1]) / LPM1_rapeseed_15[,1]) * (-100)
rel_red_LPM1_rapeseed_15_melted <- melt(rel_red_LPM1_rapeseed_15)

rapeseed_15 <- ggplot(rel_red_LPM1_rapeseed_15_melted, aes(x=variable, y=value)) + ggtitle("c) winter rapeseed \n strike level: 15°C") +
  geom_boxplot(notch = F, fill="gray", outlier.size = 0.1, coef=0) + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "black", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "black"),
        axis.text = element_text(family="Times New Roman", size=10, color="black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  ylab("Reduction in expected shortfall (%)") + xlab("")+
  scale_y_continuous(breaks = c(-40,-20,0,20,40), limits = c(-40, 55))+
  geom_hline(yintercept=c(-40,-20,0,20,40), linetype="dashed", colour="gray", alpha = 0.5)+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20", alpha=0.6)+
  geom_text(data=label.df, aes(x=Var2, y=c(55,55,55,55)), label="***", size=4)

# 
# Plot LPM1 rapeseed 20
# 

label.df <- data.frame(Var2=c("Model 1","Model 2","Model 3","Model 4"), Value = c(45,45,45,45))
wilcox.test(LPM1_rapeseed_20[,3], LPM1_rapeseed_20[,1], paired=T, alternative="l")$p.value

rel_red_LPM1_rapeseed_20 <- as.data.frame(matrix(NA, nrow=nrow(detr_rapeseed_yields), ncol=4))
colnames(rel_red_LPM1_rapeseed_20) <- c("Model 1","Model 2", "Model 3", "Model 4")
rel_red_LPM1_rapeseed_20[,1] <- ((LPM1_rapeseed_20[,2] - LPM1_rapeseed_20[,1]) / LPM1_rapeseed_20[,1]) * (-100)
rel_red_LPM1_rapeseed_20[,2] <- ((LPM1_rapeseed_20[,3] - LPM1_rapeseed_20[,1]) / LPM1_rapeseed_20[,1]) * (-100)
rel_red_LPM1_rapeseed_20[,3] <- ((LPM1_rapeseed_20[,4] - LPM1_rapeseed_20[,1]) / LPM1_rapeseed_20[,1]) * (-100)
rel_red_LPM1_rapeseed_20[,4] <- ((LPM1_rapeseed_20[,5] - LPM1_rapeseed_20[,1]) / LPM1_rapeseed_20[,1]) * (-100)
rel_red_LPM1_rapeseed_20_melted <- melt(rel_red_LPM1_rapeseed_20)

rapeseed_20 <- ggplot(rel_red_LPM1_rapeseed_20_melted, aes(x=variable, y=value)) + ggtitle("d) winter rapeseed \n strike level: 20°C") +
  geom_boxplot(notch = F, fill="gray", outlier.size = 0.1, coef=0) + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "black", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "black"),
        axis.text = element_text(family="Times New Roman", size=10, color="black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  ylab("Reduction in expected shortfall (%)") + xlab("")+
  scale_y_continuous(breaks = c(-40,-20,0,20,40), limits = c(-40, 55))+
  geom_hline(yintercept=c(-40,-20,0,20,40), linetype="dashed", colour="gray", alpha = 0.5)+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20", alpha=0.6)+
  geom_text(data=label.df, aes(x=Var2, y=c(55,55,55,55)), label="***", size=4)

LPM1 <- ggarrange(wheat_20, wheat_25,rapeseed_15,rapeseed_20, nrow=2, ncol=2)
#ggsave(plot=LPM1, file="Plots/risk-reduction/LPM1_APP.png", width = 16, height = 13, unit="cm")


# ---------------------
# Plots LPM2 
# ---------------------

# 
# Plot LPM2 wheat 20
# 

label.df <- data.frame(Var2=c("Model 1","Model 2","Model 3","Model 4"), Value = c(45,45,45,45))
wilcox.test(LPM2_wheat_20[,2], LPM2_wheat_20[,1], paired=T, alternative="l")$p.value

rel_red_LPM2_wheat_20 <- as.data.frame(matrix(NA, nrow=nrow(detr_wheat_yields), ncol=4))
colnames(rel_red_LPM2_wheat_20) <- c("Model 1","Model 2", "Model 3", "Model 4")
rel_red_LPM2_wheat_20[,1] <- ((LPM2_wheat_20[,2] - LPM2_wheat_20[,1]) / LPM2_wheat_20[,1]) * (-100)
rel_red_LPM2_wheat_20[,2] <- ((LPM2_wheat_20[,3] - LPM2_wheat_20[,1]) / LPM2_wheat_20[,1]) * (-100)
rel_red_LPM2_wheat_20[,3] <- ((LPM2_wheat_20[,4] - LPM2_wheat_20[,1]) / LPM2_wheat_20[,1]) * (-100)
rel_red_LPM2_wheat_20[,4] <- ((LPM2_wheat_20[,5] - LPM2_wheat_20[,1]) / LPM2_wheat_20[,1]) * (-100)
rel_red_LPM2_wheat_20_melted <- melt(rel_red_LPM2_wheat_20)

wheat_20 <- ggplot(rel_red_LPM2_wheat_20_melted, aes(x=variable, y=value)) + ggtitle("a) winter wheat \n strike level: 20°C") +
  geom_boxplot(notch = F, fill="gray", outlier.size = 0.1, coef=0) + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "black", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "black"),
        axis.text = element_text(family="Times New Roman", size=10, color="black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  ylab("Reduction in downside variance (%)") + xlab("")+
  scale_y_continuous(breaks = c(-40,-20,0,20,40), limits = c(-40, 55))+
  geom_hline(yintercept=c(-40,-20,0,20,40), linetype="dashed", colour="gray", alpha = 0.5)+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20", alpha=0.6)+
  geom_text(data=label.df, aes(x=Var2, y=c(55,55,55,55)), label="***", size=4)

# 
# Plot LPM2 wheat 25
# 

label.df <- data.frame(Var2=c("Model 1","Model 2","Model 3","Model 4"), Value = c(45,45,45,45))
wilcox.test(LPM2_wheat_25[,1], LPM2_wheat_25[,1], paired=T, alternative="l")$p.value

rel_red_LPM2_wheat_25 <- as.data.frame(matrix(NA, nrow=nrow(detr_wheat_yields), ncol=4))
colnames(rel_red_LPM2_wheat_25) <- c("Model 1","Model 2", "Model 3", "Model 4")
rel_red_LPM2_wheat_25[,1] <- ((LPM2_wheat_25[,2] - LPM2_wheat_25[,1]) / LPM2_wheat_25[,1]) * (-100)
rel_red_LPM2_wheat_25[,2] <- ((LPM2_wheat_25[,3] - LPM2_wheat_25[,1]) / LPM2_wheat_25[,1]) * (-100)
rel_red_LPM2_wheat_25[,3] <- ((LPM2_wheat_25[,4] - LPM2_wheat_25[,1]) / LPM2_wheat_25[,1]) * (-100)
rel_red_LPM2_wheat_25[,4] <- ((LPM2_wheat_25[,5] - LPM2_wheat_25[,1]) / LPM2_wheat_25[,1]) * (-100)
rel_red_LPM2_wheat_25_melted <- melt(rel_red_LPM2_wheat_25)

wheat_25 <- ggplot(rel_red_LPM2_wheat_25_melted, aes(x=variable, y=value)) + ggtitle("b) winter wheat \n strike level: 25°C") +
  geom_boxplot(notch = F, fill="gray", outlier.size = 0.1, coef=0) + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "black", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "black"),
        axis.text = element_text(family="Times New Roman", size=10, color="black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  ylab("Reduction in downside variance (%)") + xlab("")+
  scale_y_continuous(breaks = c(-40,-20,0,20,40), limits = c(-40, 55))+
  geom_hline(yintercept=c(-40,-20,0,20,40), linetype="dashed", colour="gray", alpha = 0.5)+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20", alpha=0.6)+
  geom_text(data=label.df, aes(x=Var2, y=c(55,55,55,55)), label="***", size=4)

# 
# Plot LPM2 rapeseed 15
# 

label.df <- data.frame(Var2=c("Model 1","Model 2","Model 3","Model 4"), Value = c(45,45,45,45))
wilcox.test(LPM2_rapeseed_15[,2], LPM2_rapeseed_15[,1], paired=T, alternative="l")$p.value

rel_red_LPM2_rapeseed_15 <- as.data.frame(matrix(NA, nrow=nrow(detr_rapeseed_yields), ncol=4))
colnames(rel_red_LPM2_rapeseed_15) <- c("Model 1","Model 2", "Model 3", "Model 4")
rel_red_LPM2_rapeseed_15[,1] <- ((LPM2_rapeseed_15[,2] - LPM2_rapeseed_15[,1]) / LPM2_rapeseed_15[,1]) * (-100)
rel_red_LPM2_rapeseed_15[,2] <- ((LPM2_rapeseed_15[,3] - LPM2_rapeseed_15[,1]) / LPM2_rapeseed_15[,1]) * (-100)
rel_red_LPM2_rapeseed_15[,3] <- ((LPM2_rapeseed_15[,4] - LPM2_rapeseed_15[,1]) / LPM2_rapeseed_15[,1]) * (-100)
rel_red_LPM2_rapeseed_15[,4] <- ((LPM2_rapeseed_15[,5] - LPM2_rapeseed_15[,1]) / LPM2_rapeseed_15[,1]) * (-100)
rel_red_LPM2_rapeseed_15_melted <- melt(rel_red_LPM2_rapeseed_15)

rapeseed_15 <- ggplot(rel_red_LPM2_rapeseed_15_melted, aes(x=variable, y=value)) + ggtitle("c) winter rapeseed \n strike level: 15°C") +
  geom_boxplot(notch = F, fill="gray", outlier.size = 0.1, coef=0) + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "black", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "black"),
        axis.text = element_text(family="Times New Roman", size=10, color="black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  ylab("Reduction in downside variance (%)") + xlab("")+
  scale_y_continuous(breaks = c(-40,-20,0,20,40), limits = c(-40, 55))+
  geom_hline(yintercept=c(-40,-20,0,20,40), linetype="dashed", colour="gray", alpha = 0.5)+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20", alpha=0.6)+
  geom_text(data=label.df, aes(x=Var2, y=c(55,55,55,55)), label="***", size=4)

# 
# Plot LPM2 rapeseed 20
# 

label.df <- data.frame(Var2=c("Model 1","Model 2","Model 3","Model 4"), Value = c(45,45,45,45))
wilcox.test(LPM2_rapeseed_20[,5], LPM2_rapeseed_20[,1], paired=T, alternative="l")$p.value

rel_red_LPM2_rapeseed_20 <- as.data.frame(matrix(NA, nrow=nrow(detr_rapeseed_yields), ncol=4))
colnames(rel_red_LPM2_rapeseed_20) <- c("Model 1","Model 2", "Model 3", "Model 4")
rel_red_LPM2_rapeseed_20[,1] <- ((LPM2_rapeseed_20[,2] - LPM2_rapeseed_20[,1]) / LPM2_rapeseed_20[,1]) * (-100)
rel_red_LPM2_rapeseed_20[,2] <- ((LPM2_rapeseed_20[,3] - LPM2_rapeseed_20[,1]) / LPM2_rapeseed_20[,1]) * (-100)
rel_red_LPM2_rapeseed_20[,3] <- ((LPM2_rapeseed_20[,4] - LPM2_rapeseed_20[,1]) / LPM2_rapeseed_20[,1]) * (-100)
rel_red_LPM2_rapeseed_20[,4] <- ((LPM2_rapeseed_20[,5] - LPM2_rapeseed_20[,1]) / LPM2_rapeseed_20[,1]) * (-100)
rel_red_LPM2_rapeseed_20_melted <- melt(rel_red_LPM2_rapeseed_20)

rapeseed_20 <- ggplot(rel_red_LPM2_rapeseed_20_melted, aes(x=variable, y=value)) + ggtitle("d) winter rapeseed \n strike level: 20°C") +
  geom_boxplot(notch = F, fill="gray", outlier.size = 0.1, coef=0) + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(family="Times New Roman", size=12, colour = "black", hjust=0.5),
        axis.title = element_text(family="Times New Roman", size=10, colour = "black"),
        axis.text = element_text(family="Times New Roman", size=10, color="black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none")+
  ylab("Reduction in downside variance (%)") + xlab("")+
  scale_y_continuous(breaks = c(-40,-20,0,20,40), limits = c(-40, 55))+
  geom_hline(yintercept=c(-40,-20,0,20,40), linetype="dashed", colour="gray", alpha = 0.5)+
  geom_hline(yintercept=0, linetype="dashed", colour="grey20", alpha=0.6)+
  geom_text(data=label.df, aes(x=Var2, y=c(55,55,55,55)), label="***", size=4)

LPM2 <- ggarrange(wheat_20, wheat_25,rapeseed_15,rapeseed_20, nrow=2, ncol=2)