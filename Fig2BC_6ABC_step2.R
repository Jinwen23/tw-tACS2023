
# load packages
{library(lme4)
  library(tidyverse)
  library(afex)
  library(standardize)
  library(reshape2)
  library(ggeffects) 
  library(emmeans)
  library(sjPlot)
  library(sjmisc)
  library(optimx)
  library(ggplot2)
  library(see)
  library(dplyr)
  library(rstatix)
  library(parameters)
  library(brms)
  library(bayesplot)
  library(posterior)
  library(ggridges)
  library(bayestestR)}

{set_theme(
  base = theme_classic(),
  axis.tickslen = 0, # hides tick marks
  axis.title.size = .9,
  axis.textsize = .9,
  legend.size = .7,
  legend.title.size = .8,
  geom.label.size = 3.5)
  theme_set(theme_classic(base_size=15))}

# load data----

rm(list = ls())
analysis_path <- dirname(rstudioapi::getActiveDocumentContext()$path)

setwd(analysis_path)
Data <- read.csv("Exp1.csv", header = TRUE)
LongData <- Data
head(LongData)

# convert to factors
LongData[,colnames(LongData)[c(1,2)]] <- data.frame(lapply(LongData[,colnames(LongData)[c(1,2)]], factor))
head(LongData)

# select data
LongData1<-na.omit(LongData)
head(LongData1)
LongData_beh<-dplyr::filter(LongData1, thres_all>0)
head(LongData_beh)

# augment threshold scale
LongData_beh$thres_all <- 100*LongData_beh$thres_all

# preprocessing----
# fw
LongData_fw <- LongData_beh
# average within blocks
LongData_fw <- as.data.frame(LongData_fw %>%
                               group_by(id,block) %>%#
                               summarise_at(vars(FW,thres_all),funs(mean(.,na.rm=TRUE))))

# remove outlier
{
  subject_id_column <- 1
  unique_subject_ids <- unique(LongData_fw[, subject_id_column])
  LongData_fwClean <- data.frame()
  LongData_fwExclude <- data.frame()
  outliersn<-data.frame()
  # Loop through each unique subject ID
  for (subject_id in unique_subject_ids) {
    # Filter data for the current subject ID
    subject_data <- LongData_fw[LongData_fw[, subject_id_column] == subject_id, ]
    # calculate median
    median_val = median(subject_data$FW, na.rm = TRUE)
    # calculate IQR
    IQR_val = IQR(subject_data$FW, na.rm = TRUE)
    # calculate the upper and lower bounds for outliers
    lower_bound = median_val - 2 * IQR_val#2.5 or 3 median outlier
    upper_bound = median_val + 2 * IQR_val
    # get the outliers
    outliersn = subject_data$FW[subject_data$FW < lower_bound | subject_data$FW > upper_bound]
    
    if(length(outliersn) > 0) {
      subject_data_no_outliers <- subject_data[-which(subject_data$FW %in% outliersn), ]
    } else {
      subject_data_no_outliers <- subject_data
    }
 
    LongData_fwClean <- rbind(LongData_fwClean, subject_data_no_outliers)
    # Extract outliers and store the result in 'LongData_fwExclude' dataframe
    subject_data_outliers <- subject_data[which(subject_data$FW %in% outliersn), ]
    LongData_fwExclude <- rbind(LongData_fwExclude, subject_data_outliers)
  }
}

# center
LongData_fwClean$FW = scale_by(FW~id,LongData_fwClean,scale=0)


# bw
LongData_bw <- LongData_beh
# average within blocks
LongData_bw <- as.data.frame(LongData_bw %>%
                               group_by(id,block) %>%#
                               summarise_at(vars(thres_all,BW),funs(mean(.,na.rm=TRUE))))

# remove outliers
{
  subject_id_column <- 1
  unique_subject_ids <- unique(LongData_bw[, subject_id_column])
  LongData_bwClean <- data.frame()
  LongData_bwExclude <- data.frame()
  outliersn<-data.frame()
  # Loop through each unique subject ID
  for (subject_id in unique_subject_ids) {
    # Filter data for the current subject ID
    subject_data <- LongData_bw[LongData_bw[, subject_id_column] == subject_id, ]
    # calculate median
    median_val = median(subject_data$BW, na.rm = TRUE)
    # calculate IQR
    IQR_val = IQR(subject_data$BW, na.rm = TRUE)
    # calculate the upper and lower bounds for outliers
    lower_bound = median_val - 2 * IQR_val#2.5 or 3 median outlier
    upper_bound = median_val + 2 * IQR_val
    # get the outliers
    outliersn = subject_data$BW[subject_data$BW < lower_bound | subject_data$BW > upper_bound]
    
    if(length(outliersn) > 0) {
      subject_data_no_outliers <- subject_data[-which(subject_data$BW %in% outliersn), ]
    } else {
      subject_data_no_outliers <- subject_data
    }
    
    LongData_bwClean <- rbind(LongData_bwClean, subject_data_no_outliers)
    # Extract outliers and store the result in 'LongData_fwExclude' dataframe
    subject_data_outliers <- subject_data[which(subject_data$BW %in% outliersn), ]
    LongData_bwExclude <- rbind(LongData_bwExclude, subject_data_outliers)
  }
}
# center
LongData_bwClean$BW = scale_by(BW~id,LongData_bwClean,scale=0)

# ap
LongData_ap <- LongData_beh
# average within blocks
LongData_ap <- as.data.frame(LongData_ap %>%
                               group_by(id,block) %>%#
                               summarise_at(vars(AP,thres_all,FW, BW),funs(mean(.,na.rm=TRUE))))

# remove outliers
{
  subject_id_column <- 1
  unique_subject_ids <- unique(LongData_ap[, subject_id_column])
  LongData_apClean <- data.frame()
  LongData_apExclude <- data.frame()
  outliersn<-data.frame()
  # Loop through each unique subject ID
  for (subject_id in unique_subject_ids) {
    # Filter data for the current subject ID
    subject_data <- LongData_ap[LongData_ap[, subject_id_column] == subject_id, ]
    # calculate median
    median_val = median(subject_data$AP, na.rm = TRUE)
    # calculate IQR
    IQR_val = IQR(subject_data$AP, na.rm = TRUE)
    # calculate the upper and lower bounds for outliers
    lower_bound = median_val - 2 * IQR_val#2.5 or 3 median outlier
    upper_bound = median_val + 2 * IQR_val
    # get the outliers
    outliersn = subject_data$AP[subject_data$AP < lower_bound | subject_data$AP > upper_bound]
    
    if(length(outliersn) > 0) {
      subject_data_no_outliers <- subject_data[-which(subject_data$AP %in% outliersn), ]
    } else {
      subject_data_no_outliers <- subject_data
    }
    
    LongData_apClean <- rbind(LongData_apClean, subject_data_no_outliers)
    # Extract outliers and store the result in 'LongData_fwExclude' dataframe
    subject_data_outliers <- subject_data[which(subject_data$AP %in% outliersn), ]
    LongData_apExclude <- rbind(LongData_apExclude, subject_data_outliers)
  }
}
# center
LongData_apClean$AP = scale_by(AP~id,LongData_apClean,scale=0)

# FW traveling wave predicts visual threshold----

BehLMM_fw = brm(thres_all  ~ 1+FW+(1|id)+(1|block),
                data=LongData_fwClean,iter = 3000,cores = parallel::detectCores(),
                control=list(adapt_delta=0.95),seed=42)
p_direction(BehLMM_fw)
summary(BehLMM_fw)
# 95%CI
fixef(BehLMM_fw)
# 95% HDI
{posterior_samples <- as_draws_matrix(BehLMM_fw)
fw_data<-as.data.frame(posterior_samples[, c("b_FW")])
hpd_intervals <- lapply(fw_data, hdi)
lapply(hpd_intervals, function(x) {
  cat(paste("95% HDI: [", round(x$CI_low, 4), ", ", round(x$CI_high, 4), "]\n", sep=""))
})}

cond_eff<-conditional_effects(BehLMM_fw, effects = "FW", method = "posterior_epred")

# plot regression line
{# Extract prediction data
pred_data <- as.data.frame(cond_eff$FW)

# Extract relevant columns from original and prediction data
original_data <- LongData_fwClean[, c("FW", "thres_all")]
pred_data <- pred_data[, c("effect1__", "estimate__", "lower__", "upper__")]
names(pred_data)[names(pred_data) == "effect1__"] <- "FW"
names(pred_data)[names(pred_data) == "estimate__"] <- "thres_all"
print(names(pred_data))

# Add a new column to distinguish the data sources
original_data$source <- "Observed data"
pred_data$source <- "Regression line"
original_data$lower__ <- NA
original_data$upper__ <- NA

# Combine original and prediction data
combined_data <- rbind(original_data, pred_data)

plotsave<-ggplot(combined_data, aes(x=FW, y=thres_all)) +
  geom_point(data=combined_data[combined_data$source=="Observed data",], aes(color=source), alpha=0.3) +
  geom_line(data=combined_data[combined_data$source=="Regression line",], aes(color=source)) +
  geom_ribbon(data=combined_data[combined_data$source=="Regression line",],
              aes(ymin=lower__, ymax=upper__, fill=source), alpha=0.3) +
  scale_color_manual(values = c("Observed data" = "black", "Regression line" = "red"), labels = c("Observed data", "Regression line (95% HDI)")) +
  scale_fill_manual(values = c("Observed data" = NA, "Regression line" = "#69b3a2")) +
  guides(color = guide_legend(title = NULL, 
                              override.aes = list(shape=c(16, NA), linetype=c("blank", "solid"), color=c("black", "red"))),
         fill = FALSE) +
  theme(plot.title=element_text(hjust=0.5,size = 13),legend.position = "top",
        axis.text.x = element_text(size = 13),axis.title.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),axis.title.y = element_text(size = 13),
        legend.text = element_text(size=13),legend.title = element_text(size=13)) +
  xlab("fTW (dB)") +
  ylab("Fitted VCT (%)")
plotsave
}

# calculate posterior probability
prob_different_from_zero <- mean(fw_data < 0)
prob_different_from_zero

# BW traveling wave predicts visual threshold----

BehLMM_bw = brm(thres_all  ~ 1+BW+(1|id)+(1|block),
              data=LongData_bwClean,iter = 3000,cores = parallel::detectCores(),control=list(adapt_delta=0.95),seed=42)
summary(BehLMM_bw)
p_direction(BehLMM_bw)
# 95%CI
fixef(BehLMM_bw)
# 95% HDI
{posterior_samples <- as_draws_matrix(BehLMM_bw)
  bw_data<-as.data.frame(posterior_samples[, c("b_BW")])
  hpd_intervals <- lapply(bw_data, hdi)
  lapply(hpd_intervals, function(x) {
    cat(paste("95% HDI: [", round(x$CI_low, 4), ", ", round(x$CI_high, 4), "]\n", sep=""))
  })}

cond_eff<-conditional_effects(BehLMM_bw, effects = "BW", method = "posterior_epred")
# plot regression line

{# Extract prediction data
  pred_data <- as.data.frame(cond_eff$BW)
  
  # Extract relevant columns from original and prediction data
  original_data <- LongData_bwClean[, c("BW", "thres_all")]
  pred_data <- pred_data[, c("effect1__", "estimate__", "lower__", "upper__")]
  names(pred_data)[names(pred_data) == "effect1__"] <- "BW"
  names(pred_data)[names(pred_data) == "estimate__"] <- "thres_all"
  print(names(pred_data))
  
  # Add a new column to distinguish the data sources
  original_data$source <- "Observed data"
  pred_data$source <- "Regression line"
  original_data$lower__ <- NA
  original_data$upper__ <- NA
  
  # Combine original and prediction data
  combined_data <- rbind(original_data, pred_data)
  
  plotsave<-ggplot(combined_data, aes(x=BW, y=thres_all)) +
    geom_point(data=combined_data[combined_data$source=="Observed data",], aes(color=source), alpha=0.3) +
    geom_line(data=combined_data[combined_data$source=="Regression line",], aes(color=source)) +
    geom_ribbon(data=combined_data[combined_data$source=="Regression line",],
                aes(ymin=lower__, ymax=upper__, fill=source), alpha=0.3) +
    scale_color_manual(values = c("Observed data" = "black", "Regression line" = "red"), labels = c("Observed data", "Regression line (95% HDI)")) +
    scale_fill_manual(values = c("Observed data" = NA, "Regression line" = "#69b3a2")) +
    guides(color = guide_legend(title = NULL, 
                                override.aes = list(shape=c(16, NA), linetype=c("blank", "solid"), color=c("black", "red"))),
           fill = FALSE) +
    theme(plot.title=element_text(hjust=0.5,size = 13),legend.position = "top",
          axis.text.x = element_text(size = 13),axis.title.x = element_text(size = 13),
          axis.text.y = element_text(size = 13),axis.title.y = element_text(size = 13),
          legend.text = element_text(size=13),legend.title = element_text(size=13)) +
    xlab("bTW (dB)") +
    ylab("Fitted VCT (%)")
  plotsave
  
}



# calculate posterior probability
prob_different_from_zero <- mean(bw_data > 0)
prob_different_from_zero


# AP predicts visual threshold----

BehLMM_ap = brm(thres_all  ~ 1+AP+(1|id)+(1|block),
                data=LongData_apClean,iter = 3000,cores = parallel::detectCores(),control=list(adapt_delta=0.95),seed=42)
summary(BehLMM_ap)
p_direction(BehLMM_ap)
# 95%CI
fixef(BehLMM_ap)
# 95% HDI
{posterior_samples <- as_draws_matrix(BehLMM_ap)
  ap_data<-as.data.frame(posterior_samples[, c("b_AP")])
  hpd_intervals <- lapply(ap_data, hdi)
  lapply(hpd_intervals, function(x) {
    cat(paste("95% HDI: [", round(x$CI_low, 4), ", ", round(x$CI_high, 4), "]\n", sep=""))
  })}

cond_eff<-conditional_effects(BehLMM_ap, effects = "AP", method = "posterior_epred")

# plot regression line
{# Extract prediction data
  pred_data <- as.data.frame(cond_eff$AP)
  
  # Extract relevant columns from original and prediction data
  original_data <- LongData_apClean[, c("AP", "thres_all")]
  pred_data <- pred_data[, c("effect1__", "estimate__", "lower__", "upper__")]
  names(pred_data)[names(pred_data) == "effect1__"] <- "AP"
  names(pred_data)[names(pred_data) == "estimate__"] <- "thres_all"
  print(names(pred_data))
  
  # Add a new column to distinguish the data sources
  original_data$source <- "Observed data"
  pred_data$source <- "Regression line"
  original_data$lower__ <- NA
  original_data$upper__ <- NA
  
  # Combine original and prediction data
  combined_data <- rbind(original_data, pred_data)
  
  plotsave<-ggplot(combined_data, aes(x=AP, y=thres_all)) +
    geom_point(data=combined_data[combined_data$source=="Observed data",], aes(color=source), alpha=0.3) +
    geom_line(data=combined_data[combined_data$source=="Regression line",], aes(color=source)) +
    geom_ribbon(data=combined_data[combined_data$source=="Regression line",],
                aes(ymin=lower__, ymax=upper__, fill=source), alpha=0.3) +
    scale_color_manual(values = c("Observed data" = "black", "Regression line" = "red"), labels = c("Observed data", "Regression line (95% HDI)")) +
    scale_fill_manual(values = c("Observed data" = NA, "Regression line" = "#69b3a2")) +
    guides(color = guide_legend(title = NULL, 
                                override.aes = list(shape=c(16, NA), linetype=c("blank", "solid"), color=c("black", "red"))),
           fill = FALSE) +
    theme(plot.title=element_text(hjust=0.5,size = 13),legend.position = "top",
          axis.text.x = element_text(size = 13),axis.title.x = element_text(size = 13),
          axis.text.y = element_text(size = 13),axis.title.y = element_text(size = 13),
          legend.text = element_text(size=13),legend.title = element_text(size=13)) +
    xlab("AP (dB)") +
    ylab("Fitted VCT (%)")
  plotsave
}

# calculate posterior probability
prob_different_from_zero <- mean(bw_data > 0)
prob_different_from_zero


# AP predicts forward traveling waves----

BehLMM_apfw = brm(FW  ~ 1+AP+(1|id)+(1|block),
                data=LongData_apClean,iter = 3000,cores = parallel::detectCores(),control=list(adapt_delta=0.95),seed=42)
summary(BehLMM_apfw)
p_direction(BehLMM_apfw)
# 95%CI
fixef(BehLMM_apfw)
# 95% HDI
{posterior_samples <- as_draws_matrix(BehLMM_apfw)
  apfw_data<-as.data.frame(posterior_samples[, c("b_AP")])
  hpd_intervals <- lapply(apfw_data, hdi)
  lapply(hpd_intervals, function(x) {
    cat(paste("95% HDI: [", round(x$CI_low, 4), ", ", round(x$CI_high, 4), "]\n", sep=""))
  })}

cond_eff<-conditional_effects(BehLMM_apfw, effects = "AP", method = "posterior_epred")
# plot regression line

{# Extract prediction data
  pred_data <- as.data.frame(cond_eff$AP)
  
  # Extract relevant columns from original and prediction data
  original_data <- LongData_apClean[, c("AP", "FW")]
  pred_data <- pred_data[, c("effect1__", "estimate__", "lower__", "upper__")]
  names(pred_data)[names(pred_data) == "effect1__"] <- "AP"
  names(pred_data)[names(pred_data) == "estimate__"] <- "FW"
  print(names(pred_data))
  
  # Add a new column to distinguish the data sources
  original_data$source <- "Observed data"
  pred_data$source <- "Regression line"
  original_data$lower__ <- NA
  original_data$upper__ <- NA
  
  # Combine original and prediction data
  combined_data <- rbind(original_data, pred_data)
  
  plotsave<-ggplot(combined_data, aes(x=AP, y=FW)) +
    geom_point(data=combined_data[combined_data$source=="Observed data",], aes(color=source), alpha=0.3) +
    geom_line(data=combined_data[combined_data$source=="Regression line",], aes(color=source)) +
    geom_ribbon(data=combined_data[combined_data$source=="Regression line",],
                aes(ymin=lower__, ymax=upper__, fill=source), alpha=0.3) +
    scale_color_manual(values = c("Observed data" = "black", "Regression line" = "red"), labels = c("Observed data", "Regression line (95% HDI)")) +
    scale_fill_manual(values = c("Observed data" = NA, "Regression line" = "#69b3a2")) +
    guides(color = guide_legend(title = NULL, 
                                override.aes = list(shape=c(16, NA), linetype=c("blank", "solid"), color=c("black", "red"))),
           fill = FALSE) +
    theme(plot.title=element_text(hjust=0.5,size = 13),legend.position = "top",
          axis.text.x = element_text(size = 13),axis.title.x = element_text(size = 13),
          axis.text.y = element_text(size = 13),axis.title.y = element_text(size = 13),
          legend.text = element_text(size=13),legend.title = element_text(size=13)) +
    xlab("AP (dB)") +
    ylab("Fitted fTW (dB)")
  plotsave

}

# calculate posterior probability
prob_different_from_zero <- mean(bw_data > 0)
prob_different_from_zero

# AP predicts backward traveling waves----

BehLMM_apbw = brm(BW  ~ 1+AP+(1|id)+(1|block),
                  data=LongData_apClean,iter = 3000,cores = parallel::detectCores(),control=list(adapt_delta=0.95),seed=42)
summary(BehLMM_apbw)
p_direction(BehLMM_apbw)
# 95%CI
fixef(BehLMM_apbw)
# 95% HDI
{posterior_samples <- as_draws_matrix(BehLMM_apbw)
  apbw_data<-as.data.frame(posterior_samples[, c("b_AP")])
  hpd_intervals <- lapply(apbw_data, hdi)
  lapply(hpd_intervals, function(x) {
    cat(paste("95% HDI: [", round(x$CI_low, 4), ", ", round(x$CI_high, 4), "]\n", sep=""))
  })}
cond_eff<-conditional_effects(BehLMM_apbw, effects = "AP", method = "posterior_epred")
# plot regression line

{# Extract prediction data
  pred_data <- as.data.frame(cond_eff$AP)
  
  # Extract relevant columns from original and prediction data
  original_data <- LongData_apClean[, c("AP", "BW")]
  pred_data <- pred_data[, c("effect1__", "estimate__", "lower__", "upper__")]
  names(pred_data)[names(pred_data) == "effect1__"] <- "AP"
  names(pred_data)[names(pred_data) == "estimate__"] <- "BW"
  print(names(pred_data))
  
  # Add a new column to distinguish the data sources
  original_data$source <- "Observed data"
  pred_data$source <- "Regression line"
  original_data$lower__ <- NA
  original_data$upper__ <- NA
  
  # Combine original and prediction data
  combined_data <- rbind(original_data, pred_data)
  
  plotsave<-ggplot(combined_data, aes(x=AP, y=BW)) +
    geom_point(data=combined_data[combined_data$source=="Observed data",], aes(color=source), alpha=0.3) +
    geom_line(data=combined_data[combined_data$source=="Regression line",], aes(color=source)) +
    geom_ribbon(data=combined_data[combined_data$source=="Regression line",],
                aes(ymin=lower__, ymax=upper__, fill=source), alpha=0.3) +
    scale_color_manual(values = c("Observed data" = "black", "Regression line" = "red"), labels = c("Observed data", "Regression line (95% HDI)")) +
    scale_fill_manual(values = c("Observed data" = NA, "Regression line" = "#69b3a2")) +
    guides(color = guide_legend(title = NULL, 
                                override.aes = list(shape=c(16, NA), linetype=c("blank", "solid"), color=c("black", "red"))),
           fill = FALSE) +
    theme(plot.title=element_text(hjust=0.5,size = 13),legend.position = "top",
          axis.text.x = element_text(size = 13),axis.title.x = element_text(size = 13),
          axis.text.y = element_text(size = 13),axis.title.y = element_text(size = 13),
          legend.text = element_text(size=13),legend.title = element_text(size=13)) +
    xlab("AP (dB)") +
    ylab("Fitted bTW (dB)")
  plotsave
 
}

