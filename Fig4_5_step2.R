
{
# load packages
library(lme4)
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
library(tidybayes)
library("RColorBrewer")}
{set_theme(
  base = theme_classic(),
  axis.tickslen = 0, # hides tick marks
  axis.title.size = .9,
  axis.textsize = .9,
  legend.size = .7,
  legend.title.size = .8,
  geom.label.size = 3.5
)}

library(tidyverse)
library(ggridges)
library(bayestestR)
theme_set(theme_classic(base_size=15))


# load and preprocess data ----

# load data
{rm(list = ls())
Data <- read.csv("Exp2.csv", header = TRUE)
LongData <- Data
head(LongData)

# preprocessing
LongData[,colnames(LongData)[c(1,2,3,4)]] <- data.frame(lapply(LongData[,colnames(LongData)[c(1,2,3,4)]], factor))
LongData$group=factor(LongData$group,levels=c(3,2,1),labels=c("Sham","-45°","45°"))
LongData$interval = factor(LongData$interval,levels=c(1,2,3),labels=c("pre-sti","during-sti","post-sti"))
LongData$thres <-100*LongData$thres
LongData$pre_thres <-100*LongData$pre_thres
LongData_TW<-na.omit(LongData)
head(LongData_TW)

# prepare data----

# behavioral data

# for fatigue
LongData_beh_all<-dplyr::filter(LongData,
                                blockIndex!=9&blockIndex!=1&interval!="pre-sti") # all data 
{# block average
  LongData_beh_all <- as.data.frame(LongData_beh_all %>%
                                      group_by(id,group,blockIndex,interval) %>%#
                                      summarise_at(vars(thres,pre_thres,block_rating,pre_rating),funs(mean(.,na.rm=TRUE))))
}
# for behavioral analysis
LongData_beh<-dplyr::filter(LongData,
                            thres>0&blockIndex!=9&blockIndex!=1&interval!="pre-sti") # all data 
{# block average
  LongData_beh <- as.data.frame(LongData_beh %>%
                                      group_by(id,group,blockIndex,interval) %>%#
                                      summarise_at(vars(thres,pre_thres,block_rating),funs(mean(.,na.rm=TRUE))))
}


# remove outliers
{
  subject_id_column <- 1
  unique_subject_ids <- unique(LongData_beh[, subject_id_column])
  LongData_Clean <- data.frame()
  LongData_Exclude <- data.frame()
  outliersn<-data.frame()
  # Loop through each unique subject ID
  for (subject_id in unique_subject_ids) {
    # Filter data for the current subject ID
    subject_data <- LongData_beh[LongData_beh[, subject_id_column] == subject_id, ]
    # calculate median
    median_val = median(subject_data$thres, na.rm = TRUE)
    # calculate IQR
    IQR_val = IQR(subject_data$thres, na.rm = TRUE)
    # calculate the upper and lower bounds for outliers
    lower_bound = median_val - 2 * IQR_val#2.5 or 3 median outlier
    upper_bound = median_val + 2 * IQR_val
    # get the outliers
    outliersn = subject_data$thres[subject_data$thres < lower_bound | subject_data$thres > upper_bound]
    
    if(length(outliersn) > 0) {
      subject_data_no_outliers <- subject_data[-which(subject_data$thres %in% outliersn), ]
    } else {
      subject_data_no_outliers <- subject_data
    }
    
    LongData_Clean <- rbind(LongData_Clean, subject_data_no_outliers)
    # Extract outliers and store the result in 'LongData_fwExclude' dataframe
    subject_data_outliers <- subject_data[which(subject_data$thres %in% outliersn), ]
    LongData_Exclude <- rbind(LongData_Exclude, subject_data_outliers)
  }
}
LongData_beh<-LongData_Clean
LongData_beh_low<-dplyr::filter(LongData_beh,
                                block_rating<4 ) # low-fatigue state
LongData_beh_high<-dplyr::filter(LongData_beh,
                                 block_rating>=4) # high-fatigue state

# behavioral data after stimulation
LongData_beh_post<-dplyr::filter(LongData_beh,interval=="post-sti")
LongData_beh_postlow<-dplyr::filter(LongData_beh_low,interval=="post-sti")
LongData_beh_posthigh<-dplyr::filter(LongData_beh_high,interval=="post-sti")

# TW data after stimulation

LongData_TW_post<-dplyr::filter(LongData_TW,interval!="pre-sti"&blockIndex!=9&blockIndex!=1
) # all TW data in post stimulation

{# block average
  LongData_TW_post <- as.data.frame(LongData_TW_post %>%
                                  group_by(id,group,blockIndex,interval) %>%#
                                  summarise_at(vars(thres,pre_thres,block_rating,pre_BW,BW_TW,pre_FW,FW_TW,pre_AP,AlphaPower),funs(mean(.,na.rm=TRUE))))
}

LongData_TW_low<-dplyr::filter(LongData_TW_post,
                               block_rating<4 ) # post-sti data in low-fatigue state
LongData_TW_high<-dplyr::filter(LongData_TW_post,
                                block_rating>=4) # post-sti data in high-fatigue state

# fw
LongData_FW_post<-LongData_TW_post
# remove outliers
{
  subject_id_column <- 1
  unique_subject_ids <- unique(LongData_FW_post[, subject_id_column])
  LongData_Clean <- data.frame()
  LongData_Exclude <- data.frame()
  outliersn<-data.frame()
  # Loop through each unique subject ID
  for (subject_id in unique_subject_ids) {
    # Filter data for the current subject ID
    subject_data <- LongData_FW_post[LongData_FW_post[, subject_id_column] == subject_id, ]
    # calculate median
    median_val = median(subject_data$FW_TW, na.rm = TRUE)
    # calculate IQR
    IQR_val = IQR(subject_data$FW_TW, na.rm = TRUE)
    # calculate the upper and lower bounds for outliers
    lower_bound = median_val - 2* IQR_val#2.5 or 3 median outlier
    upper_bound = median_val + 2 * IQR_val
    # get the outliers
    outliersn = subject_data$FW_TW[subject_data$FW_TW < lower_bound | subject_data$FW_TW > upper_bound]
    
    if(length(outliersn) > 0) {
      subject_data_no_outliers <- subject_data[-which(subject_data$FW_TW %in% outliersn), ]
    } else {
      subject_data_no_outliers <- subject_data
    }
    
    LongData_Clean <- rbind(LongData_Clean, subject_data_no_outliers)
    # Extract outliers and store the result in 'LongData_fwExclude' dataframe
    subject_data_outliers <- subject_data[which(subject_data$FW_TW %in% outliersn), ]
    LongData_Exclude <- rbind(LongData_Exclude, subject_data_outliers)
  }
}
LongData_FW<-LongData_Clean
LongData_FW_low<-dplyr::filter(LongData_FW,
                                block_rating<4 ) # post-sti data in low-fatigue state
LongData_FW_high<-dplyr::filter(LongData_FW,
                                 block_rating>=4) # post-sti data in high-fatigue state


# bw
LongData_BW_post<-LongData_TW_post
# remove outliers
{
  subject_id_column <- 1
  unique_subject_ids <- unique(LongData_BW_post[, subject_id_column])
  LongData_Clean <- data.frame()
  LongData_Exclude <- data.frame()
  outliersn<-data.frame()
  # Loop through each unique subject ID
  for (subject_id in unique_subject_ids) {
    # Filter data for the current subject ID
    subject_data <- LongData_BW_post[LongData_BW_post[, subject_id_column] == subject_id, ]
    # calculate median
    median_val = median(subject_data$BW_TW, na.rm = TRUE)
    # calculate IQR
    IQR_val = IQR(subject_data$BW_TW, na.rm = TRUE)
    # calculate the upper and lower bounds for outliers
    lower_bound = median_val - 2 * IQR_val#2.5 or 3 median outlier
    upper_bound = median_val + 2 * IQR_val
    # get the outliers
    outliersn = subject_data$BW_TW[subject_data$BW_TW < lower_bound | subject_data$BW_TW > upper_bound]
    
    if(length(outliersn) > 0) {
      subject_data_no_outliers <- subject_data[-which(subject_data$BW_TW %in% outliersn), ]
    } else {
      subject_data_no_outliers <- subject_data
    }
    
    LongData_Clean <- rbind(LongData_Clean, subject_data_no_outliers)
    # Extract outliers and store the result in 'LongData_fwExclude' dataframe
    subject_data_outliers <- subject_data[which(subject_data$BW_TW %in% outliersn), ]
    LongData_Exclude <- rbind(LongData_Exclude, subject_data_outliers)
  }
}
LongData_BW<-LongData_Clean
LongData_BW_low<-dplyr::filter(LongData_BW,
                               block_rating<4 ) # post-sti data in low-fatigue state
LongData_BW_high<-dplyr::filter(LongData_BW,
                                block_rating>=4) # post-sti data in high-fatigue state


# ap
LongData_AP_post<-LongData_TW_post
# remove outliers
{
  subject_id_column <- 1
  unique_subject_ids <- unique(LongData_AP_post[, subject_id_column])
  LongData_Clean <- data.frame()
  LongData_Exclude <- data.frame()
  outliersn<-data.frame()
  # Loop through each unique subject ID
  for (subject_id in unique_subject_ids) {
    # Filter data for the current subject ID
    subject_data <- LongData_AP_post[LongData_AP_post[, subject_id_column] == subject_id, ]
    # calculate median
    median_val = median(subject_data$AlphaPower, na.rm = TRUE)
    # calculate IQR
    IQR_val = IQR(subject_data$AlphaPower, na.rm = TRUE)
    # calculate the upper and lower bounds for outliers
    lower_bound = median_val - 2 * IQR_val#2.5 or 3 median outlier
    upper_bound = median_val + 2 * IQR_val
    # get the outliers
    outliersn = subject_data$AlphaPower[subject_data$AlphaPower < lower_bound | subject_data$AlphaPower > upper_bound]
    
    if(length(outliersn) > 0) {
      subject_data_no_outliers <- subject_data[-which(subject_data$AlphaPower %in% outliersn), ]
    } else {
      subject_data_no_outliers <- subject_data
    }
    
    LongData_Clean <- rbind(LongData_Clean, subject_data_no_outliers)
    # Extract outliers and store the result in 'LongData_fwExclude' dataframe
    subject_data_outliers <- subject_data[which(subject_data$AlphaPower %in% outliersn), ]
    LongData_Exclude <- rbind(LongData_Exclude, subject_data_outliers)
  }
}
LongData_AP<-LongData_Clean

LongData_AP_low<-dplyr::filter(LongData_AP,
                               block_rating<4 ) # post-sti data in low-fatigue state
LongData_AP_high<-dplyr::filter(LongData_AP,
                                block_rating>=4) # post-sti data in high-fatigue state
}

# Fig.3DEF: during- and after-effect on VCT ----

# frequentist LMM
Beh_ave = lmer(thres  ~ 1+pre_thres+group*blockIndex+(1+group|id),
               data=LongData_beh_ave_high,REML = FALSE,
               control = lmerControl(
                 optimizer ='optimx', optCtrl=list(method='nlminb')))
anova(Beh_ave)  
emmeans(Beh_ave,pairwise~group|blockIndex)
plot_model(Beh_ave,type="pre",terms=c("blockIndex","group"),show.data=FALSE, title = "")


# bayesian LMM
Beh_all = brm(thres  ~ 1+pre_thres+group*blockIndex+(1+group|id),
              data=LongData_beh,iter = 3000,cores = parallel::detectCores(),control=list(adapt_delta=0.95),seed=42)
Beh_low = brm(thres  ~ 1+pre_thres+group*blockIndex+(1+group|id),
              data=LongData_beh_low,iter = 3000,cores = parallel::detectCores(),control=list(adapt_delta=0.95),seed=42)
Beh_high = brm(thres  ~ 1+pre_thres+group*blockIndex+(1+group|id),
               data=LongData_beh_high,iter = 3000,cores = parallel::detectCores(),control=list(adapt_delta=0.95),seed=42)

Beh_ave <- Beh_all
Beh_ave <- Beh_low
Beh_ave <- Beh_high
summary(Beh_ave)
p_direction(Beh_ave)
fixef(Beh_ave)
emmeans(Beh_ave,pairwise~group|blockIndex)

# save statistics
{emm <-emmeans(Beh_ave,pairwise~group)
  pairwise_comparisons <- emm$contrasts
  stattable<-as.data.frame(pairwise_comparisons)
  library(bruceR)
  print_table(stattable, file="Beh_high.doc")}

# plot1
{me <- conditional_effects(Beh_ave, "blockIndex:group")
  
  plotsave<-plot(me, plot = FALSE)[[1]] +
    theme(plot.title=element_text(hjust=0.5,size = 18),legend.position = "top",
          axis.text.x = element_text(size = 13),axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 13),axis.title.y = element_text(size = 15),
          legend.text = element_text(size=13),legend.title = element_blank())+
    xlab("Block")+ylab("Fitted VCT (%)")  +
    scale_color_manual(values = c("45°"= "#c73e11", 
                                  "-45°"= "#e1a636", 
                                  "Sham"= "#1b72b3")) +
    scale_x_discrete(limits=c( "3", "4", "5", "6", "7","8"),labels = c("2", "3", "4", "5", "6", "7")) +
    scale_y_continuous(limits = c(19.5, 23.5))  # Adjust the y-axis limits here

  plotsave <- plotsave + geom_line(aes(group = group),  position = position_dodge(width = 0.35),size = 0.5)
  print(plotsave)
  
}

# calculate posterior pr|obability
{
  warp_em <- emmeans (Beh_ave,  ~ group|blockIndex,epred = TRUE)%>% 
    contrast(method = "pairwise") %>% 
    gather_emmeans_draws()
  intere <- warp_em$contrast=='Sham - 45°'&warp_em$blockIndex == '3'
  p<-mean(warp_em$.value[intere]>0)
  p}


# Fig 4A: fatigue levels----

ratings = brm(block_rating  ~ 1+pre_rating+group*blockIndex+(1+group|id),
          data=LongData_beh_all,iter = 3000,cores = parallel::detectCores(),control=list(adapt_delta=0.95),seed=42)
Beh_ave <- ratings

summary(Beh_ave)
p_direction(Beh_ave)
fixef(Beh_ave)
emmeans(Beh_ave,pairwise~group|blockIndex)

# save statistics
{emm <-emmeans(Beh_ave,pairwise~group|blockIndex)
  pairwise_comparisons <- emm$contrasts
  stattable<-as.data.frame(pairwise_comparisons)
  library(bruceR)
  print_table(stattable, file="rating.doc")}

# plot1
{me <- conditional_effects(Beh_ave, "blockIndex:group")
  
  plotsave<-plot(me, plot = FALSE)[[1]] +
    theme(plot.title=element_text(hjust=0.5,size = 18),legend.position = "top",
          axis.text.x = element_text(size = 13),axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 13),axis.title.y = element_text(size = 15),
          legend.text = element_text(size=13),legend.title = element_blank())+
    xlab("Block")+ylab("Fitted fatigue rating")  +
    scale_color_manual(values = c("45°"= "#c73e11", 
                                  "-45°"= "#e1a636", 
                                  "Sham"= "#1b72b3")) +
    scale_x_discrete(limits=c( "3", "4", "5", "6", "7","8"),labels = c("2", "3", "4", "5", "6", "7"))+
    scale_y_continuous(limits = c(2, 4))
  plotsave <- plotsave + geom_line(aes(group = group),  position = position_dodge(width = 0.35),size = 0.5)
  #ggtitle("Fatigue ratings in three stimulation groups")
  print(plotsave)
}

# calculate posterior pr|obability
{
  warp_em <- emmeans (Beh_ave,  ~ group|blockIndex,epred = TRUE)%>% 
    contrast(method = "pairwise") %>% 
    gather_emmeans_draws()
  intere <- warp_em$contrast=='Sham - (-45°)'&warp_em$blockIndex == '8'
  p<-mean(warp_em$.value[intere]<0)
  p}

# Fig.3BC: split data----
# x-axis: fatigue
plotsave <- ggplot(LongData_beh, aes(x = block_rating)) +
  geom_bar() +
  labs(x = "Fatigue rating", y = "Frequency")+
  theme(plot.title = element_text(hjust = 0.5, size=12)) +
  scale_x_continuous(breaks = seq(1, 7, by = 1))  # Set x-axis breaks from 1 to 7 with an increment of 1


# Fig.5CFI: post-stimulation VCT ----

# baye lmm
be_all_post = brm(thres  ~ 1+pre_thres+group+(1+group|id),
                  data=LongData_beh_post,iter = 3000,cores = parallel::detectCores(),control=list(adapt_delta=0.95),seed=42)

be_low_post = brm(thres  ~ 1+pre_thres+group+(1+group|id),
                  data=LongData_beh_postlow,iter = 3000,cores = parallel::detectCores(),control=list(adapt_delta=0.95),seed=42)
be_high_post = brm(thres  ~ 1+pre_thres+group*blockIndex+(1+group|id),
                   data=LongData_beh_posthigh,iter = 3000,cores = parallel::detectCores(),control=list(adapt_delta=0.95))

be_baye <-be_all_post
be_baye <-be_low_post
be_baye <-be_high_post
summary(be_baye)
p_direction(be_baye)
fixef(be_baye)
emmeans(be_baye,pairwise~group)

{emm <-emmeans(be_baye,pairwise~group)
  pairwise_comparisons <- emm$contrasts
  stattable<-as.data.frame(pairwise_comparisons)
  library(bruceR)
  print_table(stattable, file="be_baye_low.doc")}


# plot 1
{
  #this is to plot conditional effects
  cond_eff<-conditional_effects(be_baye, effects = "group", method = "posterior_epred")
  plot_data <- cond_eff$group[, c("group", "estimate__", "lower__", "upper__")]
  plot_save<-ggplot(plot_data, aes(x = group, y = estimate__, color = group)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lower__, ymax = upper__), width = 0.2, size = 2) +
    scale_color_manual(values = c("Sham" = "#1b72b3", "-45°" = "#e1a636", "45°" = "#c73e11")) +
    theme(plot.title=element_text(hjust=0.5,size = 13),legend.position = "none",
          axis.text.x = element_text(size = 13),axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 13),axis.title.y = element_text(size = 15),
          legend.text = element_text(size=13),legend.title = element_text(size=13)) +
    xlab("Stimulation condition") +
    ylab("Fitted VCT (%)")+
    scale_y_continuous(limits = c(19.5, 23.5))  # Adjust the y-axis limits here
  
  print(plot_save)
 
}


warp_em <- emmeans (be_baye,  ~ group,epred = TRUE)%>% 
  contrast(method = "pairwise") %>% 
  gather_emmeans_draws()

intere <- warp_em$contrast=='Sham - (-45°)'
p<-mean(warp_em$.value[intere]>0)
p

# Fig.5BFH: post-stimulation forward spontaneous traveling waves----

# frequentist LMM
FW_ave = lmer(FW_TW  ~ 1+pre_FW+group+(1+group|id),
               data=LongData_FW_postlow_ave,REML = FALSE,
               control = lmerControl(
                 optimizer ='optimx', optCtrl=list(method='nlminb')))
anova(FW_ave)  
emmeans(FW_ave,pairwise~group)
plot_model(FW_ave,type="pre",terms=c("group"),show.data=FALSE, title = "")


# baye lmm
FW_all_baye = brm(FW_TW  ~ 1+pre_FW+group+(1+group|id),
                  data=LongData_FW,iter = 2000,cores = parallel::detectCores(),control=list(adapt_delta=0.95),seed=42)

FW_low_baye = brm(FW_TW  ~ 1+pre_FW+group+(1+group|id),
                  data=LongData_FW_low,iter = 3000,cores = parallel::detectCores(),control=list(adapt_delta=0.95),seed=42)

FW_high_baye = brm(FW_TW  ~ 1+pre_FW+group+(1+group|id),
                   data=LongData_FW_high,iter = 3000,cores = parallel::detectCores(),control=list(adapt_delta=0.95))
FW_baye <-FW_all_baye
FW_baye <-FW_low_baye
FW_baye <-FW_high_baye
summary(FW_baye)
p_direction(FW_baye)
fixef(FW_baye)

{emm <-emmeans(FW_baye,pairwise~group)
  pairwise_comparisons <- emm$contrasts
  stattable<-as.data.frame(pairwise_comparisons)
  library(bruceR)
  print_table(stattable, file="FW_baye_high.doc")}

# 95% HDI
{posterior_samples <- as_draws_matrix(FW_baye)
  ap_data<-as.data.frame(posterior_samples[, c("b_Intercept")])
  hpd_intervals <- lapply(ap_data, hdi)
  lapply(hpd_intervals, function(x) {
    cat(paste("95% HDI: [", round(x$CI_low, 4), ", ", round(x$CI_high, 4), "]\n", sep=""))
  })}
emmeans(FW_baye,pairwise~group)

# plot 1
{
  #this is to plot conditional effects
  cond_eff<-conditional_effects(FW_baye, effects = "group", method = "posterior_epred")
  plot_data <- cond_eff$group[, c("group", "estimate__", "lower__", "upper__")]
  plot_save<-ggplot(plot_data, aes(x = group, y = estimate__, color = group)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lower__, ymax = upper__), width = 0.2, size = 2) +
    scale_color_manual(values = c("Sham" = "#1b72b3", "-45°" = "#e1a636", "45°" = "#c73e11")) +
    theme(plot.title=element_text(hjust=0.5,size = 13),legend.position = "none",
          axis.text.x = element_text(size = 13),axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 13),axis.title.y = element_text(size = 15),
          legend.text = element_text(size=13),legend.title = element_text(size=13)) +
    xlab("Stimulation condition") +
    ylab("Fitted fTW (dB)")+
    scale_y_continuous(limits = c(0.35, 0.8))  # Adjust the y-axis limits here
  
  print(plot_save)
 
}

# posterior probability
{warp_em <- emmeans (FW_baye,  ~ group,epred = TRUE)%>% 
  contrast(method = "pairwise") %>% 
  gather_emmeans_draws()

intere <- warp_em$contrast=='Sham - (-45°)'
p<-mean(warp_em$.value[intere]<0)
p}


# Fig.5ADG: post-stimulation backward spontaneous traveling waves----
# baye lmm
BW_all_baye = brm(BW_TW  ~ 1+pre_BW+group+(1+group|id),
                  data=LongData_BW,iter = 3000,cores = parallel::detectCores(),control=list(adapt_delta=0.95),seed=42)

BW_low_baye = brm(BW_TW  ~ 1+pre_BW+group+(1+group|id),
                  data=LongData_BW_low,iter = 3000,cores = parallel::detectCores(),control=list(adapt_delta=0.95),seed=42)

BW_high_baye = brm(BW_TW  ~ 1+pre_BW+group+(1+group|id),
                   data=LongData_BW_high,iter = 3000,cores = parallel::detectCores(),control=list(adapt_delta=0.95))

BW_baye <-BW_all_baye
BW_baye <-BW_low_baye
BW_baye <-BW_high_baye
summary(BW_baye)
p_direction(BW_baye)
fixef(BW_baye)
emmeans(BW_baye,pairwise~group)

# save statistics
{emm <-emmeans(BW_baye,pairwise~group)
  pairwise_comparisons <- emm$contrasts
  stattable<-as.data.frame(pairwise_comparisons)
  library(bruceR)
  print_table(stattable, file="BW_baye_high.doc")}

# plot 1
{
  #this is to plot conditional effects
  cond_eff<-conditional_effects(BW_baye, effects = "group", method = "posterior_epred")
  plot_data <- cond_eff$group[, c("group", "estimate__", "lower__", "upper__")]
  plot_save<-ggplot(plot_data, aes(x = group, y = estimate__, color = group)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lower__, ymax = upper__), width = 0.2, size = 2) +
    scale_color_manual(values = c("Sham" = "#1b72b3", "-45°" = "#e1a636", "45°" = "#c73e11")) +
    theme(plot.title=element_text(hjust=0.5,size = 13),legend.position = "none",
          axis.text.x = element_text(size = 13),axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 13),axis.title.y = element_text(size = 15),
          legend.text = element_text(size=13),legend.title = element_text(size=13)) +
    xlab("Stimulation condition") +
    ylab("Fitted bTW (dB)")+
    scale_y_continuous(limits = c(0.4, 0.8))  # Adjust the y-axis limits here
  
  print(plot_save)
 
}


warp_em <- emmeans (BW_baye,  ~ group,epred = TRUE)%>% 
  contrast(method = "pairwise") %>% 
  gather_emmeans_draws()

intere <- warp_em$contrast=='Sham - (-45°)'
p<-mean(warp_em$.value[intere]>0)
p




# Fig.6D: alpha power----

# baye lmm
AP_low = brm(AlphaPower  ~ 1+pre_AP+group+(1+group|id),
                  data=LongData_AP_low,iter = 3000,cores = parallel::detectCores(),control=list(adapt_delta=0.95),seed=42)
AP_high = brm(AlphaPower  ~ 1+pre_AP+group*blockIndex+(1+group|id),
             data=LongData_AP_high,iter = 3000,cores = parallel::detectCores(),control=list(adapt_delta=0.95),seed=42)

AP_ave<-AP_low
AP_ave<-AP_high
summary(AP_ave)
p_direction(AP_ave)
fixef(AP_ave)
emmeans(AP_ave,pairwise~group)

# plot 1
{
  #this is to plot conditional effects
  cond_eff<-conditional_effects(AP_ave, effects = "group", method = "posterior_epred")
  plot_data <- cond_eff$group[, c("group", "estimate__", "lower__", "upper__")]
  plot_save<-ggplot(plot_data, aes(x = group, y = estimate__, color = group)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lower__, ymax = upper__), width = 0.2, size = 2) +
    scale_color_manual(values = c("Sham" = "#1b72b3", "-45°" = "#e1a636", "45°" = "#c73e11")) +
    theme(plot.title=element_text(hjust=0.5,size = 13),legend.position = "none",
          axis.text.x = element_text(size = 13),axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 13),axis.title.y = element_text(size = 15),
          legend.text = element_text(size=13),legend.title = element_text(size=13)) +
    xlab("Stimulation condition") +
    ylab("Fitted AP (dB)")
  print(plot_save)
  
}


warp_em <- emmeans (AP_ave,  ~ group,epred = TRUE)%>% 
  contrast(method = "pairwise") %>% 
  gather_emmeans_draws()

intere <- warp_em$contrast=='Sham - 45°'
p<-mean(warp_em$.value[intere]>0)
p

