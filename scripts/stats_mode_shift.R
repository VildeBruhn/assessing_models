####################################################
##    Assessing models of phenotypic evolution    ##
##                                                ##
##             MODE-SHIFT STATISTICS              ##
##                                                ##
##                Marion Thaureau                 ##
##              Created 01.06.2026                ##
####################################################

# install packages
#install.packages("evoTS") #version 1.0.3
#install.packages("devtools")
#devtools::install_github("klvoje/adePEM") #version 1.1.1
#install.packages("tidyverse")
#install.packages("paleoTS") #version 0.6.1
#install.packages("gridExtra")
#install.packages("lme4")
#install.packages("scales")

# load packages
library(tidyverse)
library(gridExtra)
library(ggplot2)
library(lme4)
library(scales)

rm(list = ls())

# set working directory
PATH = "[PATH_TO_MAIN_FOLDER]"
setwd(PATH)

# import functions
source("./scripts/functions.R")


#--------------------------------
# IMPORT FILES AND PROCESS FILES
#--------------------------------
ln_data_meta_shift <- read_delim("./data/timeseries/metadata.txt", col_names = TRUE, delim = "\t")

# remove time series with less than 14 steps
ln_data_meta_shift <- subset(ln_data_meta_shift, steps >= 14)

# remove modern time series
ln_data_meta_shift <- subset(ln_data_meta_shift, period_start != "Present")

# remove Syverson
ln_data_meta_shift <- subset(ln_data_meta_shift, URL != "https://doi.org/10.1017/pab.2024.37")

# make list based on ID
ln_data_meta_shift <- lapply(split(ln_data_meta_shift,ln_data_meta_shift$tsID), function(x) as.list(x))


# load data from analyses
load("./data/ln_data_shift.Rdata")
load("./data/aicc_shift_passed.Rdata")
load("./data/adeq_shift_passed.Rdata")


# get aicc model info into metadata
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, GRW, "GRW")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, URW, "URW")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, Stasis, "stasis")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, Strict_stasis, "strict stasis")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, Decel, "decel")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, Accel, "accel")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, OU, "OU")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, OU_mov_opt_anc, "OU mov opt anc")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, OU_mov_opt, "OU mov opt")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, Stasis_Stasis, "Stasis-Stasis")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, Stasis_URW, "Stasis-URW")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, Stasis_GRW, "Stasis-GRW")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, Stasis_OU, "Stasis-OU")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, URW_URW, "URW-URW")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, URW_GRW, "URW-GRW")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, URW_OU, "URW-OU")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, GRW_GRW, "GRW-GRW")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, GRW_OU, "GRW-OU")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, OU_OU, "OU-OU")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, OU_GRW, "OU-GRW")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, OU_URW, "OU-URW")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, OU_Stasis, "OU-Stasis")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, GRW_URW, "GRW-URW")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, GRW_Stasis, "GRW-Stasis")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, URW_Stasis, "URW-Stasis")

# get the type of model (without or with shift) info into metadata
ln_data_meta_shift <- model_type(ln_data_meta_shift, GRW, "single-mode")
ln_data_meta_shift <- model_type(ln_data_meta_shift, URW, "single-mode")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Stasis, "single-mode")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Strict_stasis, "single-mode")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Decel, "single-mode")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Accel, "single-mode")
ln_data_meta_shift <- model_type(ln_data_meta_shift, OU, "single-mode")
ln_data_meta_shift <- model_type(ln_data_meta_shift, OU_mov_opt_anc, "single-mode")
ln_data_meta_shift <- model_type(ln_data_meta_shift, OU_mov_opt, "single-mode")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Stasis_Stasis, "mode-shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Stasis_URW, "mode-shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Stasis_GRW, "mode-shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Stasis_OU, "mode-shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, URW_URW, "mode-shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, URW_GRW, "mode-shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, URW_OU, "mode-shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, GRW_GRW, "mode-shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, GRW_OU, "mode-shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, OU_OU, "mode-shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, OU_GRW, "mode-shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, OU_URW, "mode-shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, OU_Stasis, "mode-shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, GRW_URW, "mode-shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, GRW_Stasis, "mode-shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, URW_Stasis, "mode-shift")

# get adequate model info into metadata
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, GRW_adeq_passed, "GRW")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, URW_adeq_passed, "URW")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, Stasis_adeq_passed, "stasis")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, Strict_stasis_adeq_passed, "strict stasis")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, Decel_adeq_passed, "decel")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, Accel_adeq_passed, "accel")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, OU_adeq_passed, "OU")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, OU_mov_opt_anc_adeq_passed, "OU mov opt anc")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, OU_mov_opt_adeq_passed, "OU mov opt")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, Stasis_Stasis_adeq_passed, "Stasis-Stasis")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, Stasis_URW_adeq_passed, "Stasis-URW")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, Stasis_GRW_adeq_passed, "Stasis-GRW")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, Stasis_OU_adeq_passed, "Stasis-OU")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, URW_URW_adeq_passed, "URW-URW")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, URW_GRW_adeq_passed, "URW-GRW")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, URW_OU_adeq_passed, "URW-OU")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, GRW_GRW_adeq_passed, "GRW-GRW")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, GRW_OU_adeq_passed, "GRW-OU")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, OU_OU_adeq_passed, "OU-OU")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, OU_GRW_adeq_passed, "OU-GRW")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, OU_URW_adeq_passed, "OU-URW")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, OU_Stasis_adeq_passed, "OU-Stasis")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, GRW_URW_adeq_passed, "GRW-URW")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, GRW_Stasis_adeq_passed, "GRW-Stasis")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, URW_Stasis_adeq_passed, "URW-Stasis")

# bind data to dataframe
unit_list <- c("popID", "total_N", "steps", "interval_MY", "trait_type",
               "lat", "lon", "model_aicc", "model_type", "model_adequate")
plot_data <- bind(ln_data_meta_shift, unit_list)

# remove time series with NA
plot_data <- plot_data %>% drop_na(model_aicc)

# ordering the model according to the number of parameters
level_order <- c("single-mode", "mode-shift")
plot_data$model_type <- factor(plot_data$model_type, levels = level_order)
plot_data$model_type <- relevel(plot_data$model_type, ref = "single-mode")


#--------------------
# STATS RELATIVE FIT
#--------------------

###### interval MY ######

# LMM
lmm_model_time <- lmer(interval_MY ~ model_type + (1| popID), data = plot_data)
summary(lmm_model_time)

# plot
#pdf("[PATH_TO_OUTPUT_FOLDER]/interval_my_shift_aicc.pdf")
intv_plot <- ggplot(plot_data, aes(x = interval_MY, y = factor(model_type, levels = level_order), fill = model_type
                                 )) + 
  geom_boxplot(alpha = 0.7) + theme_classic() + scale_y_discrete(labels = c("single-mode", "mode-shift")) + 
  scale_fill_manual(values = c("mode-shift" = "#B84400", "single-mode" = "#B7AF3B"), name = "Model type") +
  xlab("Interval (MY)") + ylab("Model type") + theme(axis.text = element_text(size = 14),
                                                legend.position = "none",
                                                axis.title = element_text(size = 18),
                                                axis.title.x = element_text(margin = margin(t = 18, r = 0, b = 0, l = 0)))
intv_plot
#dev.off()


###### steps ######

# LMM
lmm_model_steps <- lmer(steps ~ model_type + (1| popID), data = plot_data)
summary(lmm_model_steps)

# plot
#pdf("[PATH_TO_OUTPUT_FOLDER]/steps_shift_aicc.pdf")
steps_plot <- ggplot(plot_data, aes(x = log(steps), y = factor(model_type, levels = level_order), fill = model_type
)) + 
  geom_boxplot(alpha = 0.7) + theme_classic() + scale_y_discrete(labels = c("single-mode", "mode-shift")) + 
  scale_fill_manual(values = c("mode-shift" = "#B84400", "single-mode" = "#B7AF3B"), name = "Model type") + 
  xlab("ln(Steps)") + ylab("") + theme(axis.text = element_text(size = 14),
                                                     legend.position = "none",
                                                     axis.title = element_text(size = 18),
                                                     axis.title.y = element_blank(),
                                                     axis.title.x = element_text(margin = margin(t = 18, r = 0, b = 0, l = 0)))
steps_plot
#dev.off()


###### resolution ######

# LMM
plot_data$resolution = plot_data$steps/plot_data$interval_MY
lmm_model_resolution <- lmer(resolution ~ model_type + (1| popID), data = plot_data)
summary(lmm_model_resolution)

# plot
#pdf("[PATH_TO_OUTPUT_FOLDER]/resolution_shift_aicc.pdf")
res_plot <- ggplot(plot_data, aes(x = log(resolution), y = factor(model_type, levels = level_order), fill = model_type
)) + 
  geom_boxplot(alpha = 0.7) + theme_classic() + scale_y_discrete(labels = c("single-mode", "mode-shift")) + 
  scale_fill_manual(values = c("mode-shift" = "#B84400", "single-mode" = "#B7AF3B"), name = "Model type") + 
  xlab("ln(Resolution)") + ylab("") + theme(axis.text = element_text(size = 14),
                                                 legend.position = "none",
                                                 axis.title = element_text(size = 18),
                                                 axis.title.y = element_blank(),
                                                 axis.title.x = element_text(margin = margin(t = 18, r = 0, b = 0, l = 0)))
res_plot
#dev.off()


# put plots in same figure
#pdf(width = 15.0, height = 5.5, file = "[PATH_TO_OUTPUT_FOLDER]/all_shift_aicc.pdf")
grid.arrange(intv_plot,
             steps_plot,
             res_plot, nrow = 1)
#dev.off()



#-------------------
# PLOT ABSOLUTE FIT
#-------------------


# remove NA data
plot_data2 <- plot_data 
plot_data2 <- plot_data2 %>% drop_na(model_adequate)

# ordering the model according to the number of parameters
level_order <- c("single-mode", "mode-shift")
plot_data2$model_type <- factor(plot_data2$model_type, levels = level_order)
plot_data2$model_type <- relevel(plot_data2$model_type, ref = "single-mode")


###### interval MY ######

# LMM
lmm_model_time2 <- lmer(interval_MY ~ model_type + (1| popID), data = plot_data2)
summary(lmm_model_time2)

# plot
#pdf("[PATH_TO_OUTPUT_FOLDER]/interval_my_shift_adeq.pdf")
intv_plot <- ggplot(plot_data2, aes(x = interval_MY, y = factor(model_type, levels = level_order), fill = model_type
)) + 
  geom_boxplot(alpha = 0.7) + theme_classic() + scale_y_discrete(labels = c("single-mode", "mode-shift")) + 
  scale_fill_manual(values = c("mode-shift" = "#B84400", "single-mode" = "#B7AF3B"), name = "Model type") + 
  xlab("Interval (MY)") + ylab("Model type") + theme(axis.text = element_text(size = 14),
                                                     legend.position = "none",
                                                     axis.title = element_text(size = 18),
                                                     axis.title.x = element_text(margin = margin(t = 18, r = 0, b = 0, l = 0)))
intv_plot
#dev.off()

###### steps ######

# LMM
lmm_model_steps2 <- lmer(steps ~ model_type + (1| popID), data = plot_data2)
summary(lmm_model_steps2)

# plot
#pdf("[PATH_TO_OUTPUT_FOLDER]/steps_shift_adeq.pdf")
steps_plot <- ggplot(plot_data2, aes(x = log(steps), y = factor(model_type, levels = level_order), fill = model_type
)) + 
  geom_boxplot(alpha = 0.7) + theme_classic() + scale_y_discrete(labels = c("single-mode", "mode-shift")) + 
  scale_fill_manual(values = c("mode-shift" = "#B84400", "single-mode" = "#B7AF3B"), name = "Model type") + 
  xlab("ln(Steps)") + ylab("") + theme(axis.text = element_text(size = 14),
                                                 legend.position = "none",
                                                 axis.title = element_text(size = 18),
                                                 axis.title.y = element_blank(),
                                                 axis.title.x = element_text(margin = margin(t = 18, r = 0, b = 0, l = 0)))
steps_plot
#dev.off()

###### resolution ######

# LMM
plot_data2$resolution = plot_data2$steps/plot_data2$interval_MY
lmm_model_resolution2 <- lmer(resolution ~ model_type + (1| popID), data = plot_data2)
summary(lmm_model_resolution2)

# plot
#pdf("[PATH_TO_OUTPUT_FOLDER]/resolution_shift_adeq.pdf")
res_plot <- ggplot(plot_data2, aes(x = log(resolution), y = factor(model_type, levels = level_order), fill = model_type
)) + 
  geom_boxplot(alpha = 0.7) + theme_classic() + scale_y_discrete(labels = c("single-mode", "mode-shift")) + 
  scale_fill_manual(values = c("mode-shift" = "#B84400", "single-mode" = "#B7AF3B"), name = "Model type") + 
  xlab("ln(Resolution)") + ylab("") + theme(axis.text = element_text(size = 14),
                                                      legend.position = "none",
                                                      axis.title = element_text(size = 18),
                                                      axis.title.y = element_blank(),
                                                      axis.title.x = element_text(margin = margin(t = 18, r = 0, b = 0, l = 0)))
res_plot
#dev.off()



# put plots in same figure
#pdf(width = 15.0, height = 5.5, file = "[PATH_TO_OUTPUT_FOLDER]/all_shift_adeq.pdf")
grid.arrange(intv_plot,
             steps_plot,
             res_plot, nrow = 1)
#dev.off()
