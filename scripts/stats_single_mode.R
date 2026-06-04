####################################################
##    Assessing models of phenotypic evolution    ##
##                                                ##
##           SINGLE-MODE STATISTICS               ##
##                                                ##
##            Vilde Bruhn Kinneberg               ##
##            Created 01.06.2026                  ##
####################################################

# install packages
#install.packages("tidyverse")
#install.packages("gridExtra")
#install.packages("lme4")

# load packages
library(tidyverse)
library(gridExtra)
library(lme4)

rm(list = ls())

# set working directory
PATH = "[PATH_TO_MAIN_FOLDER]"
setwd(PATH)

# import functions
source("./scripts/functions.R")


#-------------------------
# IMPORT AND PROCESS FILES
#-------------------------


df <- read_delim("./data/timeseries/metadata.txt", col_names = TRUE, delim = "\t")

# remove time series with less than 7 steps
df <- subset(df, steps >= 7)

# remove modern time series
df <- subset(df, period_start != "Present")

# remove Syverson
df <- subset(df, URL != "https://doi.org/10.1017/pab.2024.37")

# make list based on ID
df <- lapply(split(df,df$tsID), function(x) as.list(x))

# load data from analyses
load("./data/aicc_single_mode_passed.Rdata")
load("./data/model_test_single_mode.Rdata")
load("./data/adeq_single_mode_passed.Rdata")


#--------------
# PREPARE DATA
#--------------


# get aicc model info into metadata
df <- model_aicc(df, GRW, "GRW")
df <- model_aicc(df, URW, "URW")
df <- model_aicc(df, stasis, "stasis")
df <- model_aicc(df, strict_stasis, "strict stasis")
df <- model_aicc(df, decel, "decel")
df <- model_aicc(df, accel, "accel")
df <- model_aicc(df, OU, "OU")
df <- model_aicc(df, OU_mov_opt_anc, "OU mov opt anc")
df <- model_aicc(df, OU_mov_opt, "OU mov opt")

# get adequate model info into metadata
df <- model_adeq(df, GRW_adeq_passed, "GRW")
df <- model_adeq(df, URW_adeq_passed, "URW")
df <- model_adeq(df, stasis_adeq_passed, "stasis")
df <- model_adeq(df, strict_stasis_adeq_passed, "strict stasis")
df <- model_adeq(df, decel_adeq_passed, "decel")
df <- model_adeq(df, accel_adeq_passed, "accel")
df <- model_adeq(df, OU_adeq_passed, "OU")
df <- model_adeq(df, OU_mov_opt_anc_adeq_passed, "OU mov opt anc")
df <- model_adeq(df, OU_mov_opt_adeq_passed, "OU mov opt")

# bind data to dataframe
unit_list <- c("popID","total_N", "steps", "interval_MY", "model_aicc", "model_adequate")
plot_data <- bind(df, unit_list)

# remove time series with NA
plot_data <- plot_data %>% drop_na(model_aicc)

# collapse strict stasis to stasis
plot_data$model_aicc <- replace(plot_data$model_aicc, plot_data$model_aicc == "strict stasis", "stasis")
plot_data$model_adequate <- replace(plot_data$model_adequate, plot_data$model_adequate == "strict stasis", "stasis")

# collapse OU mov opt anc to OU mov opt
plot_data$model_aicc <- replace(plot_data$model_aicc, plot_data$model_aicc == "OU mov opt anc", "OU mov opt")
plot_data$model_adequate <- replace(plot_data$model_adequate, plot_data$model_adequate == "OU mov opt anc", "OU mov opt")

# add parameters column
plot_data$parameters <- plot_data$model_aicc
plot_data$parameters <- replace(plot_data$parameters, plot_data$parameters == "stasis", 2)
plot_data$parameters <- replace(plot_data$parameters, plot_data$parameters == "URW", 2)
plot_data$parameters <- replace(plot_data$parameters, plot_data$parameters == "GRW", 3)
plot_data$parameters <- replace(plot_data$parameters, plot_data$parameters == "accel", 3)
plot_data$parameters <- replace(plot_data$parameters, plot_data$parameters == "decel", 3)
plot_data$parameters <- replace(plot_data$parameters, plot_data$parameters == "OU", 4)
plot_data$parameters <- replace(plot_data$parameters, plot_data$parameters == "OU mov opt", 5)

# ordering the model according to the number of parameters
level_order <- c("stasis", "URW", "GRW", "accel", "decel", "OU", "OU mov opt")
plot_data$model_aicc <- factor(plot_data$model_aicc, levels = level_order)

# set colors
col_val <- c("#B7AF3B", "#BFE0E1", "#EADC72", "#B84400")


#---------------
# RELATIVE FIT
#---------------


###### interval MY ######
intv <- plot_data

# LMM
lmm_intv <- lmer(interval_MY ~ 1 + model_aicc + (1| popID), data = intv)
summary(lmm_intv)

# plot
#pdf("[PATH_TO_OUTPUT_FOLDER]/interval_single_mode_aicc.pdf")
intv_plot <- ggplot(intv, aes(x = interval_MY, y = model_aicc,
                    fill = parameters)) +
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("Stasis", "URW", 
                                                                 "GRW", "Accel.", "Decel.", 
                                                                 "OU", "OU mov. opt.")) + 
  scale_fill_discrete(name = "Parameters", labels = c("2", "3", "4",
                                                      "5"), palette = col_val) +
  xlab("Interval (MY)") + ylab("Model") + theme(axis.text = element_text(size = 14),
                                                axis.title = element_text(size = 18),
                                                axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)),
                                                legend.title = element_text(size = 14), legend.text = element_text(size = 14))
intv_plot
#dev.off()


###### steps ######
steps <- plot_data

# LMM
lmm_steps <- lmer(steps ~ model_aicc + (1| popID), data = steps)
summary(lmm_steps)

# plot
#pdf("[PATH_TO_OUTPUT_FOLDER]/steps_single_mode_aicc.pdf")
steps_plot <- ggplot(steps, aes(x = log(steps), y = model_aicc, fill = parameters)) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("Stasis", "URW", 
                                                                 "GRW", "Accel.", "Decel.", 
                                                                 "OU", "OU mov. opt.")) +
  scale_fill_discrete(name = "Parameters", labels = c("2", "3", "4",
                                                      "5"), palette = col_val) +
  xlab("ln(Steps)") + ylab("Model") + theme(axis.text = element_text(size = 14),
                                        axis.title = element_text(size = 18),
                                        axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)),
                                        legend.title = element_text(size = 14), legend.text = element_text(size = 14))
steps_plot
#dev.off()


###### resolution ######
res <- plot_data
res$resolution <- res$steps/res$interval_MY

# LMM
lmm_res <- lmer(resolution ~ 1 + model_aicc + (1| popID), data = res)
summary(lmm_res)

# LMM remove outliers GRW
res_no <- res[-c(519,520),]
lmm_res_no <- lmer(resolution ~ 1 + model_aicc + (1| popID), data = res_no)
summary(lmm_res_no)

# plot
#pdf("[PATH_TO_OUTPUT_FOLDER]/resolution_single_mode_aicc.pdf")
res_plot <- ggplot(res_no, aes(x = log(resolution), y = model_aicc, fill = parameters)) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("Stasis", "URW", 
                                                                 "GRW", "Accel.", "Decel.", 
                                                                 "OU", "OU mov. opt.")) +
  scale_fill_discrete(name = "Parameters", labels = c("2", "3", "4",
                                                      "5"), palette = col_val) +
  xlab("ln(Resolution)") + ylab("Model") + theme(axis.text = element_text(size = 14),
                                                        axis.title = element_text(size = 18),
                                                        axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)),
                                                        legend.title = element_text(size = 14), legend.text = element_text(size = 14))
res_plot
#dev.off()



# put plots in same figure
#pdf(width = 15.0, height = 5.5, file = "[PATH_TO_OUTPUT_FOLDER]/all_single_mode_aicc.pdf")
grid.arrange(intv_plot,
             steps_plot,
             res_plot, nrow = 1)
#dev.off()


#---------------
# ABSOLUTE FIT
#---------------


# remove NA data
plot_data2 <- plot_data 
plot_data2 <- plot_data2 %>% drop_na(model_adequate)

# ordering the model according to the number of parameters
plot_data2$model_adequate <- factor(plot_data2$model_adequate, levels = level_order)


###### interval MY ######
intv2 <- plot_data2

# LMM
lmm_intv2 <- lmer(interval_MY ~ 1+ model_adequate + (1| popID), data = intv2)
summary(lmm_intv2)

# plot
#pdf("[PATH_TO_OUTPUT_FOLDER]/interval_single_mode_adeq.pdf")
intv_plot2 <- ggplot(intv2, aes(x = interval_MY, y = model_adequate,
                                 fill = parameters)) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("Stasis", "URW", 
                                                                 "GRW", "Accel.", "Decel.", 
                                                                 "OU", "OU mov. opt.")) + 
  scale_fill_discrete(name = "Parameters", labels = c("2", "3", "4",
                                                      "5"), palette = col_val) +
  xlab("Interval (MY)") + ylab("Model") + theme(axis.text = element_text(size = 14),
                                                axis.title = element_text(size = 18),
                                                axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)),
                                                legend.title = element_text(size = 14), legend.text = element_text(size = 14))
intv_plot2
#dev.off()


###### steps ######
steps2 <- plot_data2

# LMM
lmm_steps2 <- lmer(steps ~ model_adequate + (1| popID), data = steps2)
summary(lmm_steps2)

# plot
#pdf("[PATH_TO_OUTPUT_FOLDER]/steps_single_mode_adeq.pdf")
steps_plot2 <- ggplot(steps2, aes(x = log(steps), y = model_adequate, fill = parameters)) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("Stasis", "URW", 
                                                                 "GRW", "Accel.", "Decel.", 
                                                                 "OU", "OU mov. opt.")) +
  scale_fill_discrete(name = "Parameters", labels = c("2", "3", "4",
                                                      "5"), palette = col_val) +
  xlab("ln(Steps)") + ylab("Model") + theme(axis.text = element_text(size = 14),
                                            axis.title = element_text(size = 18),
                                            axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)),
                                            legend.title = element_text(size = 14), legend.text = element_text(size = 14))
steps_plot2
#dev.off()


###### resolution ######
res2 <- plot_data2
res2$resolution <- res2$steps/res2$interval_MY

# LMM
lmm_res2 <- lmer(resolution ~ 1+ model_aicc + (1| popID), data = res2)
summary(lmm_res2)

# LMM remove outliers GRW
res2_no <- res2[-c(414,415),]
lmm_res2_no <- lmer(resolution ~ 1 + model_aicc + (1| popID), data = res2_no)
summary(lmm_res2_no)

# plot
#pdf("[PATH_TO_OUTPUT_FOLDER]/resolution_single_mode_adeq.pdf")
res_plot2 <- ggplot(res2, aes(x = log(resolution), y = model_adequate, fill = parameters)) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("Stasis", "URW", 
                                                                 "GRW", "Accel.", "Decel.", 
                                                                 "OU", "OU mov. opt.")) +
  scale_fill_discrete(name = "Parameters", labels = c("2", "3", "4",
                                                      "5"), palette = col_val) +
  xlab("ln(Resolution)") + ylab("Model") + theme(axis.text = element_text(size = 14),
                                                 axis.title = element_text(size = 18),
                                                 axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)),
                                                 legend.title = element_text(size = 14), legend.text = element_text(size = 14))
                                                 
res_plot2
#dev.off()


# put plots in same figure
#pdf(width = 15.0, height = 5.5, file = "[PATH_TO_OUTPUT_FOLDER]/all_single_mode_adeq.pdf")
grid.arrange(intv_plot2,
             steps_plot2,
             res_plot2, nrow = 1)
#dev.off()
