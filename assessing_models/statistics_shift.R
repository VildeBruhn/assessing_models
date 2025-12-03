##############################################################
##   Assessing models of evolution (mode shift included)    ##
##                       STATISTICS                         ##
##############################################################

# install packages
#install.packages("evoTS") #version 1.0.3
#install.packages("devtools")
#devtools::install_github("klvoje/adePEM")
#install.packages("tidyverse")
#install.packages("paleoTS") #version 0.6.1
#install.packages("wesanderson")

# load packages
library(tidyverse)
library(devtools)
library(wesanderson)
library(gridExtra)
library(ggplot2)
library(lme4)
library(lmerTest)
library(dplyr)

rm(list = ls())

# set working directory
setwd("/Users/vildebruhnkinneberg/Documents/GitHub/assessing_models_evolution/assessing_models/")
source("./assessing_models_uni_functions.R")


#--------------
# IMPORT FILES
#--------------


# load data
load("./ln_data_meta_shift.Rdata")
load("./ln_data_shift.Rdata")

# load relative fit time series
load("./aicc_shift_passed.Rdata")

# load adequate time series 
load("./adeq_shift_passed.Rdata")


#--------------
# PREPARE DATA
#--------------


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
ln_data_meta_shift <- model_type(ln_data_meta_shift, GRW, "no shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, URW, "no shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Stasis, "no shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Strict_stasis, "no shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Decel, "no shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Accel, "no shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, OU, "no shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, OU_mov_opt_anc, "no shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, OU_mov_opt, "no shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Stasis_Stasis, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Stasis_URW, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Stasis_GRW, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Stasis_OU, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, URW_URW, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, URW_GRW, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, URW_OU, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, GRW_GRW, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, GRW_OU, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, OU_OU, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, OU_GRW, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, OU_URW, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, OU_Stasis, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, GRW_URW, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, GRW_Stasis, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, URW_Stasis, "mode shift")

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
level_order <- c("no shift", "mode shift")
plot_data$model_type <- factor(plot_data$model_type, levels = level_order)
plot_data$model_type <- relevel(plot_data$model_type, ref = "no shift")


#--------------------
# STATS RELATIVE FIT
#--------------------

###### interval MY ######
intv_my <- plot_data[c("model_type", "interval_MY")]

# LMM
#lmm_model_time <- lmer(interval_MY ~ 1+ model_type + (1| popID), data = plot_data)
lmm_model_time <- lm(interval_MY ~ model_type, data = plot_data)
summary(lmm_model_time)

# plot
pdf("./results_paleoTS_v0.6.1/plot/interval_my_uni_aicc.pdf")
intv_plot <- ggplot(intv_my, aes(x = interval_MY, y = factor(model_aicc, levels = level_order),
                                 fill = parameters)) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("Stasis", "URW", 
                                                                 "GRW", "Accel.", "Decel.", 
                                                                 "OU", "OU mov. opt.")) + 
  scale_fill_discrete(name = "Parameters", labels = c("2", "3", "4",
                                                      "5"), palette = col_val2) +
  xlab("Interval (MY)") + ylab("Model") + theme(axis.text = element_text(size = 10),
                                                axis.title = element_text(size = 13),
                                                axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)))
intv_plot
dev.off()

###### steps ######
steps <- plot_data[c("model_type", "steps")]

# LMM
#lmm_model_steps <- lmer(steps ~ model_aicc + (1| popID), data = plot_data)
lmm_model_steps <- lm(steps ~ model_type, data = plot_data)
summary(lmm_model_steps)

# plot
pdf("./results_paleoTS_v0.6.1/plot/steps_uni_aicc.pdf")
steps_plot <- ggplot(steps, aes(x = log(steps), y = factor(model_aicc, levels = level_order), fill = parameters)) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("Stasis", "URW", 
                                                                 "GRW", "Accel.", "Decel.", 
                                                                 "OU", "OU mov. opt.")) +
  scale_fill_discrete(name = "Parameters", labels = c("2", "3", "4",
                                                      "5"), palette = col_val2) +
  xlab("ln(Steps)") + ylab("Model") + theme(axis.text = element_text(size = 10),
                                            axis.title = element_text(size = 13),
                                            axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)))
steps_plot
dev.off()

###### resolution ######
res <- plot_data[c("model_aicc", "steps", "interval_MY")]
res$resolution <- res$steps/res$interval_MY

# LMM
plot_data$resolution = plot_data$steps/plot_data$interval_MY
#lmm_model_resolution <- lmer(resolution ~ 1+ model_aicc + (1| popID), data = plot_data)
lmm_model_resolution <- lm(resolution ~ model_type, data = plot_data)
summary(lmm_model_resolution)

# plot
pdf("./results_paleoTS_v0.6.1/plot/resolution_uni_aicc.pdf")
res_plot <- ggplot(res, aes(x = log(resolution), y = factor(model_aicc, levels = level_order), fill = parameters)) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("Stasis", "URW", 
                                                                 "GRW", "Accel.", "Decel.", 
                                                                 "OU", "OU mov. opt.")) +
  scale_fill_discrete(name = "Parameters", labels = c("2", "3", "4",
                                                      "5"), palette = col_val2) +
  xlab("ln(Resolution)") + ylab("Model") + theme(axis.text = element_text(size = 10),
                                                 axis.title = element_text(size = 13),
                                                 axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)))
res_plot
dev.off()



# put plots in same figure
pdf(width = 15.0, height = 5.5, file = "[PATH_TO_RESULTS_FOLDER]/empirical.pdf")
grid.arrange(intv_plot,
             steps_plot,
             res_plot, nrow = 1)
dev.off()



# put the output of the linear mixed effects models in the same table
models <- rownames(summary(lmm_model_time)$coefficients)
interval <- as.data.frame(summary(lmm_model_time)$coefficients)
steps <- as.data.frame(summary(lmm_model_steps)$coefficients)
resolution <- as.data.frame(summary(lmm_model_resolution)$coefficients)

interval_df <- data.frame(term = models, i.Est = interval$Estimate,
                          i.SE = interval$"Std", i.p.value = interval$"Pr")
interval_df$i.Est[2:nrow(interval_df)] <- interval_df$i.Est[1] + interval_df$i.Est[2:nrow(interval_df)]

steps_df <- data.frame(s.Est = steps$Estimate,
                       s.SE = steps$"Std", s.p.value = steps$"Pr")
steps_df$s.Est[2:nrow(steps_df)] <- steps_df$s.Est[1] + steps_df$s.Est[2:nrow(steps_df)]

resolution_df <- data.frame(r.Est = resolution$Estimate,
                            r.SE = resolution$"Std", r.p.value = resolution$"Pr")
resolution_df$r.Est[2:nrow(resolution_df)] <- resolution_df$r.Est[1] + resolution_df$r.Est[2:nrow(resolution_df)]

lmm_result_table <- cbind(interval_df, steps_df, resolution_df)
lmm_result_table[, -1] <- round(lmm_result_table[, -1], 3)

write.csv(lmm_result_table, file = "./results_paleoTS_v0.6.1/table_lmm_shift.pdf", row.names = FALSE)



#-------------------
# PLOT ABSOLUTE FIT
#-------------------


# remove NA data
plot_data2 <- plot_data 
plot_data2 <- plot_data2 %>% drop_na(model_adequate)

# ordering the model according to the number of parameters
level_order <- c("no shift", "mode shift")
plot_data$model_type <- factor(plot_data$model_type, levels = level_order)
plot_data$model_type <- relevel(plot_data$model_type, ref = "no shift")


###### interval MY ######
intv_my <- plot_data2[c("model_adequate", "interval_MY", "parameters")]

# LMM
#lmm_model_time2 <- lmer(interval_MY ~ 1+ model_type + (1| popID), data = plot_data2)
lmm_model_time2 <- lm(interval_MY ~ model_type, data = plot_data2)
summary(lmm_model_time2)

# plot
level_order <- c("stasis", "URW", "accel", "decel", "OU", "OU mov opt")
pdf("./results_paleoTS_v0.6.1/plot/interval_my_uni_adeq.pdf")
intv_plot <- ggplot(intv_my, aes(x = interval_MY, y = factor(model_adequate, levels = level_order),
                                 fill = parameters)) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("Stasis", "URW", 
                                                                 "Accel.", "Decel.", 
                                                                 "OU", "OU mov. opt.")) + 
  scale_fill_discrete(name = "Parameters", labels = c("2", "3", "4",
                                                      "5"), palette = col_val2) +
  xlab("Interval (MY)") + ylab("Model") + theme(axis.text = element_text(size = 10),
                                                axis.title = element_text(size = 13),
                                                axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)))
intv_plot
dev.off()

###### steps ######
steps <- plot_data2[c("model_adequate", "steps", "parameters")]

# LMM
#lmm_model_steps2 <- lmer(steps ~ model_type + (1| popID), data = plot_data2)
lmm_model_steps2 <- lm(steps ~ model_type, data = plot_data2)
summary(lmm_model_steps2)

# plot
level_order <- c("stasis", "URW", "accel", "decel", "OU", "OU mov opt")
pdf("./results_paleoTS_v0.6.1/plot/steps_uni_adeq.pdf")
steps_plot <- ggplot(steps, aes(x = log(steps), y = factor(model_adequate, levels = level_order), fill = parameters)) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("Stasis", "URW", 
                                                                 "Accel.", "Decel.", 
                                                                 "OU", "OU mov. opt.")) +
  scale_fill_discrete(name = "Parameters", labels = c("2", "3", "4",
                                                      "5"), palette = col_val2) +
  xlab("ln(Steps)") + ylab("Model") + theme(axis.text = element_text(size = 10),
                                            axis.title = element_text(size = 13),
                                            axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)))
steps_plot
dev.off()

###### resolution ######
res <- plot_data2[c("model_adequate", "steps", "interval_MY", "parameters")]
res$resolution <- res$steps/res$interval_MY

# LMM
plot_data2$resolution = plot_data2$steps/plot_data2$interval_MY
#lmm_model_resolution2 <- lmer(resolution ~ 1+ model_type + (1| popID), data = plot_data2)
lmm_model_resolution2 <- lm(resolution ~ model_type, data = plot_data2)
summary(lmm_model_resolution2)

# plot
level_order <- c("stasis", "URW", "accel", "decel", "OU", "OU mov opt")
pdf("./results_paleoTS_v0.6.1/plot/resolution_uni_adeq.pdf")
res_plot <- ggplot(res, aes(x = log(resolution), y = factor(model_adequate, levels = level_order), fill = parameters)) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("Stasis", "URW", 
                                                                 "Accel.", "Decel.", 
                                                                 "OU", "OU mov. opt.")) +
  scale_fill_discrete(name = "Parameters", labels = c("2", "3", "4",
                                                      "5"), palette = col_val2) +
  xlab("ln(Resolution)") + ylab("Model") + theme(axis.text = element_text(size = 10),
                                                 axis.title = element_text(size = 13),
                                                 axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)))
res_plot
dev.off()



# put plots in same figure
pdf(width = 15.0, height = 5.5, file = "./results_paleoTS_v0.6.1/plot/interval_steps_res_adeq.pdf")
grid.arrange(intv_plot,
             steps_plot,
             res_plot, nrow = 1)
dev.off()

# put the output of the linear mixed effects models in the same table
models2 = rownames(summary(lmm_model_time2)$coefficients)
interval2 = as.data.frame(summary(lmm_model_time2)$coefficients)
steps2 = as.data.frame(summary(lmm_model_steps2)$coefficients)
resolution2 = as.data.frame(summary(lmm_model_resolution2)$coefficients)

interval_df2 <- data.frame(term = models2, i.Est = interval2$Estimate,
                           i.SE = interval2$"Std", i.p.value = interval2$"Pr")

steps_df2 <- data.frame(s.Est = steps2$Estimate,
                        s.SE = steps2$"Std", s.p.value = steps2$"Pr")

resolution_df2 <- data.frame(r.Est = resolution2$Estimate,
                             r.SE = resolution2$"Std", r.p.value = resolution2$"Pr")

lmm_result_table2 <- cbind(interval_df2, steps_df2, resolution_df2)

write.csv(lmm_result_table2, file = "./results_paleoTS_v0.6.1/table_lmm_shift_adeq.pdf", row.names = FALSE)


