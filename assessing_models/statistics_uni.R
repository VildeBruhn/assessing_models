#########################################
##   Assessing models of evolution     ##
##            STATISTICS               ##
#########################################

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
load("./ln_data_meta_uni.Rdata")
load("./ln_data_uni.Rdata")

# load relative fit time series
load("./aicc_uni_passed.Rdata")

# load adequate time series 
load("./adeq_uni_passed.Rdata")


#--------------
# PREPARE DATA
#--------------


# get aicc model info into metadata
ln_data_meta <- model_aicc(ln_data_meta, GRW, "GRW")
ln_data_meta <- model_aicc(ln_data_meta, URW, "URW")
ln_data_meta <- model_aicc(ln_data_meta, stasis, "stasis")
ln_data_meta <- model_aicc(ln_data_meta, strict_stasis, "strict stasis")
ln_data_meta <- model_aicc(ln_data_meta, decel, "decel")
ln_data_meta <- model_aicc(ln_data_meta, accel, "accel")
ln_data_meta <- model_aicc(ln_data_meta, OU, "OU")
ln_data_meta <- model_aicc(ln_data_meta, OU_mov_opt_anc, "OU mov opt anc")
ln_data_meta <- model_aicc(ln_data_meta, OU_mov_opt, "OU mov opt")

# get adequate model info into metadata
ln_data_meta <- model_adeq(ln_data_meta, GRW_adeq_passed, "GRW")
ln_data_meta <- model_adeq(ln_data_meta, URW_adeq_passed, "URW")
ln_data_meta <- model_adeq(ln_data_meta, stasis_adeq_passed, "stasis")
ln_data_meta <- model_adeq(ln_data_meta, strict_stasis_adeq_passed, "strict stasis")
ln_data_meta <- model_adeq(ln_data_meta, decel_adeq_passed, "decel")
ln_data_meta <- model_adeq(ln_data_meta, accel_adeq_passed, "accel")
ln_data_meta <- model_adeq(ln_data_meta, OU_adeq_passed, "OU")
ln_data_meta <- model_adeq(ln_data_meta, OU_mov_opt_anc_adeq_passed, "OU mov opt anc")
ln_data_meta <- model_adeq(ln_data_meta, OU_mov_opt_adeq_passed, "OU mov opt")

# bind data to dataframe
unit_list <- c("popID", "total_N", "steps", "interval_MY", "trait_type", "microfossil",
               "lat", "lon", "sediment", "model_aicc", "model_adequate", 
               "environment")
plot_data <- bind(ln_data_meta, unit_list)

# collapse strict stasis to stasis
plot_data$model_aicc <- replace(plot_data$model_aicc, plot_data$model_aicc == "strict stasis", "stasis")
plot_data$model_adequate <- replace(plot_data$model_adequate, plot_data$model_adequate == "strict stasis", "stasis")

# collapse OU mov opt anc to OU mov opt
plot_data$model_aicc <- replace(plot_data$model_aicc, plot_data$model_aicc == "OU mov opt anc", "OU mov opt")
plot_data$model_adequate <- replace(plot_data$model_adequate, plot_data$model_adequate == "OU mov opt anc", "OU mov opt")

# collapse all marine into one
#plot_data$environment <- replace(plot_data$environment, plot_data$environment == "marine benthic", "marine")
#plot_data$environment <- replace(plot_data$environment, plot_data$environment == "marine pelagic", "marine")
#plot_data$environment <- replace(plot_data$environment, plot_data$environment == "marine planktic", "marine")

# add parameters column
plot_data$parameters <- plot_data$model_aicc
plot_data$parameters <- replace(plot_data$parameters, plot_data$parameters == "stasis", 2)
plot_data$parameters <- replace(plot_data$parameters, plot_data$parameters == "URW", 2)
plot_data$parameters <- replace(plot_data$parameters, plot_data$parameters == "GRW", 3)
plot_data$parameters <- replace(plot_data$parameters, plot_data$parameters == "accel", 3)
plot_data$parameters <- replace(plot_data$parameters, plot_data$parameters == "decel", 3)
plot_data$parameters <- replace(plot_data$parameters, plot_data$parameters == "OU", 4)
plot_data$parameters <- replace(plot_data$parameters, plot_data$parameters == "OU mov opt", 5)

# remove time series with NA
plot_data <- plot_data %>% drop_na(model_aicc)
#plot_data <- plot_data %>% drop_na(environment)
<<<<<<< HEAD

# ordering the model according to the number of parameters
level_order <- c("stasis", "URW", "GRW", "accel", "decel", "OU", "OU mov opt")
plot_data$model_aicc <- factor(plot_data$model_aicc, levels = level_order)
plot_data$model_aicc <- relevel(plot_data$model_aicc, ref = "stasis")
=======
>>>>>>> 15cd3634ef9dd1df70f22064c075a40aa283aa18

# set colors
col_val1 <- c("#85B7B9", "#DCCB4E")
col_val2 <- c("#F11B00","#E5A208", "#ADC397", "#3A9AB2")
col_val_extra <- c("#3A9AB2", "#85B7B9", "#ADC397", "#DCCB4E", "#E5A208", "#ED6E04", "#F11B00")


#-------------------
# PLOT RELATIVE FIT
#-------------------

###### micro vs. macro ######
#micro_macro <- plot_data
#level_order <- c("stasis", "URW", "GRW", "accel", "decel", "OU", "OU mov opt")
#pdf("./results_paleoTS_v0.6.1/plot/micro_macro_aicc.pdf")
#ggplot(micro_macro, aes(x = factor(model_aicc, levels = level_order), fill = microfossil)) + geom_bar() +
#  theme_classic() + scale_x_discrete(labels = c("Stasis", "URW", "GRW", "Accel.", "Decel.", "OU", "OU mov. opt.")) +
#  labs(fill = "") + scale_fill_discrete(name = "", labels = c("Macrofossils", "Microfossils"), palette = col_val1) +
#  xlab("Model") + ylab("Count") + theme(legend.text = element_text(size = 15))
#dev.off()

<<<<<<< HEAD
###### environment ######
#env <- plot_data
#level_order <- c("stasis", "URW", "GRW", "accel", "decel", "OU", "OU mov opt")
#pdf("./results_paleoTS_v0.6.1/plot/environment_aicc.pdf")
#ggplot(env, aes(x = factor(model_aicc, levels = level_order), fill = environment)) + geom_bar() +
#  theme_classic() + 
#  scale_x_discrete(labels = c("Stasis", "URW", "GRW", "Accel.", "Decel.", "OU", "OU mov. opt.")) +
#  labs(fill = "") + scale_fill_discrete(name = "Environment", labels = c("Fluvial", "Lacustrine", "Marine",
#                                                                         "Terrestrial"), palette = col_val2) +
#  xlab("Model") + ylab("Count") + theme(legend.text = element_text(size = 10), 
#                                        legend.title = element_text(size = 12),
#                                        axis.title = element_text(size = 12))
#dev.off()

=======
# micro vs. macro
#micro_macro <- plot_data
#level_order <- c("stasis", "URW", "GRW", "accel", "decel", "OU", "OU mov opt")
#pdf("./results_paleoTS_v0.6.1/plot/micro_macro_aicc.pdf")
#ggplot(micro_macro, aes(x = factor(model_aicc, levels = level_order), fill = microfossil)) + geom_bar() +
#  theme_classic() + scale_x_discrete(labels = c("Stasis", "URW", "GRW", "Accel.", "Decel.", "OU", "OU mov. opt.")) +
#  labs(fill = "") + scale_fill_discrete(name = "", labels = c("Macrofossils", "Microfossils"), palette = col_val1) +
#  xlab("Model") + ylab("Count") + theme(legend.text = element_text(size = 15))
#dev.off()

# environment
#env <- plot_data
#level_order <- c("stasis", "URW", "GRW", "accel", "decel", "OU", "OU mov opt")
#pdf("./results_paleoTS_v0.6.1/plot/environment_aicc.pdf")
#ggplot(env, aes(x = factor(model_aicc, levels = level_order), fill = environment)) + geom_bar() +
#  theme_classic() + 
#  scale_x_discrete(labels = c("Stasis", "URW", "GRW", "Accel.", "Decel.", "OU", "OU mov. opt.")) +
#  labs(fill = "") + scale_fill_discrete(name = "Environment", labels = c("Fluvial", "Lacustrine", "Marine",
#                                                                        "Terrestrial"), palette = col_val2) +
#  xlab("Model") + ylab("Count") + theme(legend.text = element_text(size = 10), 
#                                        legend.title = element_text(size = 12),
#                                        axis.title = element_text(size = 12))
#dev.off()
>>>>>>> 15cd3634ef9dd1df70f22064c075a40aa283aa18

###### interval MY ######
intv_my <- plot_data[c("model_aicc", "interval_MY", "parameters")]


# stats 
intv_my_mean <- aggregate(intv_my[, 2], list(intv_my$model_aicc), mean)
intv_my_mean <- as.data.frame(intv_my_mean)
names(intv_my_mean) <- c("model", "mean interval (MY)")
intv_my_median <- aggregate(intv_my[, 2], list(intv_my$model_aicc), median)
intv_my_median <- as.data.frame(intv_my_median)
names(intv_my_median) <- c("model", "median interval (MY)")

sink(file = "./results_paleoTS_v0.6.1/interval_my_results_aicc.txt")
intv_my_mean
intv_my_median
sink()

# LMM
#lmm_model_time <- lmer(interval_MY ~ 1+ model_aicc + (1| popID), data = plot_data)
lmm_model_time <- lm(interval_MY ~ model_aicc, data = plot_data)
summary(lmm_model_time)

# plot
level_order <- c("stasis", "URW", "GRW", "accel", "decel", "OU", "OU mov opt")
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
steps <- plot_data[c("model_aicc", "steps", "parameters")]

# stats
steps_mean <- aggregate(steps[, 2], list(steps$model_aicc), mean)
steps_mean <- as.data.frame(steps_mean)
names(steps_mean) <- c("model", "mean steps")
steps_median <- aggregate(steps[, 2], list(steps$model_aicc), median)
steps_median <- as.data.frame(steps_median)
names(steps_median) <- c("model", "median steps")

sink(file = "./results_paleoTS_v0.6.1/steps_uni_results_aicc.txt")
steps_mean
steps_median
sink()

# LMM
#lmm_model_steps <- lmer(steps ~ model_aicc + (1| popID), data = plot_data)
lmm_model_steps <- lm(steps ~ model_aicc, data = plot_data)
summary(lmm_model_steps)

# plot
level_order <- c("stasis", "URW", "GRW", "accel", "decel", "OU", "OU mov opt")
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
res <- plot_data[c("model_aicc", "steps", "interval_MY", "parameters")]
res$resolution <- res$steps/res$interval_MY

# stats
res_mean <- aggregate(res[, 5], list(res$model_aicc), mean)
res_mean <- as.data.frame(res_mean)
names(res_mean) <- c("model", "mean res")
res_median <- aggregate(res[, 5], list(res$model_aicc), median)
res_median <- as.data.frame(res_median)
names(res_median) <- c("model", "median res")

sink(file = "./results_paleoTS_v0.6.1/resolution_uni_results_aicc.txt")
res_mean
res_median
sink()

# LMM
plot_data$resolution = plot_data$steps/plot_data$interval_MY
#lmm_model_resolution <- lmer(resolution ~ 1+ model_aicc + (1| popID), data = plot_data)
lmm_model_resolution <- lm(resolution ~ model_aicc, data = plot_data)
summary(lmm_model_resolution)

# plot
level_order <- c("stasis", "URW", "GRW", "accel", "decel", "OU", "OU mov opt")
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
pdf(width = 15.0, height = 5.5, file = "./results_paleoTS_v0.6.1/plot/interval_steps_res_aicc.pdf")
grid.arrange(intv_plot,
             steps_plot,
             res_plot, nrow = 1)
dev.off()



# put the output of the linear mixed effects models in the same table
models = rownames(summary(lmm_model_time)$coefficients)
interval = as.data.frame(summary(lmm_model_time)$coefficients)
steps = as.data.frame(summary(lmm_model_steps)$coefficients)
resolution = as.data.frame(summary(lmm_model_resolution)$coefficients)

interval_df <- data.frame(term = models, i.Est = interval$Estimate,
                          i.SE = interval$"Std", i.p.value = interval$"Pr")

steps_df <- data.frame(s.Est = steps$Estimate,
                       s.SE = steps$"Std", s.p.value = steps$"Pr")

resolution_df <- data.frame(r.Est = resolution$Estimate,
                            r.SE = resolution$"Std", r.p.value = resolution$"Pr")

lmm_result_table <- cbind(interval_df, steps_df, resolution_df)

write.csv(lmm_result_table, file = "./results_paleoTS_v0.6.1/table_lmm_uni.pdf", row.names = FALSE)



#-------------------
# PLOT ABSOLUTE FIT
#-------------------


# remove NA data
plot_data2 <- plot_data 
plot_data2 <- plot_data2 %>% drop_na(model_adequate)

<<<<<<< HEAD
# ordering the model according to the number of parameters
level_order <- c("stasis", "URW", "GRW", "accel", "decel", "OU", "OU mov opt")
plot_data2$model_aicc <- factor(plot_data2$model_aicc, levels = level_order)
plot_data2$model_aicc <- relevel(plot_data2$model_aicc, ref = "stasis")

=======
>>>>>>> 15cd3634ef9dd1df70f22064c075a40aa283aa18
# micro vs. macro
#micro_macro <- plot_data2[c("model_adequate", "microfossil")]
#level_order <- c("stasis", "URW", "GRW", "accel", "decel", "OU", "OU mov opt")
#pdf("./results_paleoTS_v0.6.1/plot/micro_macro_adeq.pdf")
#ggplot(micro_macro, aes(x = factor(model_adequate, levels = level_order), fill = microfossil)) + geom_bar() +
#  theme_classic() + scale_x_discrete(labels = c("Stasis", "URW", "GRW", "Accel.", "Decel.", "OU", "OU mov. opt.")) +
#  labs(fill = "") + scale_fill_discrete(name = "", labels = c("Macrofossils", "Microfossils"), palette = col_val1) +
#  xlab("Model") + ylab("Count") + theme(legend.text = element_text(size = 15))
#dev.off()

# environment
#env <- plot_data2[c("model_adequate", "environment")]
#level_order <- c("stasis", "URW", "GRW", "accel", "decel", "OU", "OU mov opt")
#pdf("./results_paleoTS_v0.6.1/plot/environment_adeq.pdf")
#ggplot(env, aes(x = factor(model_adequate, levels = level_order), fill = environment)) + geom_bar() +
#  theme_classic() + 
#  scale_x_discrete(labels = c("Stasis", "URW", "GRW", "Accel.", "Decel.", "OU", "OU mov. opt.")) +
#  labs(fill = "") + scale_fill_discrete(name = "Environment", labels = c("Lacustrine", "Marine",
#                                                                         "Terrestrial"), palette = col_val2[2:4]) +
#  xlab("Model") + ylab("Count") + theme(legend.text = element_text(size = 10), 
#                                        legend.title = element_text(size = 12),
#                                        axis.title = element_text(size = 12))
# dev.off()

###### interval MY ######
intv_my <- plot_data2[c("model_adequate", "interval_MY", "parameters")]


# stats
intv_my_mean <- aggregate(intv_my[, 2], list(intv_my$model_adequate), mean)
intv_my_mean <- as.data.frame(intv_my_mean)
names(intv_my_mean) <- c("model", "mean interval (MY)")
intv_my_median <- aggregate(intv_my[, 2], list(intv_my$model_adequate), median)
intv_my_median <- as.data.frame(intv_my_median)
names(intv_my_median) <- c("model", "median interval (MY)")

sink(file = "./results_paleoTS_v0.6.1/interval_my_results_adeq.txt")
intv_my_mean
intv_my_median
sink()

# LMM
#lmm_model_time2 <- lmer(interval_MY ~ 1+ model_aicc + (1| popID), data = plot_data2)
lmm_model_time2 <- lm(interval_MY ~ model_aicc, data = plot_data2)
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

# stats
steps_mean <- aggregate(steps[, 2], list(steps$model_adequate), mean)
steps_mean <- as.data.frame(steps_mean)
names(steps_mean) <- c("model", "mean steps")
steps_median <- aggregate(steps[, 2], list(steps$model_adequate), median)
steps_median <- as.data.frame(steps_median)
names(steps_median) <- c("model", "median steps")

sink(file = "./results_paleoTS_v0.6.1/steps_uni_results_adeq.txt")
steps_mean
steps_median
sink()

# LMM
#lmm_model_steps2 <- lmer(steps ~ model_aicc + (1| popID), data = plot_data2)
lmm_model_steps2 <- lm(steps ~ model_aicc, data = plot_data2)
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

# stats
res_mean <- aggregate(res[, 5], list(res$model_adequate), mean)
res_mean <- as.data.frame(res_mean)
names(res_mean) <- c("model", "mean res")
res_median <- aggregate(res[, 5], list(res$model_adequate), median)
res_median <- as.data.frame(res_median)
names(res_median) <- c("model", "median res")

sink(file = "./results_paleoTS_v0.6.1/resolution_uni_results_adeq.txt")
res_mean
res_median
sink()

# LMM
plot_data2$resolution = plot_data2$steps/plot_data2$interval_MY
#lmm_model_resolution2 <- lmer(resolution ~ 1+ model_aicc + (1| popID), data = plot_data2)
lmm_model_resolution2 <- lm(resolution ~ model_aicc, data = plot_data2)
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

write.csv(lmm_result_table2, file = "./results_paleoTS_v0.6.1/table_lmm_uni_adeq.pdf", row.names = FALSE)
