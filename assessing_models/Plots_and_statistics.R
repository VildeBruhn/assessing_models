###########################################
##                                       ##
##  Plots and statistics on the results  ##
##                                       ##
###########################################

#R version 4.2.1
#evoTS version 1.0.3
#adePEM new models version
#paleoTS version 0.6.2

rm(list = ls())

library(lme4)
library(lmerTest)
library(readr)
library(dplyr)
library(tidyverse)
library(gridExtra)
library(ggplot2)

########################
# Models without shift #
########################

#Import results of the test (no models with double processes, no adequacy)
load("./model_test_uni.Rdata")

#Import metadata and extract the ones for time series tested
metadata <- read_delim("./timeseries/metadata.txt", col_names = TRUE, delim = "\t")
metadata_model_test <- metadata[metadata$tsID %in% names(model_test), ] 
metadata_model_test$tsID = as.character(metadata_model_test$tsID)

# missing_ids <- setdiff(names(model_test), metadata_model_test$tsID)
# model_test <- model_test[!names(model_test) %in% missing_ids]

#Add the model which has the minimal AICc for each time series to the metadata
model_aicc_min <- lapply(model_test, function(aicc_min) {
  which.min(as.numeric(unlist(aicc_min$AICc)))
})

metadata_model_test$best_model_numeric <- model_aicc_min[metadata_model_test$tsID]

metadata_model_test$best_model <- 0
metadata_model_test <- metadata_model_test %>%
  mutate(best_model = case_when(
    best_model_numeric == 1 ~ 'GRW',
    best_model_numeric == 2 ~ 'URW',
    best_model_numeric == 3 ~ 'Stasis',
    best_model_numeric == 4 ~ 'Stasis',
    best_model_numeric == 5 ~ 'Decel',
    best_model_numeric == 6 ~ 'Accel',
    best_model_numeric == 7 ~ 'OU',
    best_model_numeric == 8 ~ 'OU_mov_opt',
    best_model_numeric == 9 ~ 'OU_mov_opt',
    TRUE ~ as.character(best_model)
  ))

model_order <- c("Stasis", "URW", "GRW", "Accel", "Decel", "OU", "OU_mov_opt")
metadata_model_test$best_model <- factor(metadata_model_test$best_model, levels = model_order)

histogram_main_results = ggplot(metadata_model_test, aes(x = best_model, fill = best_model)) +
  geom_bar(color = "black") +
  scale_fill_manual(values = rainbow(length(unique(metadata_model_test$best_model)))) +
  labs(title = "Time series best fitted for each model",
       x = "Models",
       y = "Time series") +
  scale_y_continuous(limits = c(0, 175)) +
  theme(
  panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
  plot.title = element_text(size = 38),
  axis.title = element_text(size = 38),
  axis.text = element_text(size = 36),
  legend.text = element_text(size = 36)
)

# save the graphs
png("./results_paleoTS_v0.6.1/plot/histogram_main_results.png", width = 2000, height = 1500)
grid.arrange(histogram_main_results, nrow = 1)
dev.off()


#-----------------------------------------------------
# Correlation between best model and number of steps 
#-----------------------------------------------------

lmm_model_steps <- lmer(steps ~ best_model + (1| popID), data = metadata_model_test)
summary(lmm_model_steps)

# Create the Boxplot
boxplot_model_steps = ggplot(metadata_model_test, aes(x = best_model, y = steps, fill = best_model)) +
  geom_boxplot() +
  labs(title = "Comparison of time series resolution (number of steps) with their best model",
       x = "Best model",
       y = "Number of steps") +
  theme_classic() +
  scale_y_log10() +
  scale_fill_manual(values = rainbow(length(unique(metadata_model_test$best_model)))) +
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
    plot.title = element_text(size = 38),
    axis.title = element_text(size = 38),
    axis.text = element_text(size = 36),
    legend.text = element_text(size = 36)
  )


#--------------------------------------------------------
# Correlation between best model and total time interval 
#--------------------------------------------------------

lmm_model_time <- lmer(interval_MY ~ 1+ best_model + (1| popID), data = metadata_model_test)
summary(lmm_model_time)

# Create the Boxplot
boxplot_model_time = ggplot(metadata_model_test, aes(x = best_model, y = interval_MY, fill = best_model)) +
  geom_boxplot() +
  labs(title = "Comparison of time series length (total time interval) with their best model",
       x = "Best model",
       y = "Total time interval in MY") +
  theme_classic() +
  scale_fill_manual(values = rainbow(length(unique(metadata_model_test$best_model)))) +
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
    plot.title = element_text(size = 38),
    axis.title = element_text(size = 38),
    axis.text = element_text(size = 36),
    legend.text = element_text(size = 36)
  )


#------------------------------------------------
# Correlation between best model and resolution
#------------------------------------------------

metadata_model_test$resolution = metadata_model_test$steps/metadata_model_test$interval_MY

lmm_model_time <- lmer(resolution ~ 1+ best_model + (1| popID), data = metadata_model_test)
summary(lmm_model_time)

# Create the Boxplot
boxplot_model_resolution = ggplot(metadata_model_test, aes(x = best_model, y = resolution, fill = best_model)) +
  geom_boxplot() +
  labs(title = "Comparison of time series length (total time interval) with their best model",
       x = "Best model",
       y = "Resolution (step/interval)") +
  theme_classic() +
  scale_y_log10() +
  scale_fill_manual(values = rainbow(length(unique(metadata_model_test$best_model)))) +
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
    plot.title = element_text(size = 38),
    axis.title = element_text(size = 38),
    axis.text = element_text(size = 36),
    legend.text = element_text(size = 36)
  )

# save the graphs
png("./results_paleoTS_v0.6.1/plot/boxplot_model.png", width = 6000, height = 1500)
grid.arrange(boxplot_model_steps, boxplot_model_time, boxplot_model_resolution, nrow = 1)
dev.off()




#####################
# Models with shift #
#####################

#Import results of the test (no models with double processes, no adequacy)
load("./results_paleoTS_v0.6.1/Results_fit_shiftmodels.RData")

#Import metadata and extract the ones for time series tested
metadata <- read_delim("./timeseries/metadata.txt", col_names = TRUE, delim = "\t")
metadata_shiftmodel_test <- metadata[metadata$tsID %in% names(model_shift_results), ] 
metadata_shiftmodel_test$tsID = as.character(metadata_shiftmodel_test$tsID)

# missing_ids <- setdiff(names(model_shift_results), metadata_shiftmodel_test$tsID)
# model_shift_results <- model_shift_results[!names(model_shift_results) %in% missing_ids]

#Add the model which has the minimal AICc for each time series to the metadata
metadata_shiftmodel_test <- metadata_shiftmodel_test %>% filter(tsID %in% names(aicc_min))
metadata_shiftmodel_test$best_model_numeric <- unlist(aicc_min[metadata_shiftmodel_test$tsID])

metadata_shiftmodel_test$best_model <- 0
metadata_shiftmodel_test <- metadata_shiftmodel_test %>%
  mutate(best_model = case_when(
    best_model_numeric == 1 ~ 'GRW',
    best_model_numeric == 2 ~ 'URW',
    best_model_numeric == 3 ~ 'Stasis',
    best_model_numeric == 4 ~ 'Strict_stasis',
    best_model_numeric == 5 ~ 'Decel',
    best_model_numeric == 6 ~ 'Accel',
    best_model_numeric == 7 ~ 'OU',
    best_model_numeric == 8 ~ 'OU_mov_opt_anc',
    best_model_numeric == 9 ~ 'OU_mov_opt',
    best_model_numeric == 10 ~ 'Stasis_Stasis',
    best_model_numeric == 11 ~ 'Stasis_URW',
    best_model_numeric == 12 ~ 'Stasis_GRW',
    best_model_numeric == 13 ~ 'Stasis_OU',
    best_model_numeric == 14 ~ 'URW_URW',
    best_model_numeric == 15 ~ 'URW_GRW',
    best_model_numeric == 16 ~ 'URW_OU',
    best_model_numeric == 17 ~ 'GRW_GRW',
    best_model_numeric == 18 ~ 'GRW_OU',
    best_model_numeric == 19 ~ 'OU_OU',
    best_model_numeric == 20 ~ 'OU_GRW',
    best_model_numeric == 21 ~ 'OU_URW',
    best_model_numeric == 22 ~ 'OU_Stasis',
    best_model_numeric == 23 ~ 'GRW_URW',
    best_model_numeric == 24 ~ 'GRW_Stasis',
    best_model_numeric == 25 ~ 'URW_Stasis',
    TRUE ~ as.character(best_model)
  ))

#Add a column to indicate if the best model is with or without shift
metadata_shiftmodel_test$shift <- ifelse(metadata_shiftmodel_test$best_model_numeric >= 1 & metadata_shiftmodel_test$best_model_numeric <= 9, "without shift", "with shift")


#--------------------------------------------------------
# Correlation between type of model and number of steps 
#--------------------------------------------------------

lmm_shift_steps <- lmer(steps ~ shift + (1| popID), data = metadata_shiftmodel_test)
summary(lmm_shift_steps)

shift_order <- c("without shift", "with shift")
metadata_shiftmodel_test$shift <- factor(metadata_shiftmodel_test$shift, levels = shift_order)

# Create the Boxplot
boxplot_shiftmodel_steps = ggplot(metadata_shiftmodel_test, aes(x = shift, y = steps, fill = shift)) +
  geom_boxplot() +
  labs(title = "Comparison of time series resolution (number of steps) with their best model",
       x = "Best model",
       y = "Number of steps") +
  theme_classic() +
  scale_y_log10() +
  scale_fill_manual(values = c("without shift" = "darkgreen", "with shift" = "purple")) +
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
    plot.title = element_text(size = 38),
    axis.title = element_text(size = 38),
    axis.text = element_text(size = 36),
    legend.text = element_text(size = 36)
  )


#-----------------------------------------------------------
# Correlation between type of model and total time interval 
#-----------------------------------------------------------

lmm_shift_time <- lmer(interval_MY ~ 1+ shift + (1| popID), data = metadata_shiftmodel_test)
summary(lmm_shift_time)

# Create the Boxplot
boxplot_shiftmodel_time = ggplot(metadata_shiftmodel_test, aes(x = shift, y = interval_MY, fill = shift)) +
  geom_boxplot() +
  labs(title = "Comparison of time series length (total time interval) with their best model",
       x = "Best model",
       y = "Total time interval in MY") +
  theme_classic() +
  scale_fill_manual(values = c("without shift" = "darkgreen", "with shift" = "purple")) +
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
    plot.title = element_text(size = 38),
    axis.title = element_text(size = 38),
    axis.text = element_text(size = 36),
    legend.text = element_text(size = 36)
  )


#---------------------------------------------------
# Correlation between type of model and resolution
#---------------------------------------------------

metadata_shiftmodel_test$resolution = metadata_shiftmodel_test$steps/metadata_shiftmodel_test$interval_MY

lmm_shift_time <- lmer(resolution ~ 1+ shift + (1| popID), data = metadata_shiftmodel_test)
summary(lmm_shift_time)

# Create the Boxplot
boxplot_shiftmodel_resolution = ggplot(metadata_shiftmodel_test, aes(x = shift, y = resolution, fill = shift)) +
  geom_boxplot() +
  labs(title = "Comparison of time series length (total time interval) with their best model",
       x = "Best model",
       y = "Resolution (step/interval)") +
  theme_classic() +
  scale_y_log10() +
  scale_fill_manual(values = c("without shift" = "darkgreen", "with shift" = "purple")) +
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
    plot.title = element_text(size = 38),
    axis.title = element_text(size = 38),
    axis.text = element_text(size = 36),
    legend.text = element_text(size = 36)
  )

# save the graphs
png("./results_paleoTS_v0.6.1/plot/boxplot_shiftmodel.png", width = 6000, height = 1500)
grid.arrange(boxplot_shiftmodel_steps, boxplot_shiftmodel_time, boxplot_shiftmodel_resolution, nrow = 1)
dev.off()
