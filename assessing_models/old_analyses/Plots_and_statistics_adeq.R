###################################################################
##                                                               ##
##  Plots and statistics on the results for adequate time series ##
##                                                               ##
###################################################################

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
load("./adeq_uni_passed.Rdata")
model_test_adeq = c(GRW_adeq_passed, URW_adeq_passed, stasis_adeq_passed, strict_stasis_adeq_passed, 
                    decel_adeq_passed, accel_adeq_passed, OU_adeq_passed, OU_mov_opt_anc_adeq_passed, 
                    OU_mov_opt_adeq_passed)

#Import metadata and extract the ones for time series tested
metadata <- read_delim("./timeseries/metadata.txt", col_names = TRUE, delim = "\t")
metadata_model_test_adeq <- metadata[metadata$tsID %in% names(model_test_adeq), ] 
metadata_model_test_adeq$tsID = as.character(metadata_model_test_adeq$tsID)

# missing_ids <- setdiff(names(model_test), metadata_model_test$tsID)
# model_test <- model_test[!names(model_test) %in% missing_ids]

#Add the model which has the minimal AICc for each time series to the metadata
model_aicc_min <- lapply(model_test, function(aicc_min) {
  which.min(as.numeric(unlist(aicc_min$AICc)))
})

model_aicc_min <- model_aicc_min[names(model_aicc_min) %in% names(model_test_adeq)] 

metadata_model_test_adeq$best_model_numeric <- model_aicc_min[metadata_model_test_adeq$tsID]

metadata_model_test_adeq$best_model <- 0
metadata_model_test_adeq <- metadata_model_test_adeq %>%
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


#-----------------------------------------------------
# Correlation between best model and number of steps 
#-----------------------------------------------------

lmm_model_steps <- lmer(steps ~ best_model + (1| popID), data = metadata_model_test_adeq)
summary(lmm_model_steps)

model_order <- c("Stasis", "URW", "GRW", "Accel", "Decel", "OU", "OU_mov_opt")
metadata_model_test_adeq$best_model <- factor(metadata_model_test_adeq$best_model, levels = model_order)

# Create the Boxplot
boxplot_model_steps_adeq = ggplot(metadata_model_test_adeq, aes(x = best_model, y = steps, fill = best_model)) +
  geom_boxplot() +
  labs(title = "Comparison of time series resolution (number of steps) with their best model",
       x = "Best model",
       y = "Number of steps") +
  theme_classic() +
  scale_y_log10() +
  scale_fill_manual(values = rainbow(length(unique(metadata_model_test_adeq$best_model)))) +
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

lmm_model_time <- lmer(interval_MY ~ 1+ best_model + (1| popID), data = metadata_model_test_adeq)
summary(lmm_model_time)

# Create the Boxplot
boxplot_model_time_adeq = ggplot(metadata_model_test_adeq, aes(x = best_model, y = interval_MY, fill = best_model)) +
  geom_boxplot() +
  labs(title = "Comparison of time series length (total time interval) with their best model",
       x = "Best model",
       y = "Total time interval in MY") +
  theme_classic() +
  scale_fill_manual(values = rainbow(length(unique(metadata_model_test_adeq$best_model)))) +
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

metadata_model_test_adeq$resolution = metadata_model_test_adeq$steps/metadata_model_test_adeq$interval_MY

lmm_model_time <- lmer(resolution ~ 1+ best_model + (1| popID), data = metadata_model_test_adeq)
summary(lmm_model_time)

# Create the Boxplot
boxplot_model_resolution_adeq = ggplot(metadata_model_test_adeq, aes(x = best_model, y = resolution, fill = best_model)) +
  geom_boxplot() +
  labs(title = "Comparison of time series length (total time interval) with their best model",
       x = "Best model",
       y = "Resolution (step/interval)") +
  theme_classic() +
  scale_y_log10() +
  scale_fill_manual(values = rainbow(length(unique(metadata_model_test_adeq$best_model)))) +
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
    plot.title = element_text(size = 38),
    axis.title = element_text(size = 38),
    axis.text = element_text(size = 36),
    legend.text = element_text(size = 36)
  )

# save the graphs
png("./results_paleoTS_v0.6.1/plot/boxplot_model_adeq.png", width = 6000, height = 1500)
grid.arrange(boxplot_model_steps_adeq, boxplot_model_time_adeq, boxplot_model_resolution_adeq, nrow = 1)
dev.off()




#####################
# Models with shift #
#####################

#Import results of the test (no models with double processes, no adequacy)
load("./model_test_shift.Rdata")
load("./ln_data_meta_shift.Rdata")
load("./adeq_shift_passed.Rdata")

#load("results_paleoTS_v0.6.1/Results_fit_shiftmodels.RData")
#load("results_paleoTS_v0.6.1/Results_fit_adequacy_shiftmodels.RData")

model_test_shift_adeq = c(GRW_adeq_passed, URW_adeq_passed, stasis_adeq_passed, strict_stasis_adeq_passed, decel_adeq_passed,
accel_adeq_passed, OU_adeq_passed, OU_mov_opt_anc_adeq_passed, OU_mov_opt_adeq_passed, Stasis_Stasis_adeq_passed,
Stasis_URW_adeq_passed, Stasis_GRW_adeq_passed, Stasis_OU_adeq_passed, URW_URW_adeq_passed, URW_GRW_adeq_passed,
URW_OU_adeq_passed, GRW_GRW_adeq_passed, GRW_OU_adeq_passed, OU_OU_adeq_passed, OU_GRW_adeq_passed, OU_URW_adeq_passed,
OU_Stasis_adeq_passed, GRW_URW_adeq_passed, GRW_Stasis_adeq_passed, URW_Stasis_adeq_passed)

#Import metadata and extract the ones for time series tested
ln_data_meta_shift_adeq <- ln_data_meta_shift[ln_data_meta_shift_adeq$tsID %in% names(model_test_shift_adeq), ] 
ln_data_meta_shift_adeq$tsID = as.character(ln_data_meta_shift_adeq$tsID)

# missing_ids <- setdiff(names(model_shift_results), ln_data_meta_shift_adeq$tsID)
# model_shift_results <- model_shift_results[!names(model_shift_results) %in% missing_ids]

#Add the model which has the minimal AICc for each time series to the metadata
ln_data_meta_shift_adeq <- ln_data_meta_shift_adeq %>% filter(tsID %in% names(aicc_min))
ln_data_meta_shift_adeq$best_model_numeric <- unlist(aicc_min[ln_data_meta_shift_adeq$tsID])

ln_data_meta_shift_adeq$best_model <- 0
ln_data_meta_shift_adeq <- ln_data_meta_shift_adeq %>%
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
ln_data_meta_shift_adeq$shift <- ifelse(ln_data_meta_shift_adeq$best_model_numeric >= 1 & ln_data_meta_shift_adeq$best_model_numeric <= 9, "without shift", "with shift")

#--------------------------------------------------------
# Correlation between type of model and number of steps 
#--------------------------------------------------------

lmm_shift_steps <- lmer(steps ~ shift + (1| popID), data = ln_data_meta_shift_adeq)
summary(lmm_shift_steps)

shift_order <- c("without shift", "with shift")
ln_data_meta_shift_adeq$shift <- factor(ln_data_meta_shift_adeq$shift, levels = shift_order)

# Create the Boxplot
boxplot_shiftmodel_steps_adeq = ggplot(ln_data_meta_shift_adeq, aes(x = shift, y = steps, fill = shift)) +
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

lmm_shift_time <- lmer(interval_MY ~ 1+ shift + (1| popID), data = ln_data_meta_shift_adeq)
summary(lmm_shift_time)

# Create the Boxplot
boxplot_shiftmodel_time_adeq = ggplot(ln_data_meta_shift_adeq, aes(x = shift, y = interval_MY, fill = shift)) +
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

ln_data_meta_shift_adeq$resolution = ln_data_meta_shift_adeq$steps/ln_data_meta_shift_adeq$interval_MY

lmm_shift_time <- lmer(resolution ~ 1+ shift + (1| popID), data = ln_data_meta_shift_adeq)
summary(lmm_shift_time)

# Create the Boxplot
boxplot_shiftmodel_resolution_adeq = ggplot(ln_data_meta_shift_adeq, aes(x = shift, y = resolution, fill = shift)) +
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
png("./results_paleoTS_v0.6.1/plot/boxplot_shiftmodel_adeq.png", width = 6000, height = 1500)
grid.arrange(boxplot_shiftmodel_steps_adeq, boxplot_shiftmodel_time_adeq, boxplot_shiftmodel_resolution_adeq, nrow = 1)
dev.off()


#------------------------------------------------------------------
# How many time series have a best model better than 2 AICc units
#------------------------------------------------------------------

#Create a metadata with time series which passed and did not pass the adequacy tests
metadata_shiftmodel_test <- metadata[metadata$tsID %in% names(model_shift_results), ] 
metadata_shiftmodel_test$tsID = as.character(metadata_shiftmodel_test$tsID)

#Add adequacy status to the metadata
metadata_shiftmodel_test$adequacy_status <- ifelse(metadata_shiftmodel_test$tsID %in% ln_data_meta_shift_adeq$tsID, "adequate", "inadequate")

#Extract DeltaAICc for the second best model
aicc_gap <- list()
aicc_mingap <- list()

for (i in 1:length(aicc)) {
  for (j in 1:nrow(aicc[[1]][])) {
    if (aicc[[i]][j,] != aicc[[i]][aicc_min[[i]],]) {
      aicc_gap[[j]] = abs(aicc[[i]][j,] - aicc[[i]][aicc_min[[i]],])
    } 
  }
  aicc_mingap[i] <- min(unlist(aicc_gap))
  names(aicc_mingap)[i] <- names(aicc_min)[i]
}

missing_ids <- setdiff(names(aicc_mingap), metadata_shiftmodel_test$tsID)
aicc_mingap <- aicc_mingap[!names(aicc_mingap) %in% missing_ids]

metadata_shiftmodel_test$deltaAICc <- unlist(aicc_mingap[metadata_shiftmodel_test$tsID])

#How many time series have a deltaAICc inferior or equal to 2
TS_deltaAICc2_c <- sum(metadata_shiftmodel_test$deltaAICc <= 2)
TS_deltaAICc2_p <- TS_deltaAICc2_c/nrow(metadata_shiftmodel_test)*100

adequacy_order <- c("adequate", "inadequate")
metadata_shiftmodel_test$adequacy_status <- factor(metadata_shiftmodel_test$adequacy_status, levels = adequacy_order)

# Create the Boxplot
boxplot_adequacy_deltaaicc = ggplot(metadata_shiftmodel_test, aes(x = adequacy_status, y = deltaAICc, fill = adequacy_status)) +
  geom_boxplot() +
  labs(title = "Comparison of the deltaAICc between time series which passed and failed adequacy tests",
       x = "deltaAICc gap (second best model - first best model)",
       y = "Adequacy status") +
  theme_classic() +
  scale_fill_manual(values = c("adequate" = "darkgreen", "inadequate" = "purple")) +
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
    plot.title = element_text(size = 38),
    axis.title = element_text(size = 38),
    axis.text = element_text(size = 36),
    legend.text = element_text(size = 36)
  )

# save the graphs
png("./results_paleoTS_v0.6.1/plot/boxplot_adequacy_deltaAICc.png", width = 2000, height = 1500)
grid.arrange(boxplot_adequacy_deltaaicc, nrow = 1)
dev.off()



#------------------------------------
# TABLE FOR RESULTS WITH ALL MODELS
#------------------------------------

# get the number of parameters for each model
k_noshift <- model_noshift_results[[1]]$K
names(k_noshift) <- row.names(model_noshift_results[[1]])

k_shift <- sapply(model_shift_results[[1]], function(model) model$K)
names(k_shift) <- sapply(model_shift_results[[1]], function(model) model$modelName)

k_list = c(k_noshift, k_shift)
names(k_list) = c(names(k_noshift), names(k_shift))
new_names <- c("GRW" = "GRW", "URW" = "URW", "Stasis" = "Stasis", "StrictStasis" = "Strict_stasis", "Decel" = "Decel", "Accel" = "Accel", "OU" = "OU",
               "OU model with moving optimum (ancestral state at optimum)" = "OU_mov_opt_anc", "OU model with moving optimum" = "OU_mov_opt",
               "Stasis-Stasis" = "Stasis_Stasis", "Stasis-URW" = "Stasis_URW", "Stasis-GRW" = "Stasis_GRW", "Stasis-OU" = "Stasis_OU", "URW-URW" = "URW_URW",
               "URW-GRW" = "URW_GRW", "URW-OU" = "URW_OU", "GRW-GRW" = "GRW_GRW", "GRW-OU" = "GRW_OU", "OU-OU" = "OU_OU", "URW-Stasis" = "URW_Stasis",
               "GRW-Stasis" = "GRW_Stasis", "OU-Stasis" = "OU_Stasis", "GRW-URW" = "GRW_URW", "OU-URW" = "OU_URW", "OU-GRW" = "OU_GRW")
  
names(k_list) <- new_names[names(k_list)]
models_list <- names(k_list)

#Relative count - number of time series that are fitted the best by each model before adequacy
model_relative_count <- c()
for (model in model_list) {
  model_TS <- mget(model)
  model_relative_count[model] <- length(model_TS[[model]])
}

#############################################################
model_relative_count <- c()
for (model in model_list) {
  model_relative_count[model] <- length(model)
}

length(model_TS[[GRW_URW]])
length(URW_OU[["URW_OU"]])
length(GRW_URW)
############################################################

#Relative percentage - time series that are fitted the best by each model over the total number of time series before adequacy
model_relative_percentage <- c()
for (model in model_list) {
  model_relative_percentage[model] <- model_relative_count[[model]]/nrow(metadata_shiftmodel_test)*100
}

#Adequacy passed percentage for each model
model_adeq_passed_percentage <- c()
for (model in model_list) {
  model_TS <- mget(model)
  model_TS_adeq_pass_name <- paste0(model, "_adeq_passed")
  model_TS_adeq_pass <- mget(model_TS_adeq_pass_name)
  model_adeq_passed_percentage[model] <- length(model_TS_adeq_pass[[1]])/length(model_TS[[model]])*100
}

#Absolute count - number of time series that are fitted the best by each model after adequacy
model_absolute_count_adeq <- c()
for (model in model_list) {
  model_TS_adeq_pass_name <- paste0(model, "_adeq_passed")
  model_TS_adeq_pass <- mget(model_TS_adeq_pass_name)
  model_absolute_count_adeq[model] <- length(model_TS_adeq_pass[[1]])
}

#Absolute percentage - time series that are fitted the best by each model over the total number of time series which passed adequacy
model_absolute_percentage_adeq <- c()
for (model in model_list) {
  model_TS_adeq_pass_name <- paste0(model, "_adeq_passed")
  model_TS_adeq_pass <- mget(model_TS_adeq_pass_name)
  model_absolute_percentage_adeq[model] <- length(model_TS_adeq_pass[[1]])/nrow(ln_data_meta_shift_adeq)*100
}

#Make output table
resultshift_table <- data.frame(
  model = model_list, nbr_parameters = k_list, relative_count = model_relative_count, relative_percentage = model_relative_percentage,
  adequacy_passed_percentage = model_adeq_passed_percentage, absolute_count_adeq = model_absolute_count_adeq, absolute_percentage_adeq = model_absolute_percentage_adeq)
  
#Write to file
sink(file = "./results_paleoTS_v0.6.1/table_result_shift.txt")
resultshift_table
paste("Total number of time series investigated:", length(model_noshift_results))
paste("Total number of time series which passed adequacy tests:", Total_adeq_passed)
paste("Percentage of time series which passed adequacy tests:", (Total_adeq_passed*100)/length(model_noshift_results))
sink()