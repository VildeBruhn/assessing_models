##################
## Code for LMM ##
##################


#R version 4.2.1
#evoTS version 1.0.3
#adePEM new models version
#paleoTS version 0.6.2

rm(list = ls())

library(lme4)
library(lmerTest)
library(readr)
library(dplyr)

########################
# Models without shift #
########################

# importing metadata
metadata <- read_delim("./timeseries/metadata.txt", col_names = TRUE, delim = "\t")

# importing metadata with trait category
load("./timeseries/metadata_trait_2.RData")

# load adequate timeseries 
load("./adeq_uni_passed.Rdata")

# getting the list of the time series which passed the adequacy tests
list_adeq_passed = c(GRW_adeq_passed, URW_adeq_passed, stasis_adeq_passed, strict_stasis_adeq_passed, 
                     decel_adeq_passed, accel_adeq_passed, OU_adeq_passed, OU_mov_opt_anc_adeq_passed, 
                     OU_mov_opt_adeq_passed)

# keep in the metadata only the time series which passed the adequacy
metadata_adeq_passed_trait <- metadata_trait_organized[metadata_trait_organized$tsID %in% names(list_adeq_passed), ] # Only 309 time series passed the tests and can be categorized in trait and size
metadata_adeq_passed <- metadata[metadata$tsID %in% names(list_adeq_passed), ]

# change the TS id into characters 
metadata_adeq_passed$tsID = as.character(metadata_adeq_passed$tsID)

# add the best model to the metadata for each time series
metadata_adeq_passed$best_model <- NA

for (i in 1:nrow(metadata_adeq_passed)) {
  for (j in 1:length(GRW_adeq_passed)) {
    if (metadata_adeq_passed$tsID[i] == names(GRW_adeq_passed[j])){
      metadata_adeq_passed$best_model[i] = "GRW"
    }
  }
} 

for (i in 1:nrow(metadata_adeq_passed)) {
  for (j in 1:length(URW_adeq_passed)) {
    if (metadata_adeq_passed$tsID[i] == names(URW_adeq_passed[j])){
      metadata_adeq_passed$best_model[i] = "URW"
    }
  }
} 

for (i in 1:nrow(metadata_adeq_passed)) {
  for (j in 1:length(stasis_adeq_passed)) {
    if (metadata_adeq_passed$tsID[i] == names(stasis_adeq_passed[j])){
      metadata_adeq_passed$best_model[i] = "stasis"
    }
  }
} 

for (i in 1:nrow(metadata_adeq_passed)) {
 for (j in 1:length(strict_stasis_adeq_passed)) {
  if (metadata_adeq_passed$tsID[i] == names(strict_stasis_adeq_passed[j])){
  metadata_adeq_passed$best_model[i] = "strict_stasis"
  }
 }
} 

for (i in 1:nrow(metadata_adeq_passed)) {
  for (j in 1:length(decel_adeq_passed)) {
    if (metadata_adeq_passed$tsID[i] == names(decel_adeq_passed[j])){
      metadata_adeq_passed$best_model[i] = "decel"
    }
  }
}

for (i in 1:nrow(metadata_adeq_passed)) {
  for (j in 1:length(accel_adeq_passed)) {
    if (metadata_adeq_passed$tsID[i] == names(accel_adeq_passed[j])){
      metadata_adeq_passed$best_model[i] = "accel"
    }
  }
}

for (i in 1:nrow(metadata_adeq_passed)) {
  for (j in 1:length(OU_adeq_passed)) {
    if (metadata_adeq_passed$tsID[i] == names(OU_adeq_passed[j])){
      metadata_adeq_passed$best_model[i] = "OU"
    }
  }
}

for (i in 1:nrow(metadata_adeq_passed)) {
  for (j in 1:length(OU_mov_opt_anc_adeq_passed)) {
    if (metadata_adeq_passed$tsID[i] == names(OU_mov_opt_anc_adeq_passed[j])){
      metadata_adeq_passed$best_model[i] = "OU_mov_opt_anc"
    }
  }
}

for (i in 1:nrow(metadata_adeq_passed)) {
  for (j in 1:length(OU_mov_opt_adeq_passed)) {
    if (metadata_adeq_passed$tsID[i] == names(OU_mov_opt_adeq_passed[j])){
      metadata_adeq_passed$best_model[i] = "OU_mov_opt"
    }
  }
}

# transfer the best_model column to the matrix with trait category
#metadata_adeq_passed_trait$best_model <- metadata_adeq_passed$best_model[metadata_adeq_passed$tsID %in% metadata_adeq_passed_trait$tsID]


#--------------------------------------------------------------------------------------------
# Correlation between best fitted model and total time interval within adequate time series
#--------------------------------------------------------------------------------------------

lmm_model_time <- lmer(interval_MY ~ 1+ best_model + (1| popID), data = metadata_adeq_passed)
summary(lmm_model_time)

# Create the Boxplot
box_bestfit_length = ggplot(metadata_adeq_passed, aes(x = interaction(best_model), y = interval_MY, fill = best_model)) +
  geom_boxplot() +
  labs(title = "Comparison of time series length (total time interval) with their best model (single process - adequacy passed) ",
       x = "Best model",
       y = "Total time interval in MY") +
  theme_classic() +
  scale_fill_manual(values = rainbow(length(unique(metadata_adeq_passed$best_model)))) +
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
    plot.title = element_text(size = 38),
    axis.title = element_text(size = 38),
    axis.text = element_text(size = 36),
    legend.text = element_text(size = 36)
  )

# save the graphs
png("./results_paleoTS_v0.6.1/plot/results_bestfit_length_box.png", width = 2400, height = 1800)
grid.arrange(box_bestfit_length, ncol = 1)
dev.off()

#--------------------------------------------------------------------------------------------
# Correlation between best fitted model and number of steps within adequate time series
#--------------------------------------------------------------------------------------------

lmm_model_steps <- lmer(steps ~ best_model + (1| popID), data = metadata_adeq_passed)
summary(lmm_model_steps)

# Create the Boxplot
box_bestfit_steps = ggplot(metadata_adeq_passed, aes(x = interaction(best_model), y = steps, fill = best_model)) +
  geom_boxplot() +
  labs(title = "Comparison of time series resolution (number of steps) with their best model (single process - adequacy passed) ",
       x = "Best model",
       y = "Number of steps") +
  theme_classic() +
  scale_fill_manual(values = rainbow(length(unique(metadata_adeq_passed$best_model)))) +
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
    plot.title = element_text(size = 38),
    axis.title = element_text(size = 38),
    axis.text = element_text(size = 36),
    legend.text = element_text(size = 36)
  )

# save the graphs
png("./results_paleoTS_v0.6.1/plot/results_bestfit_steps_box.png", width = 2400, height = 1800)
grid.arrange(box_bestfit_steps, ncol = 1)
dev.off()


# Create a new column with the number of parameters of the best model
metadata_adeq_passed <- metadata_adeq_passed %>%
  mutate(K = case_when(
    best_model == "strict_stasis" ~ 1,
    best_model == "stasis" ~ 2,
    best_model == "URW" ~ 2,
    best_model == "GRW" ~ 3,
    best_model == "accel" ~ 3,
    best_model == "decel" ~ 3,
    best_model == "OU" ~ 4,
    best_model == "OU_mov_opt_anc" ~ 4,
    best_model == "OU_mov_opt" ~ 4))


# Create the Boxplot
plot_bestfitK_steps = ggplot(metadata_adeq_passed, aes(x = K, y = steps, fill = best_model)) +
  geom_boxplot() +
  labs(title = "Comparison of time series resolution (number of steps) with the number of parameter in the best model (single process - adequacy passed) ",
       x = "Number of parameters in the best model (K)",
       y = "Number of steps") +
  theme_classic() +
  scale_fill_manual(values = rainbow(length(unique(metadata_adeq_passed$best_model)))) +
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
    plot.title = element_text(size = 38),
    axis.title = element_text(size = 38),
    axis.text = element_text(size = 36),
    legend.text = element_text(size = 36)
  )

# save the graphs
png("./results_paleoTS_v0.6.1/plot/results_bestfitK_steps_plot.png", width = 2400, height = 1800)
grid.arrange(plot_bestfitK_steps, ncol = 1)
dev.off()



#####################
# Models with shift #
#####################

# load the results
load('./results_paleoTS_v0.6.1/Results_fit_adequacy_shiftmodels.RData')

# names of the time series which have 20 or more steps
TS_long_names <- lapply(result_list, names) 

TS_long_best_model <- list()
for (model_name in names(TS_long_names)) {
  TS <- TS_long_names[[model_name]]
  for (TS in TS) {
    TS_long_best_model[[TS]] <- model_name
  }
}

# adding a column for the best model
metadatalong2 <- metadata[metadata$tsID %in% names(TS_long_best_model), ]

metadatalong2$best_model <- aicc_min

# adding a column for the fit shift/no shift models
metadatalong2$shift_or_noshift <- ifelse(metadatalong2$best_model >= 1 & metadatalong2$best_model <= 9, "without_shift", "with_shift")

# changing the best model into character
metadatalong2 <- metadatalong2 %>%
  mutate(best_model = case_when(
    aicc_min == 1 ~ 'GRW',
    aicc_min == 2 ~ 'URW',
    aicc_min == 3 ~ 'Stasis',
    aicc_min == 4 ~ 'Strict_stasis',
    aicc_min == 5 ~ 'Decel',
    aicc_min == 6 ~ 'Accel',
    aicc_min == 7 ~ 'OU',
    aicc_min == 8 ~ 'OU_mov_opt_anc',
    aicc_min == 9 ~ 'OU_mov_opt',
    aicc_min == 10 ~ 'Stasis_Stasis',
    aicc_min == 11 ~ 'Stasis_URW',
    aicc_min == 12 ~ 'Stasis_GRW',
    aicc_min == 13 ~ 'Stasis_OU',
    aicc_min == 14 ~ 'URW_URW',
    aicc_min == 15 ~ 'URW_GRW',
    aicc_min == 16 ~ 'URW_OU',
    aicc_min == 17 ~ 'GRW_GRW',
    aicc_min == 18 ~ 'GRW_OU',
    aicc_min == 19 ~ 'OU_OU',
    aicc_min == 20 ~ 'OU_GRW',
    aicc_min == 21 ~ 'OU_URW',
    aicc_min == 22 ~ 'OU_Stasis',
    aicc_min == 23 ~ 'GRW_URW',
    aicc_min == 24 ~ 'GRW_Stasis',
    aicc_min == 25 ~ 'URW_Stasis',
    TRUE ~ as.character(best_model)
  ))

# get list of time series which passed adequacy
list_adeq_passed_shift <- c(names(GRW_adeq_passed), names(URW_adeq_passed), names(stasis_adeq_passed), names(strict_stasis_adeq_passed),
                  names(decel_adeq_passed), names(accel_adeq_passed), names(OU_adeq_passed), names(OU_mov_opt_anc_adeq_passed),
                  names(OU_mov_opt_adeq_passed), names(Stasis_Stasis_adeq_passed), names(Stasis_URW_adeq_passed), names(Stasis_GRW_adeq_passed),
                  names(Stasis_OU_adeq_passed), names(URW_URW_adeq_passed), names(URW_GRW_adeq_passed), names(URW_OU_adeq_passed),
                  names(GRW_GRW_adeq_passed), names(GRW_OU_adeq_passed), names(OU_OU_adeq_passed), names(OU_GRW_adeq_passed),
                  names(OU_URW_adeq_passed), names(OU_Stasis_adeq_passed), names(GRW_URW_adeq_passed), names(GRW_Stasis_adeq_passed), names(URW_Stasis_adeq_passed))

# adding a column for the fit shift/no shift models
metadatalong2$adequacy_status <- ifelse(metadatalong2$tsID %in% list_adeq_passed_shift[], "adequate", "non_adequate")

# getting the AICc gap between the two first best fitting models
aicc_gap <- list()
aicc_mingap <- list()
aicc_mingap_list <- list()

for (i in 1:length(aicc)) {
  for (j in 1:nrow(aicc[[1]])) {
    if (aicc[[i]][j,] != aicc[[i]][aicc_min[[i]],]) {
      aicc_gap[[j]] = abs(aicc[[i]][j,] - aicc[[i]][aicc_min[[i]],])
    } 
  }
  aicc_mingap <- min(unlist(aicc_gap, use.names = FALSE))
  sublist_name <- names(aicc)[i]
  aicc_mingap_list[[sublist_name]] <- aicc_mingap
  aicc_gap <- list()
}

# getting the number of best fitting models equal to a gap of 2 aicc
nbr_good_fit <- 0
nbr_good_fit2 <- list()
nbr_good_fit2_list <- list()

for (i in 1:length(aicc)) {
  for (j in 1:nrow(aicc[[1]])) {
    if (aicc[[i]][j, ] - aicc[[i]][aicc_min[[i]], ] <= 2)  {
      nbr_good_fit = nbr_good_fit + 1
    }
  }
  nbr_good_fit2 <- min(unlist(nbr_good_fit, use.names = FALSE))
  sublist_name <- names(aicc)[i]
  nbr_good_fit2_list[[sublist_name]] <- nbr_good_fit2
  nbr_good_fit <- 0
}

# getting the number of best fitting models equal to a gap of 5 aicc
nbr_good_fit5 <- list()
nbr_good_fit5_list <- list()

for (i in 1:length(aicc)) {
  for (j in 1:nrow(aicc[[1]])) {
    if (aicc[[i]][j, ] - aicc[[i]][aicc_min[[i]], ] <= 5)  {
      nbr_good_fit = nbr_good_fit + 1
    }
  }
  nbr_good_fit5 <- min(unlist(nbr_good_fit, use.names = FALSE))
  sublist_name <- names(aicc)[i]
  nbr_good_fit5_list[[sublist_name]] <- nbr_good_fit5
  nbr_good_fit <- 0
}

# adding a column for the aicc gap between the two best fitting models
metadatalong2$tsID = as.character(metadatalong2$tsID)
metadatalong2$aicc_bestmodels_gap <- sapply(metadatalong2$tsID, function(tsID) if (tsID %in% names(aicc_mingap_list)) aicc_mingap_list[[(tsID)]] else NA)

# adding the number of best fitting models for a gap in deltaaicc of 2, and 5.
metadatalong2$nbr_bestmodels_gap2 <- sapply(metadatalong2$tsID, function(tsID) if (tsID %in% names(nbr_good_fit2_list)) nbr_good_fit2_list[[(tsID)]] else NA)
metadatalong2$nbr_bestmodels_gap5 <- sapply(metadatalong2$tsID, function(tsID) if (tsID %in% names(nbr_good_fit5_list)) nbr_good_fit5_list[[(tsID)]] else NA)

# get a metadata table with only adequate time series
metadatalong2_adeq_passed <- metadatalong2[metadatalong2$tsID %in% list_adeq_passed_shift, ]


#---------------------------------------------------------------------------------------------------------------------
# Correlation between adequacy status and aicc gap between the two best fitted models for all the long time series
#---------------------------------------------------------------------------------------------------------------------

lmm_adequacy_aiccgap <- lmer(aicc_bestmodels_gap ~ adequacy_status + (1| popID), data = metadatalong2)
summary(lmm_adequacy_aiccgap)

# Create the Boxplot
box_adequacy_aiccgap = ggplot(metadatalong2, aes(x = adequacy_status, y = aicc_bestmodels_gap, fill = adequacy_status)) +
  geom_boxplot() +
  labs(title = "Comparison of time series aicc gap of best models between adequate and non-adequate time series",
       x = "Adequate or non-adequate",
       y = "aicc gap between the two best models") +
  theme_classic() +
  scale_fill_manual(values = rainbow(length(unique(metadatalong2$adequacy_status)))) +
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
    plot.title = element_text(size = 38),
    axis.title = element_text(size = 38),
    axis.text = element_text(size = 36),
    legend.text = element_text(size = 36)
  )

# save the graphs
png("./results_paleoTS_v0.6.1/plot/result_adequacy_aiccgap_boxplot.png", width = 2600, height = 1800)
grid.arrange(box_adequacy_aiccgap, ncol = 1)
dev.off()


#-----------------------------------------------------------------------------------------------------------------
# Correlation between best fitted model and aicc gap between the two best fitted models for adequate time series
#-----------------------------------------------------------------------------------------------------------------

lmm_model_aiccgap <- lmer(aicc_bestmodels_gap ~ best_model + (1| popID), data = metadatalong2_adeq_passed)
print(lmm_model_aiccgap, correlation = TRUE)
summary(lmm_model_aiccgap)

# Create the Boxplot
plot_bestfit_aiccgap_shift = ggplot(metadatalong2_adeq_passed, aes(x = best_model, y = aicc_bestmodels_gap, fill = best_model)) +
  geom_boxplot() +
  labs(title = "Comparison of time series aicc gap of best models with which model is the best fit",
       x = "best model",
       y = "aicc gap between the two best models") +
  theme_classic() +
  scale_fill_manual(values = rainbow(length(unique(metadatalong2_adeq_passed$best_model)))) +
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
    plot.title = element_text(size = 38),
    axis.title = element_text(size = 38),
    axis.text = element_text(size = 36),
    legend.text = element_text(size = 36)
  )

# save the graphs
png("./results_paleoTS_v0.6.1/plot/results_shift_bestfit_aiccgap_plot.png", width = 10000, height = 1800)
grid.arrange(plot_bestfit_aiccgap_shift, ncol = 1)
dev.off()


#-----------------------------------------------------------------------------------------------------------------
# Correlation between best fitted model and how many fits are good (threshold = 2)
#-----------------------------------------------------------------------------------------------------------------

lmm_model_goodfit2 <- lmer(nbr_bestmodels_gap2 ~ best_model + (1| popID), data = metadatalong2)
summary(lmm_model_goodfit2)

# Create the Boxplot
plot_bestmodel_goodfit2_shift = ggplot(metadatalong2, aes(x = best_model, y = nbr_bestmodels_gap2, fill = best_model)) +
  geom_boxplot() +
  labs(title = "Comparison of the number of good fits (threshold = 2) with which model is the best fit",
       x = "best model",
       y = "how many models are a good fit (threshold = 2)") +
  theme_classic() +
  scale_fill_manual(values = rainbow(length(unique(metadatalong2$best_model)))) +
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
    plot.title = element_text(size = 38),
    axis.title = element_text(size = 38),
    axis.text = element_text(size = 36),
    legend.text = element_text(size = 36)
  )

# save the graphs
png("./results_paleoTS_v0.6.1/plot/results_shift_bestmodel_goodfit2_plot.png", width = 10000, height = 1800)
grid.arrange(plot_bestmodel_goodfit2_shift, ncol = 1)
dev.off()


#-----------------------------------------------------------------------------------------------------------------
# Correlation between best fitted model and how many fits are good (threshold = 5)
#-----------------------------------------------------------------------------------------------------------------

lmm_model_goodfit5 <- lmer(nbr_bestmodels_gap5 ~ best_model + (1| popID), data = metadatalong2)
summary(lmm_model_goodfit5)

# Create the Boxplot
plot_bestmodel_goodfit5_shift = ggplot(metadatalong2, aes(x = best_model, y = nbr_bestmodels_gap5, fill = best_model)) +
  geom_boxplot() +
  labs(title = "Comparison of the number of good fits (threshold = 5) with which model is the best fit",
       x = "best model",
       y = "how many models are a good fit (threshold = 5)") +
  theme_classic() +
  scale_fill_manual(values = rainbow(length(unique(metadatalong2$best_model)))) +
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
    plot.title = element_text(size = 38),
    axis.title = element_text(size = 38),
    axis.text = element_text(size = 36),
    legend.text = element_text(size = 36)
  )

# save the graphs
png("./results_paleoTS_v0.6.1/plot/results_shift_bestmodel_goodfit5_plot.png", width = 10000, height = 1800)
grid.arrange(plot_bestmodel_goodfit5_shift, ncol = 1)
dev.off()

#----------------------------------------------------------------------------------------------------------------
# Correlation between best fitted model (with or without shift) and number of steps within adequate time series
#----------------------------------------------------------------------------------------------------------------

lmm_shift_steps <- lmer(steps ~ shift_or_noshift + (1| popID), data = metadatalong2_adeq_passed)
summary(lmm_shift_steps)

# Create the Boxplot
box_shift_steps = ggplot(metadatalong2_adeq_passed, aes(x = shift_or_noshift, y = steps, fill = shift_or_noshift)) +
  geom_boxplot() +
  labs(title = "Comparison of time series resolution (steps) between time series best fitted by models with or without shift",
       x = "Shift or no shift",
       y = "Number of steps in the time series") +
  theme_classic() +
  scale_y_log10() +
  scale_fill_manual(values = rainbow(length(unique(metadatalong2_adeq_passed$shift_or_noshift)))) +
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
    plot.title = element_text(size = 38),
    axis.title = element_text(size = 38),
    axis.text = element_text(size = 36),
    legend.text = element_text(size = 36)
  )

# save the graphs
png("./results_paleoTS_v0.6.1/plot/result_steps_shift_boxplot.png", width = 2600, height = 1800)
grid.arrange(box_shift_steps, ncol = 1)
dev.off()

# Create a new column with the number of parameters of the best model
metadatalong2_adeq_passed <- metadatalong2_adeq_passed %>%
  mutate(K = case_when(
    best_model == "Strict_stasis" ~ 1,
    best_model == "Stasis" ~ 2,
    best_model == "URW" ~ 2,
    best_model == "GRW" ~ 3,
    best_model == "Accel" ~ 3,
    best_model == "Decel" ~ 3,
    best_model == "OU" ~ 4,
    best_model == "OU_mov_opt_anc" ~ 4,
    best_model == "OU_mov_opt" ~ 4,
    best_model == 'Stasis_Stasis' ~ 4,
    best_model == 'Stasis_URW' ~ 4,
    best_model == 'Stasis_GRW' ~ 5,
    best_model == 'Stasis_OU' ~ 6,
    best_model == 'URW_URW' ~ 4,
    best_model == 'URW_GRW' ~ 5,
    best_model == 'URW_OU' ~ 6,
    best_model == 'GRW_GRW' ~ 6,
    best_model == 'GRW_OU' ~ 7,
    best_model == 'OU_OU' ~ 8,
    best_model == 'OU_GRW' ~ 7,
    best_model == 'OU_URW' ~ 6,
    best_model == 'OU_Stasis' ~ 6,
    best_model == 'GRW_URW' ~ 5,
    best_model == 'GRW_Stasis' ~ 5,
    best_model == 'URW_Stasis' ~ 4))

lmm_model_steps <- lmer(steps ~ best_model + (1| popID), data = metadatalong2_adeq_passed)
summary(lmm_model_steps)

lmm_model_steps <- lmer(steps ~ K + (1| popID), data = metadatalong2_adeq_passed)
summary(lmm_model_steps)

# Create the Boxplot
plot_bestfitK_steps_shift = ggplot(metadatalong2_adeq_passed, aes(x = K, y = steps, fill = best_model)) +
  geom_boxplot() +
  labs(title = "Comparison of time series resolution (number of steps) with the number of parameter in the best model (single process - adequacy passed) ",
       x = "Number of parameters in the best model (K)",
       y = "Number of steps") +
  theme_classic() +
  scale_y_log10() +
  scale_fill_manual(values = rainbow(length(unique(metadatalong2_adeq_passed$best_model)))) +
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
    plot.title = element_text(size = 38),
    axis.title = element_text(size = 38),
    axis.text = element_text(size = 36),
    legend.text = element_text(size = 36)
  )

# save the graphs
png("./results_paleoTS_v0.6.1/plot/results_shift_bestfitK_steps_plot.png", width = 10000, height = 1800)
grid.arrange(plot_bestfitK_steps_shift, ncol = 1)
dev.off()

# Create the plot
linearplot_bestfitK_steps_shift = ggplot(metadatalong2_adeq_passed, aes(x = steps, y = K, color = best_model)) +
  geom_point(size = 10) +
  labs(title = "Comparison of time series resolution (number of steps) with the number of parameter in the best model (double process - adequacy passed) ",
       x = "Number of steps",
       y = "Number of parameters in the best model (K)") +
  theme_classic() +
  scale_x_log10() +
  scale_color_manual(values = rainbow(length(unique(metadatalong2_adeq_passed$best_model)))) + 
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
    plot.title = element_text(size = 38),
    axis.title = element_text(size = 38),
    axis.text = element_text(size = 36),
    legend.text = element_text(size = 36)
  )

# save the graphs
png("./results_paleoTS_v0.6.1/plot/results_shift_bestfitK_steps_linearplot.png", width = 2400, height = 1800)
grid.arrange(linearplot_bestfitK_steps_shift, ncol = 1)
dev.off()


#----------------------------------------------------------------------------------------------------------------
# Correlation between best fitted model (with or without shift) and time interval within adequate time series
#----------------------------------------------------------------------------------------------------------------

lmm_shift_time <- lmer(interval_MY ~ shift_or_noshift + (1| popID), data = metadatalong2_adeq_passed)
summary(lmm_shift_time)

# Create the Boxplot
box_shift_length = ggplot(metadatalong2_adeq_passed, aes(x = shift_or_noshift, y = interval_MY, fill = shift_or_noshift)) +
  geom_boxplot() +
  labs(title = "Comparison of time series length (total time interval) between time series best fitted by models with or without shift",
       x = "Shift or no shift",
       y = "Total time interval in MY") +
  theme_classic() +
  scale_fill_manual(values = rainbow(length(unique(metadatalong2_adeq_passed$shift_or_noshift)))) +
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
    plot.title = element_text(size = 38),
    axis.title = element_text(size = 38),
    axis.text = element_text(size = 36),
    legend.text = element_text(size = 36)
  )

# save the graphs
png("./results_paleoTS_v0.6.1/plot/result_length_shift_boxplot.png", width = 2600, height = 1800)
grid.arrange(box_shift_length, ncol = 1)
dev.off()


lmm_model_steps <- lmer(interval_MY ~ best_model + (1| popID), data = metadatalong2_adeq_passed)
summary(lmm_model_steps)

lmm_model_steps <- lmer(interval_MY ~ K + (1| popID), data = metadatalong2_adeq_passed)
summary(lmm_model_steps)

# Create the Boxplot
box_shift_bestfit_length = ggplot(metadatalong2_adeq_passed, aes(x = K, y = interval_MY, fill = best_model)) +
  geom_boxplot() +
  labs(title = "Comparison of time series length (total time interval) with their best model (single process - adequacy passed) ",
       x = "Best model",
       y = "Total time interval in MY") +
  theme_classic() +
  scale_fill_manual(values = rainbow(length(unique(metadatalong2_adeq_passed$best_model)))) +
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
    plot.title = element_text(size = 38),
    axis.title = element_text(size = 38),
    axis.text = element_text(size = 36),
    legend.text = element_text(size = 36)
  )

# save the graphs
png("./results_paleoTS_v0.6.1/plot/results_shift_bestfit_length_box.png", width = 10000, height = 1800)
grid.arrange(box_shift_bestfit_length, ncol = 1)
dev.off()



####################
# Save the outputs #
####################

# save the test results
sink(file = "./results_paleoTS_v0.6.1/lmm_results.txt")
cat(rep("-", 60), "\n\n")
cat(rep("-", 60), "\n\n")
paste("Tests performed on the", nrow(metadata_adeq_passed), "time series used in the analysis with no shift models (steps >= 10), which passed adequacy test.")
cat(rep("-", 60), "\n\n")
paste("Correlation between best fitted model and total time interval within adequate time series.")
print(summary(lmm_model_time))
cat(rep("-", 60), "\n\n")
paste("Correlation between best fitted model and number of steps within adequate time series.")
print(summary(lmm_model_steps))
cat(rep("-", 60), "\n\n")
cat(rep("-", 60), "\n\n")
paste("Tests performed on the", nrow(metadatalong2), "time series used in the analysis with shift models (steps >= 20), before adequacy test.")
cat(rep("-", 60), "\n\n")
paste("Correlation between adequacy status and aicc gap between the two best fitted models for all the long time series.")
print(summary(lmm_adequacy_aiccgap))
cat(rep("-", 60), "\n\n")
cat(rep("-", 60), "\n\n")
paste("Tests performed on the", nrow(metadatalong2_adeq_passed), "time series used in the analysis with shift models (steps >= 20), which passed adequacy test.")
cat(rep("-", 60), "\n\n")
paste("Correlation between best fitted model and aicc gap between the two best fitted models for adequate time series.")
print(lmm_model_aiccgap, correlation = TRUE)
cat(rep("-", 60), "\n\n")
paste("Correlation between best fitted model (with or without shift) and number of steps within adequate time series.")
print(summary(lmm_shift_steps))
cat(rep("-", 60), "\n\n")
paste("Correlation between best fitted model (with or without shift) and time interval within adequate time series.")
print(summary(lmm_shift_time))
sink()