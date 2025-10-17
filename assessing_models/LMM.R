##################
## Code for LMM ##
##################


#R version 4.2.1
#evoTS version 1.0.3
#adePEM new models version
#paleoTS version 0.6.2

rm(list = ls())

library(lme4)
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
metadata_adeq_passed_trait$best_model <- metadata_adeq_passed$best_model[metadata_adeq_passed$tsID %in% metadata_adeq_passed_trait$tsID]


#--------------------------------------------------------------------------------------------
# Correlation between best fitted model and total time interval within adequate time series
#--------------------------------------------------------------------------------------------

lmm_model_time <- lmer(interval_MY ~ best_model + (1| popID), data = metadata_adeq_passed)
summary(lmm_model_time)
plot(lmm_model_time)
print(metadata_adeq_passed)

#--------------------------------------------------------------------------------------------
# Correlation between best fitted model and number of steps within adequate time series
#--------------------------------------------------------------------------------------------

lmm_model_steps <- lmer(steps ~ best_model + (1| popID), data = metadata_adeq_passed)
summary(lmm_model_steps)
plot(lmm_model_steps)


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

# adding a column for the aicc gap between the two best fitting models
metadatalong2$tsID = as.character(metadatalong2$tsID)
metadatalong2$aicc_bestmodels_gap <- sapply(metadatalong2$tsID, function(tsID) if (tsID %in% names(aicc_mingap_list)) aicc_mingap_list[[(tsID)]] else NA)

# get a metadata table with only adequate time series
metadatalong2_adeq_passed <- metadatalong2[metadatalong2$tsID %in% list_adeq_passed_shift, ]


#---------------------------------------------------------------------------------------------------------------------
# Correlation between adequacy status and aicc gap between the two best fitted models for all the long time series
#---------------------------------------------------------------------------------------------------------------------

lmm_adequacy_aiccgap <- lmer(aicc_bestmodels_gap ~ adequacy_status + (1| popID), data = metadatalong2)
summary(lmm_adequacy_aiccgap)
plot(lmm_adequacy_aiccgap)

glmm_adequacy_aiccgap <- glmer(aicc_bestmodels_gap ~ adequacy_status + (1| popID), family = binomial, data = metadatalong2)
summary(glmm_adequacy_aiccgap)
plot(glmm_adequacy_aiccgap)

#-----------------------------------------------------------------------------------------------------------------
# Correlation between best fitted model and aicc gap between the two best fitted models for adequate time series
#-----------------------------------------------------------------------------------------------------------------

lmm_model_aiccgap <- lmer(aicc_bestmodels_gap ~ best_model + (1| popID), data = metadatalong2_adeq_passed)
print(lmm_model_aiccgap, correlation = TRUE)
plot(lmm_model_aiccgap)


#----------------------------------------------------------------------------------------------------------------
# Correlation between best fitted model (with or without shift) and number of steps within adequate time series
#----------------------------------------------------------------------------------------------------------------

lmm_shift_steps <- lmer(steps ~ shift_or_noshift + (1| popID), data = metadatalong2_adeq_passed)
summary(lmm_shift_steps)
plot(lmm_shift_steps)


#----------------------------------------------------------------------------------------------------------------
# Correlation between best fitted model (with or without shift) and time interval within adequate time series
#----------------------------------------------------------------------------------------------------------------

lmm_shift_time <- lmer(interval_MY ~ shift_or_noshift + (1| popID), data = metadatalong2_adeq_passed)
summary(lmm_shift_time)
plot(lmm_shift_time)



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