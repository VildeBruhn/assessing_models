######################
## Code for Boxplot ##
######################

#R version 4.2.1
#evoTS version 1.0.3
#adePEM new models version
#paleoTS version 0.6.2


library(ggplot2)
library(readr)
library(dplyr)
library(gridExtra)


# load the results
load('./results_paleoTS_v0.6.1/Results_fit_adequacy_shiftmodels.RData')
metadata_trait <- read_delim("./timeseries/metadata_envi.txt", col_names = TRUE, delim = "\t")
metadata_envi <- read_delim("./timeseries/metadata_envi.txt", col_names = TRUE, delim = "\t")

#------------------------------------------------------------------------------------------------------------ 
# Correlation between multiple models fitting best one time series (AICc gap < 2) and the adequacy status
#------------------------------------------------------------------------------------------------------------

# Remove problematic time series from metadatalong 
metadatalong_clear <- metadatalong[!(metadatalong$tsID %in% pblm_TS), ]
# remove time series that can not be processed by the loglikelihood function from metadatalong
metadatalong_clear <- metadatalong_clear[!metadatalong_clear$tsID %in% removed_TS, ]


# checking number of adequate time series
Adequacy_passed <- c(GRW_adeq_passed, URW_adeq_passed, stasis_adeq_passed,
                     strict_stasis_adeq_passed, decel_adeq_passed, accel_adeq_passed,
                     OU_adeq_passed, OU_mov_opt_anc_adeq_passed, OU_mov_opt_adeq_passed,
                     Stasis_Stasis_adeq_passed, Stasis_URW_adeq_passed, Stasis_GRW_adeq_passed,
                     Stasis_OU_adeq_passed, URW_URW_adeq_passed, URW_GRW_adeq_passed, URW_OU_adeq_passed,
                     GRW_GRW_adeq_passed, GRW_OU_adeq_passed, OU_OU_adeq_passed, OU_GRW_adeq_passed, OU_URW_adeq_passed,
                     OU_Stasis_adeq_passed, GRW_URW_adeq_passed, GRW_Stasis_adeq_passed, URW_Stasis_adeq_passed)

length(Adequacy_passed)

# adding a new column for adequacy status
metadatalong_clear$tsID <- as.character(metadatalong_clear$tsID)

inadequate_series <- metadatalong_clear[!(
  metadatalong_clear$tsID %in% unlist(
    lapply(Filter(Negate(is.null), list(
      GRW_adeq_passed, URW_adeq_passed, stasis_adeq_passed,
      strict_stasis_adeq_passed, decel_adeq_passed, accel_adeq_passed,
      OU_adeq_passed, OU_mov_opt_anc_adeq_passed, OU_mov_opt_adeq_passed,
      Stasis_Stasis_adeq_passed, Stasis_URW_adeq_passed, Stasis_GRW_adeq_passed,
      Stasis_OU_adeq_passed, URW_URW_adeq_passed, URW_GRW_adeq_passed, URW_OU_adeq_passed,
      GRW_GRW_adeq_passed, GRW_OU_adeq_passed, OU_OU_adeq_passed, OU_GRW_adeq_passed,
      OU_URW_adeq_passed, OU_Stasis_adeq_passed, GRW_URW_adeq_passed,
      GRW_Stasis_adeq_passed, URW_Stasis_adeq_passed
    )), function(lst) names(lst))
  )
), ]

metadatalong_clear$adequacy_status <- ifelse(metadatalong_clear$tsID %in% inadequate_series$tsID, "inadequate", "adequate")


# adding model_fit (TS best described by one or multiple models) to metadata
multiple_models <- metadatalong_clear[metadatalong_clear$tsID %in% names(aicc_filtered), ]
metadatalong_clear$model_fit <- ifelse(metadatalong_clear$tsID %in% multiple_models$tsID, "multiple_models", "unique_model")

# splitting adequate and inadequate time series
adeq_TS = metadatalong_clear[metadatalong_clear$adequacy_status == "adequate", ]
inadeq_TS = metadatalong_clear[metadatalong_clear$adequacy_status == "inadequate", ]

# chi-squared test (447 ts for this test)
contingency_table <- table(metadatalong_clear$adequacy_status, metadatalong_clear$model_fit)
chi_squared_test <- chisq.test(contingency_table)

# list of the minimal aicc gap for each time series (183 ts for this test)
aicc_gap <- list()
aicc_mingap <- list()

for (i in 1:length(aicc)) {
  for (j in 1:nrow(aicc[[1]][])) {
    if (aicc[[i]][j,] != aicc[[i]][aicc_min[[i]],]) {
      aicc_gap[[j]] = abs(aicc[[i]][j,] - aicc[[i]][aicc_min[[i]],])
    } 
  }
  aicc_mingap[i] <- min(unlist(aicc_gap, use.names = FALSE))
  aicc_gap <- list()
}

aicc_mingap_names = names(aicc)
names(aicc_mingap) <- aicc_mingap_names

metadatalong_clear$aicc_mingap = unlist(aicc_mingap)

# anova test
anova_test <- aov(aicc_mingap ~ metadatalong_clear$adequacy_status, data = metadatalong_clear)
summary.lm(anova_test)

# save the test results
sink(file = "./results_paleoTS_v0.6.1/correlation_adequacy_multiplefit_with_shift.txt")
paste("Comparison of adequacy state (passed or failed) in relation to the number of good fit (if AICC gap < 2).")
paste("Tests performed on the", length(metadatalong_clear$tsID), "time series used in the analysis with shift models.")
cat(rep("-", 60), "\n\n")
cat("Chi-squared test between adequacy and models described by more than one model (AICC gap < 2):\n") 
print(chi_squared_test)
cat(rep("-", 60), "\n\n")
cat("Anova test based on the minimum gap between two Aiccs of each time series and their adequacy state:\n")
print(summary.lm(anova_test))
sink()

# boxplot
png("./results_paleoTS_v0.6.1/plot/correlation_adequacy_multiplefit_with_shift_boxplot.png")
ggplot(metadatalong_clear, aes(x = adequacy_status, y = aicc_mingap, fill = adequacy_status)) +
  geom_boxplot() +
  labs(title = "Boxplot of Number of Data Points by Model Type",
       x = "Adequacy state (passed of failed)",
       y = "Gap in the Aicc of the best models") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()



#------------------------------------------------------------------------------------------
# Correlation time series length (number of datapoints) and complex models as best fits
#------------------------------------------------------------------------------------------

data_length <- data.frame(tsID = names(a), best_model = unlist(a))
data_length$nbr_datapoint <- metadatalong$steps[match(data_length$tsID, metadatalong$tsID)]

# Remove problematic time series from metadatalong 
data_length_clear <- data_length[!(data_length$tsID %in% pblm_TS), ]
# remove time series that can not be processed by the loglikelihood function from metadatalong
data_length_clear <- data_length_clear[!data_length_clear$tsID %in% removed_TS, ]

data_length_clear$model_type <- ifelse(data_length_clear$best_model >= 1 & data_length_clear$best_model <= 9, "without_shift", "with_shift")


TS_with_shift <- data_length_clear$nbr_datapoint[data_length_clear$model_type == "with_shift"]  # Models with shift
TS_without_shift <- data_length_clear$nbr_datapoint[data_length_clear$model_type == "without_shift"]  # Models without shift

# T test
t_test_result <- t.test(TS_with_shift, TS_without_shift)
print(t_test_result)

# anova test
data_for_anova <- data.frame(
  nbr_datapoint = data_length_clear$nbr_datapoint,
  model_type = data_length_clear$model_type
)

anova_result <- aov(nbr_datapoint ~ model_type, data = data_for_anova)
summary.lm(anova_result)

# save the test results
sink(file = "./results_paleoTS_v0.6.1/correlation_length_modelfit_with_shift.txt")
paste("Comparison of datapoint number in relation to the best fit being a model with or without a shift.")
paste("Tests performed on the", length(data_length_clear$tsID), "time series used in the analysis with shift models.")
cat(rep("-", 60), "\n\n")
cat("Results of the T test:\n")
print(t_test_result)
cat(rep("-", 60), "\n\n")
cat("Anova test based on the minimum gap between two Aiccs of each time series and their adequacy state:\n")
print(summary.lm(anova_result))
sink()

# boxplot
png("./results_paleoTS_v0.6.1/plot/correlation_length_modelfit_with_shift_boxplot.png")
ggplot(data_for_anova, aes(x = model_type, y = nbr_datapoint, fill = model_type)) +
  geom_boxplot() +
  labs(title = "Boxplot of Number of Data Points by Model Type",
       x = "Model Type",
       y = "Number of Data Points") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()


#---------------------------------------------------
# Plot of the model fits depending on environment
#---------------------------------------------------

# Get all the time series which are microfossils, and add column of environment to the dataset
metadatashort2 <- metadatashort[metadatashort$tsID %in% metadata_envi$tsID, ]
metadata_envi$tsID <- as.character(metadata_envi$tsID)
metadatashort2$tsID <- as.character(metadatashort2$tsID)
metadatashort2$marine_envi <- left_join(metadatashort2, metadata_envi, by = "tsID")$Marine_environment

# extract the support results
model_test_envi <- model_test[names(model_test) %in% metadatashort2$tsID]

metadatashort2$support_Stasis_envi <- lapply(model_test_envi, function(df) {
  return(df["Stasis", "Akaike.wt"])
})
metadatashort2$support_Stasis_envi <- ifelse(metadatashort2$support_Stasis_envi > 0, 
                                             metadatashort2$support_Stasis_envi, NA)

metadatashort2$support_URW_envi <- lapply(model_test_envi, function(df) {
  return(df["URW", "Akaike.wt"])
})
metadatashort2$support_URW_envi <- ifelse(metadatashort2$support_URW_envi > 0, 
                                          metadatashort2$support_URW_envi, NA)

metadatashort2$support_GRW_envi <- lapply(model_test_envi, function(df) {
  return(df["GRW", "Akaike.wt"])
})
metadatashort2$support_GRW_envi <- ifelse(metadatashort2$support_GRW_envi > 0, 
                                          metadatashort2$support_GRW_envi, NA)

metadatashort2$support_Strict_stasis_envi <- lapply(model_test_envi, function(df) {
  return(df["StrictStasis", "Akaike.wt"])
})
metadatashort2$support_Strict_stasis_envi <- ifelse(metadatashort2$support_Strict_stasis_envi > 0, 
                                                    metadatashort2$support_Strict_stasis_envi, NA)

metadatashort2$support_Decel_envi <- lapply(model_test_envi, function(df) {
  return(df["Decel", "Akaike.wt"])
})
metadatashort2$support_Decel_envi <- ifelse(metadatashort2$support_Decel_envi > 0, 
                                            metadatashort2$support_Decel_envi, NA)

metadatashort2$support_Accel_envi <- lapply(model_test_envi, function(df) {
  return(df["Accel", "Akaike.wt"])
})
metadatashort2$support_Accel_envi <- ifelse(metadatashort2$support_Accel_envi > 0, 
                                            metadatashort2$support_Accel_envi, NA)

metadatashort2$support_OU_envi <- lapply(model_test_envi, function(df) {
  return(df["OU", "Akaike.wt"])
})
metadatashort2$support_OU_envi <- ifelse(metadatashort2$support_OU_envi > 0, 
                                         metadatashort2$support_OU_envi, NA)

metadatashort2$support_OU_mov_opt_anc_envi <- lapply(model_test_envi, function(df) {
  return(df["OU model with moving optimum (ancestral state at optimum)", "Akaike.wt"])
})
metadatashort2$support_OU_mov_opt_anc_envi <- ifelse(metadatashort2$support_OU_mov_opt_anc_envi > 0, 
                                                     metadatashort2$support_OU_mov_opt_anc_envi, NA)

metadatashort2$support_OU_mov_opt_envi <- lapply(model_test_envi, function(df) {
  return(df["OU model with moving optimum", "Akaike.wt"])
})
metadatashort2$support_OU_mov_opt_envi <- ifelse(metadatashort2$support_OU_mov_opt_envi > 0, 
                                                 metadatashort2$support_OU_mov_opt_envi, NA)

# convert the supports into numeric
metadatashort2$support_Stasis_envi <- as.numeric(as.character(metadatashort2$support_Stasis_envi))

# remove data with marine envi = unsure
metadatashort2_filtered <- metadatashort2 %>%
  filter(!marine_envi %in% c("unsure", "unsure benthic"))

# remove data with trait category = number
metadatashort2_filtered <- metadatashort2_filtered %>%
  filter(!trait_category %in% "number")

# classify trait_category in shape or size 
metadatashort2_filtered <- metadatashort2_filtered %>%
  mutate(size_or_shape = case_when(
    trait_category %in% c("size", "area") ~ "size",    # Assign "size" if ll is "sh"
    trait_category %in% c("shape", "angle") ~ "shape",))

metadatashort2_filtered <- metadatashort2_filtered %>%
  filter(!is.na(size_or_shape))

# Create the Boxplot
box_stasis = ggplot(metadatashort2_filtered, aes(x = interaction(marine_envi, size_or_shape), y = support_Stasis_envi, fill = marine_envi)) +
  geom_boxplot() +
  labs(title = "Comparison of Size and Shape Traits",
       x = "Trait (Species Type - Trait Type)",
       y = "Support for stasis") +
  theme_minimal()

# save the graphs
png("./results_paleoTS_v0.6.1/plot/results_support_envi_trait_boxplot.png", width = 2400, height = 1800)
grid.arrange(box_stasis, ncol = 1)
dev.off()

