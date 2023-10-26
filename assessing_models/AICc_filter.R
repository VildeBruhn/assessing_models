rm(list = ls())

library(data.table)

# set working directory
setwd("/Users/vildeki/GitHub/assessing_models/")

#####################################
## Fit models and find best (AICc) ##
#####################################

# test all possible univariate models from evoTS on time series
#model_test <- mclapply(ln_data, fit.all.univariate, pool = TRUE)
load("model_test.Rdata")

# extract AICc values on all results
aicc <- lapply(model_test, function(x) x[(names(x) %in% c("AICc"))])

# get best AICc values within range
aicc_count <- lapply(aicc, function(x) {
  min_aicc <- min(x$AICc)
  within_range <- between(x$AICc, min_aicc, min_aicc + 2)
  sum(within_range)
})

# filter out time series that only have one AICc value within range
aicc_count <- unlist(aicc_count)
aicc_filter <- aicc_count[aicc_count == 1]

# get results compared to before filtering
sink(file = "./results/AICc_filter.txt")
paste("Count time series before filtering:", length(aicc_count))
paste("Count time series after filtering:", length(aicc_filter))
paste("Percentage time series after filtering:", (length(aicc_filter)/length(aicc_count))*100)
sink()