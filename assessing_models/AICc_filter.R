#########################################
## Evolutionary rates and time scaling ##
#########################################

#evoTS GitHub version
#adePEM new models version

rm(list = ls())

library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(adePEM)

source("/Users/vildeki/GitHub/assessing_models/assessing_models_functions.R")

# set working directory
setwd("/Users/vildeki/GitHub/assessing_models/")

# -------------------------
# Set up for parallel runs
# -------------------------

n_cores <- parallel::detectCores() - 1

# create the cluster
my_cluster <- parallel::makeCluster(
  n_cores, 
  type = "FORK"
)

# register it to be used
doParallel::registerDoParallel(cl = my_cluster)

#-----------------
# IMPORT FILES
#-----------------

# import
timeseries <- read_delim("./timeseries/timeseries.txt", col_names = TRUE, delim = "\t")
metadata <- read_delim("./timeseries/metadata.txt", col_names = TRUE, delim = "\t")

# join dataframes
df <- left_join(timeseries, metadata, by = c("tsID"))

# remove time series with less than 10 steps
df <- subset(df, steps >= 10)

# make list based on ID
df <- lapply(split(df,df$tsID), function(x) as.list(x))

# process data
ln_data_meta <- dt(df, "tsID")
ln_data <- lapply(ln_data_meta, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$N, tt = x$tt, oldest = "first")
})

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