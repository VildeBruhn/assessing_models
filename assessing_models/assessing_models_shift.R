##########################################################################
## Evolutionary rates and time series scaling with shift model included ##
##########################################################################

#paleoTS.v.0.5.3
#evoTS GitHub version

rm(list = ls())

library(parallel)
library(doParallel)
library(adePEM)
library(evoTS)
library(paleoTS)

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

#---------------------------------------------------------------------------------------------
# EXCLUDING THE SHORTEST TIMESERIES (Keep timeseries only if containing more than 10 steps)
#---------------------------------------------------------------------------------------------

metadatalong <- matrix(nrow = 0, ncol = ncol(metadata))
timeserieslong <- matrix(nrow = 0, ncol = ncol(timeseries))
tsIDremoved <- c()

# Filter rows of the metadata file based on step number 
for (i in 1:nrow(metadata)) {                          
  if (metadata[i, "steps"] >= 10) {   
    metadatalong <- rbind(metadatalong, metadata[i, ])  # Add the filtered row to metadatalong
  } else {
    tsIDremoved <- c(tsIDremoved, metadata[i, "tsID"]) # Save the ID of removed timeseries
  }
}

# Filter rows of the timeseries file based on tsID number 
for (j in 1:nrow(timeseries)) {
    tsID <- timeseries[j, "tsID"]
    matching_tsID <- tsID %in% tsIDremoved # Check if tsID exists in tsIDremoved
    if (!matching_tsID) {
      timeserieslong <- rbind(timeserieslong, timeseries[j, ])  # Add the filtered row to timeserieslong if tsID does not exist in tsIDremoved
    }
}

# join dataframes
dflong <- left_join(timeserieslong, metadatalong, by = c("tsID"))

# make list based on ID
dflong <- lapply(split(dflong,dflong$tsID), function(x) as.list(x))

# process data
ln_data_metalong <- dt(dflong, "tsID")
ln_datalong <- lapply(ln_data_metalong, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$N, tt = x$tt, oldest = "first")
})
                        
#####################################################
## Fit models including shift and find best (AICc) ##
#####################################################

# test all possible univariate models from evoTS on every timeseries
model_test_noshift <- mclapply(ln_datalong, fit.all.univariate, pool = TRUE)

# test all possible shift models from evoTS on time series
fit_mode_shift <- function(ln_datalong) {
  models_list <- c("Stasis", "URW", "GRW", "OU")
  store_results <- list()
  k <- 0
  for (i in 1:4) {
    model1 <- models_list[i]  
    for (j in 1:4) {
      model2 <- models_list[j]
        fit_result <- fit.mode.shift(ln_datalong, model1, model2, minb = 5)
        k = k + 1
        store_results[[k]] <- fit_result
    } 
  }
  return(store_results)
}

model_shift_results <- mclapply(ln_datalong, fit_mode_shift)

#rename the sublist with the names of the models
model_shift_results_sample <- lapply(model_shift_results_sample, function(x) {
    names(x) <- c("Stasis-Stasis","Stasis-URW", "Stasis-GRW", "Stasis-OU", "URW-URW", "URW-GRW", 
                  "URW-OU","GRW-GRW", "GRW-OU", "OU-OU", "OU-GRW", "OU-URW", "OU-Stasis", 
                  "GRW-URW","GRW-Stasis", "URW-Stasis")
    return(x)
})
                 
# Save the results
save.image(file <- 'results_fit_models.RData')

                 
##############################################################
# Combining results of no shift and shift in the same file
##############################################################









             
