##########################################################################
## Evolutionary rates and time series scaling with shift model included ##
##########################################################################

#paleoTS.v.0.5.3
#evoTS GitHub version

rm(list = ls())

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

# make list based on ID
df <- lapply(split(df,df$tsID), function(x) as.list(x))

# process data
ln_data_meta <- dt(df, "tsID")
ln_data <- lapply(ln_data_meta, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$N, tt = x$tt, oldest = "first")
})
             
# exclude too short timeseries            




             
#####################################################
## Fit models including shift and find best (AICc) ##
#####################################################

#-----------------------------------------------------------------------------------------------
# need to define the minimum of samples that can be considered to be a segment of the model
#------------------------------------------------------------------------------------------------


















             
