#########################################
## Evolutionary rates and time scaling ##
#########################################

#paleoTS.v.0.5.3
#evoTS GitHub version

rm(list = ls())

library(parallel)
library(doParallel)
library(adePEM)
library(paleoTS)

#source("/Users/vildeki/GitHub/assessing_models/assessing_models_functions.R")

# set working directory
setwd("/Users/kjetillv/GitHub/assessing_models/")

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

load(file = "OU_test.Rdata")
load(file = "OU_mov_opt_anc_test.Rdata")
load(file = "OU_mov_opt_test.Rdata")

###################
## Test adequacy ##
###################

OU_adeq <- lapply(OU_test[1:2], adeq_OU, plot = FALSE) # function added manually
OU_mov_opt_anc_adeq <- mclapply(OU_mov_opt_anc_test, adeq_OU, plot = FALSE) # function added manually
OU_mov_opt_adeq <- mclapply(OU_mov_opt_test, adeq_OU, plot = FALSE) # function added manually