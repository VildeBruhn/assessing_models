#########################################
## Evolutionary rates and time scaling ##
#########################################

#paleoTS.v.0.5.3
#evoTS GitHub version

rm(list = ls())

library(adePEM)
library(evoTS)
library(paleoTS)

# set working directory
setwd("/Users/kjetillv/GitHub/assessing_models/")

#-----------------
# IMPORT FILES
#-----------------

load(file = "OU_test.Rdata")

###################
## Test adequacy ##
###################

OU_adeq <- lapply(OU_test, fit3adequacy.OU, plot = FALSE)
