#########################################
## Evolutionary rates and time scaling ##
#########################################

#evoTS GitHub version

rm(list = ls())

library(adePEM)
library(evoTS)
library(paleoTS)

# set working directory
setwd("/Users/kjetillv/GitHub/assessing_models/")
setwd("/Users/vildeki/GitHub/assessing_models/")
#-----------------
# IMPORT FILES
#-----------------

load(file = "OU_issue.Rdata")

###################
## Test adequacy ##
###################

OU_adeq <- lapply(issue, fit3adequacy.OU, plot = FALSE)
