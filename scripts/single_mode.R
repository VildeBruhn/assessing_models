####################################################
##    Assessing models of phenotypic evolution    ##
##                                                ##
##                 SINGLE-MODE                    ##
##                                                ##
##            Vilde Bruhn Kinneberg               ##
##            Created 01.06.2026                  ##
####################################################


# install packages
#install.packages("evoTS") #version 1.0.3
#install.packages("devtools")
#devtools::install_github("klvoje/adePEM") #version 1.1.1
#install.packages("tidyverse")
#install.packages("paleoTS") #version 0.6.1

# clean environment
rm(list = ls())

# load packages
library(evoTS)
library(adePEM)
library(tidyverse)
library(paleoTS)1

# set working directory
PATH = "[PATH_TO_MAIN_FOLDER]"
setwd(PATH)

# import functions
source("./scripts/functions.R")


#-------------------------
# IMPORT AND PROCESS FILES
#-------------------------


# import
timeseries <- read_delim("./data/timeseries/timeseries.txt", col_names = TRUE, delim = "\t")
metadata <- read_delim("./data/timeseries/metadata.txt", col_names = TRUE, delim = "\t")

# join data frames
df <- left_join(timeseries, metadata, by = c("tsID"))

# remove time series with less than 10 steps
df <- subset(df, steps >= 7)

# remove modern time series
df <- subset(df, period_start != "Present")

# remove Syverson
df <- subset(df, URL != "https://doi.org/10.1017/pab.2024.37")

# make list based on ID
df <- lapply(split(df,df$tsID), function(x) as.list(x))

# process data (warnings ok)
ln_data_meta <- dt(df, "tsID")
ln_data <- lapply(ln_data_meta, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$N, tt = x$tt, oldest = "first")
})

# Convert the time vector to unit length
ln_data_meta <- lapply(ln_data_meta, function(x) {
  x$tt <- x$tt/(max(x$tt))
  x
})
ln_data <- lapply(ln_data, function(x) {
  x$tt <- x$tt/(max(x$tt))
  x
})


#-----------------------
# FIT SINGLE-MODE MODELS
#-----------------------

# example of how to run model test (takes time, load model test used
# in article below)
#model_test <- list()
#for(i in 1:length(ln_data)){
#  try(model_test[[i]] <- fit.all.univariate(ln_data[[i]]))
#}
## add time series IDs
#names_list <- names(ln_data)
#names(model_test) <- names_list
## remove time series that cannot be processed by the loglikelihood function
#ln_data = ln_data[-which(sapply(model_test, is.null))]
#ln_data_meta = ln_data_meta[-which(sapply(model_test, is.null))]
#model_test = model_test[-which(sapply(model_test, is.null))]


## load model test and data used in article
load("./data/model_test_single_mode.Rdata")
load("./data/ln_data_single_mode.Rdata")
load("./data/ln_data_meta_single_mode.Rdata")


#----------------------------
# ASSESS RELATIVE FIT (AICc)
#----------------------------


# extract AICc values on all results
aicc <- lapply(model_test, function(x) x[(names(x) %in% c("AICc"))])

# check which AICc value is the lowest
aicc <- lapply(aicc, function(x) {
  which.min(as.numeric(unlist(x)))
})


#----------------------------
# MAKE RELATIVE FIT DATA SET
#----------------------------


# filter time series according to best AICc (or load data used in article below)
data_aicc <- mapply(c, ln_data, aicc, SIMPLIFY = FALSE) 
  #adds index of lowest AICc as column in ln_data

GRW <- Filter(function(x) x[[10]] == 1, data_aicc)
GRW <- lapply(GRW, function(x) { x[[10]] <- NULL; x })
GRW <- lapply(GRW, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt)
})

URW <- Filter(function(x) x[[10]] == 2, data_aicc)
URW <- lapply(URW, function(x) { x[[10]] <- NULL; x })
URW <- lapply(URW, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt)
})

stasis <- Filter(function(x) x[[10]] == 3, data_aicc)
stasis <- lapply(stasis, function(x) { x[[10]] <- NULL; x })
stasis<- lapply(stasis, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt)
})

strict_stasis <- Filter(function(x) x[[10]] == 4, data_aicc)
strict_stasis <- lapply(strict_stasis, function(x) { x[[10]] <- NULL; x })
strict_stasis <- lapply(strict_stasis, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt)
})

decel <- Filter(function(x) x[[10]] == 5, data_aicc)
decel <- lapply(decel, function(x) { x[[10]] <- NULL; x })
decel <- lapply(decel, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt)
})

accel <- Filter(function(x) x[[10]] == 6, data_aicc)
accel <- lapply(accel, function(x) { x[[10]] <- NULL; x })
accel <- lapply(accel, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt)
})

OU <- Filter(function(x) x[[10]] == 7, data_aicc)
OU <- lapply(OU, function(x) { x[[10]] <- NULL; x })
OU <- lapply(OU, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt)
})

OU_mov_opt_anc <- Filter(function(x) x[[10]] == 8, data_aicc)
OU_mov_opt_anc <- lapply(OU_mov_opt_anc, function(x) { x[[10]] <- NULL; x })
OU_mov_opt_anc<- lapply(OU_mov_opt_anc, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt)
})

OU_mov_opt <- Filter(function(x) x[[10]] == 9, data_aicc)
OU_mov_opt <- lapply(OU_mov_opt, function(x) { x[[10]] <- NULL; x })
OU_mov_opt <- lapply(OU_mov_opt, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt)
})


## load data used in article
load("./data/aicc_single_mode_passed.Rdata")


#--------------------------------
# ASSESS ABSOLUTE FIT (ADEQUACY) 
#--------------------------------


# Test adequacy

# GRW
GRW_adeq <- lapply(GRW, fit3adequacy.trend, plot = FALSE)

# URW (warnings ok)
URW_adeq <- lapply(URW, fit3adequacy.RW, plot = FALSE)

# stasis
stasis_adeq <- lapply(stasis, fit4adequacy.stasis, plot = FALSE) 

# strict stasis
strict_stasis_adeq <- lapply(strict_stasis, fit4adequacy.stasis, plot = FALSE)

# decelerated
decel_adeq <- lapply(decel, fit3adequacy.decel, plot = FALSE)

# accelerated
# reverse accelerated to become decelerated
accel_decel <- accel

accel_decel <- lapply(accel_decel, function(x) {
  x$mm <- rev(x$mm)
  x$vv <- rev(x$vv)
  x$nn <- rev(x$nn)
  x$tt <- rev(x$tt)
  for (i in 1:length(x$tt)){
    x$tt[i] <- 1 - x$tt[i]
  }
  return(x)
})

accel_adeq <- lapply(accel_decel, fit3adequacy.decel, plot = FALSE)

# example of how to test adequacy for OU models (this will take time, 
# load tests used in article below)
#OU_adeq <- list()
#for(i in 1:length(OU)){
#  OU_adeq[[i]] <- tryCatch({
#    fit3adequacy.OU(OU[[i]], plot = FALSE)
#  }, error = function(x) {
#    cat("An error occurred:", x$message, "\n")
#    return(OU_adeq[[i]] <- NA)
#  })
#}
## set names
#names(OU_adeq) <- names(OU)
## remove time series that cannot be processed by the loglikelihood function
#OU = Filter(function(x) length(x) > 1, OU_adeq)
#OU_adeq <- Filter(function(x) length(x) > 1, OU_adeq)

## same approach for OU_mov_opt and OU_mov_opt_anc

# load OU data run on HPC cluster
load("./data/OU_single_mode.Rdata")
load("./data/OU_mov_opt_single_mode.Rdata")
load("./data/OU_mov_opt_anc_single_mode.Rdata")
load("./data/OU_single_mode_adeq.Rdata")
load("./data/OU_mov_opt_single_mode_adeq.Rdata")
load("./data/OU_mov_opt_anc_single_mode_adeq.Rdata")


#----------------------------
# MAKE ABSOLUTE FIT DATA SET
#----------------------------


# get only adequate time series (or load data used in article below)
GRW_adeq_passed <- adequate3tests(GRW_adeq)
URW_adeq_passed <- adequate3tests(URW_adeq)
stasis_adeq_passed <- adequate4tests(stasis_adeq)
strict_stasis_adeq_passed <- adequate4tests(strict_stasis_adeq)
decel_adeq_passed <- adequate3tests(decel_adeq)
accel_adeq_passed <- adequate3tests(accel_adeq)
OU_adeq_passed <- adequate2tests(OU_adeq)
OU_mov_opt_anc_adeq_passed <- adequate2tests(OU_mov_opt_anc_adeq)
OU_mov_opt_adeq_passed <- adequate2tests(OU_mov_opt_adeq)

## load all adequacy data used in article
load("./data/adeq_single_mode_passed.Rdata")


