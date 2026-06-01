##########################################################################
## Evolutionary rates and time series scaling with shift model included ##
##########################################################################

#paleoTS.v.0.5.3
#evoTS GitHub version

.libPaths("/cluster/home/marionth/R")


rm(list = ls())

library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(adePEM)
library(evoTS)
library(paleoTS)
library(gridExtra)
library(ggplot2)
library(tidyverse)

source("/cluster/home/marionth/Project6/assessing_models/assessing_models_uni_functions.R")

# set working directory
setwd("/cluster/home/marionth/Project6/assessing_models/")

# -------------------------
# Set up for parallel runs
# -------------------------


#number of cores that I want to use (also defined in the batch script)
registerDoParallel(52)

#To calculate the time needed to process each iteration
start_time <- Sys.time()


#---------------------------------------
# Testing the adequacy of the OU models
#---------------------------------------

#R version 4.2.1
#evoTS version 1.0.3
#adePEM new models version
#paleoTS version 0.6.2

load('OU_shift.RData')


OU_adeq <- vector("list", length(OU))
if (length(OU) > 0) { 
  for (i in seq_along(OU)) {
    OU_adeq[[i]] <- tryCatch(
      fit3adequacy.OU(OU[[i]], plot = FALSE),
      error = function(e) NA
    )
  }
}
names(OU_adeq) <- names(OU)
adeq_issues_OU <- setdiff(names(OU), names(OU_adeq))
OU_adeq <- Filter(function(x) !is.na(x)[1], OU_adeq)

OU_mov_opt_anc_adeq <- vector("list", length(OU_mov_opt_anc))
if (length(OU_mov_opt_anc) > 0) { 
  for (i in seq_along(OU_mov_opt_anc)) {
    OU_mov_opt_anc_adeq[[i]] <- tryCatch(
      fit3adequacy.OU(OU_mov_opt_anc[[i]], plot = FALSE),
      error = function(e) NA
    )
  }
}
names(OU_mov_opt_anc_adeq) <- names(OU_mov_opt_anc)
adeq_issues_OU_mov_opt_anc <- which(sapply(OU_mov_opt_anc_adeq, function(x) is.na(x)[1]))
adeq_issues_OU = c(adeq_issues_OU, adeq_issues_OU_mov_opt_anc)
OU_mov_opt_anc_adeq <- Filter(function(x) !is.na(x)[1], OU_mov_opt_anc_adeq)

OU_mov_opt_adeq <- vector("list", length(OU_mov_opt))
if (length(OU_mov_opt) > 0) { 
  for (i in seq_along(OU_mov_opt)) {
    OU_mov_opt_adeq[[i]] <- tryCatch(
      fit3adequacy.OU(OU_mov_opt[[i]], plot = FALSE),
      error = function(e) NA
    )
  }
}
names(OU_mov_opt_adeq) <- names(OU_mov_opt)
adeq_issues_OU_mov_opt <- which(sapply(OU_mov_opt_adeq, function(x) is.na(x)[1]))
adeq_issues_OU = c(adeq_issues_OU, adeq_issues_OU_mov_opt)
OU_mov_opt_adeq <- Filter(function(x) !is.na(x)[1], OU_mov_opt_adeq)

Stasis_OU_subset2_adeq <- vector("list", length(Stasis_OU_subset2))
if (length(Stasis_OU_subset2) > 0) { 
  for (i in seq_along(Stasis_OU_subset2)) {
    Stasis_OU_subset2_adeq[[i]] <- tryCatch(
      fit3adequacy.OU(Stasis_OU_subset2[[i]], plot = FALSE),
      error = function(e) NA
    )
  }
}
names(Stasis_OU_subset2_adeq) <- names(Stasis_OU_subset2)
adeq_issues_Stasis_OU_subset2 <- which(sapply(Stasis_OU_subset2_adeq, function(x) is.na(x)[1]))
adeq_issues_OU = c(adeq_issues_OU, adeq_issues_Stasis_OU_subset2)
Stasis_OU_subset2_adeq <- Filter(function(x) !is.na(x)[1], Stasis_OU_subset2_adeq)

URW_OU_subset2_adeq <- vector("list", length(URW_OU_subset2))
if (length(URW_OU_subset2) > 0) { 
  for (i in seq_along(URW_OU_subset2)) {
    URW_OU_subset2_adeq[[i]] <- tryCatch(
      fit3adequacy.OU(URW_OU_subset2[[i]], plot = FALSE),
      error = function(e) NA
    )
  }
}
names(URW_OU_subset2_adeq) <- names(URW_OU_subset2)
adeq_issues_URW_OU_subset2 <- which(sapply(URW_OU_subset2_adeq, function(x) is.na(x)[1]))
adeq_issues_OU = c(adeq_issues_OU, adeq_issues_URW_OU_subset2)
URW_OU_subset2_adeq <- Filter(function(x) !is.na(x)[1], URW_OU_subset2_adeq)

GRW_OU_subset2_adeq <- vector("list", length(GRW_OU_subset2))
if (length(GRW_OU_subset2) > 0) { 
  for (i in seq_along(GRW_OU_subset2)) {
    GRW_OU_subset2_adeq[[i]] <- tryCatch(
      fit3adequacy.OU(GRW_OU_subset2[[i]], plot = FALSE),
      error = function(e) NA
    )
  }
}
names(GRW_OU_subset2_adeq) <- names(GRW_OU_subset2)
adeq_issues_GRW_OU_subset2 <- which(sapply(GRW_OU_subset2_adeq, function(x) is.na(x)[1]))
adeq_issues_OU = c(adeq_issues_OU, adeq_issues_GRW_OU_subset2)
GRW_OU_subset2_adeq <- Filter(function(x) !is.na(x)[1], GRW_OU_subset2_adeq)

OU_OU_subset1_adeq <- vector("list", length(OU_OU_subset1))
if (length(OU_OU_subset1) > 0) { 
  for (i in seq_along(OU_OU_subset1)) {
    OU_OU_subset1_adeq[[i]] <- tryCatch(
      fit3adequacy.OU(OU_OU_subset1[[i]], plot = FALSE),
      error = function(e) NA
    )
  }
}
names(OU_OU_subset1_adeq) <- names(OU_OU_subset1)
adeq_issues_OU_OU_subset1 <- which(sapply(OU_OU_subset1_adeq, function(x) is.na(x)[1]))
adeq_issues_OU = c(adeq_issues_OU, adeq_issues_OU_OU_subset1)
OU_OU_subset1_adeq <- Filter(function(x) !is.na(x)[1], OU_OU_subset1_adeq)

OU_OU_subset2_adeq <- vector("list", length(OU_OU_subset2))
if (length(OU_OU_subset2) > 0) { 
  for (i in seq_along(OU_OU_subset2)) {
    OU_OU_subset2_adeq[[i]] <- tryCatch(
      fit3adequacy.OU(OU_OU_subset2[[i]], plot = FALSE),
      error = function(e) NA
    )
  }
}
names(OU_OU_subset2_adeq) <- names(OU_OU_subset2)
adeq_issues_OU_OU_subset2 <- which(sapply(OU_OU_subset2_adeq, function(x) is.na(x)[1]))
adeq_issues_OU = c(adeq_issues_OU, adeq_issues_OU_OU_subset2)
OU_OU_subset2_adeq <- Filter(function(x) !is.na(x)[1], OU_OU_subset2_adeq)

OU_GRW_subset1_adeq <- vector("list", length(OU_GRW_subset1))
if (length(OU_GRW_subset1) > 0) { 
  for (i in seq_along(OU_GRW_subset1)) {
    OU_GRW_subset1_adeq[[i]] <- tryCatch(
      fit3adequacy.OU(OU_GRW_subset1[[i]], plot = FALSE),
      error = function(e) NA
    )
  }
}
names(OU_GRW_subset1_adeq) <- names(OU_GRW_subset1)
adeq_issues_OU_GRW_subset1 <- which(sapply(OU_GRW_subset1_adeq, function(x) is.na(x)[1]))
adeq_issues_OU = c(adeq_issues_OU, adeq_issues_OU_GRW_subset1)
OU_GRW_subset1_adeq <- Filter(function(x) !is.na(x)[1], OU_GRW_subset1_adeq)

OU_URW_subset1_adeq <- vector("list", length(OU_URW_subset1))
if (length(OU_URW_subset1) > 0) { 
  for (i in seq_along(OU_URW_subset1)) {
    OU_URW_subset1_adeq[[i]] <- tryCatch(
      fit3adequacy.OU(OU_URW_subset1[[i]], plot = FALSE),
      error = function(e) NA
    )
  }
}
names(OU_URW_subset1_adeq) <- names(OU_URW_subset1)
adeq_issues_OU_URW_subset1 <- which(sapply(OU_URW_subset1_adeq, function(x) is.na(x)[1]))
adeq_issues_OU = c(adeq_issues_OU, adeq_issues_OU_URW_subset1)
OU_URW_subset1_adeq <- Filter(function(x) !is.na(x)[1], OU_URW_subset1_adeq)

OU_Stasis_subset1_adeq <- vector("list", length(OU_Stasis_subset1))
if (length(OU_Stasis_subset1) > 0) { 
  for (i in seq_along(OU_Stasis_subset1)) {
    OU_Stasis_subset1_adeq[[i]] <- tryCatch(
      fit3adequacy.OU(OU_Stasis_subset1[[i]], plot = FALSE),
      error = function(e) NA
    )
  }
}
names(OU_Stasis_subset1_adeq) <- names(OU_Stasis_subset1)
adeq_issues_OU_Stasis_subset1 <- which(sapply(OU_Stasis_subset1_adeq, function(x) is.na(x)[1]))
adeq_issues_OU = c(adeq_issues_OU, adeq_issues_OU_Stasis_subset1)
OU_Stasis_subset1_adeq <- Filter(function(x) !is.na(x)[1], OU_Stasis_subset1_adeq)

save(adeq_issues_OU, OU_adeq, OU_mov_opt_anc_adeq, OU_mov_opt_adeq, Stasis_OU_subset2_adeq, 
     URW_OU_subset2_adeq, GRW_OU_subset2_adeq, OU_OU_subset1_adeq, OU_OU_subset2_adeq, 
     OU_GRW_subset1_adeq, OU_URW_subset1_adeq, OU_Stasis_subset1_adeq, file="OU_shift_adeq.RData")
