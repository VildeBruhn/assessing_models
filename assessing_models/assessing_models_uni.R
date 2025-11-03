#########################################
##    Assessing models of evolution    ##
##          UNIVARIATE MODELS          ##
#########################################


# install packages
#install.packages("evoTS") #version 1.0.3
#install.packages("devtools")
#devtools::install_github("klvoje/adePEM")
#install.packages("tidyverse")
#install.packages("data.table")
#install.packages("paleoTS") #version 0.6.1

# clean environment
rm(list = ls())

# load packages
library(evoTS) # version 1.0.3
library(adePEM) # version 1.1.1
library(tidyverse)
library(data.table)
library(paleoTS) # version 0.6.1
library(doParallel)

# load functions
source("/Users/markusof/Dropbox/Familiemappe/Vilde_jobb/assessing_models/assessing_models_uni_functions.R")

# set working directory
setwd("/Users/markusof/Dropbox/Familiemappe/Vilde_jobb/assessing_models/")


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


#--------------
# IMPORT FILES
#--------------


# import
timeseries <- read_delim("./timeseries/timeseries.txt", col_names = TRUE, delim = "\t")
metadata <- read_delim("./timeseries/metadata.txt", col_names = TRUE, delim = "\t")

# join dataframes
df <- left_join(timeseries, metadata, by = c("tsID"))

# remove time series with less than 10 steps
df <- subset(df, steps >= 7)

# remove modern time series
df <- subset(df, period_start != "Present")

# make list based on ID
df <- lapply(split(df,df$tsID), function(x) as.list(x))

# process data
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
# FIT UNIVARIATE MODELS
#-----------------------

# example of how to run model test (takes time, load model test used
# in article below)
model_test <- list()
for(i in 1:length(ln_data)){
  try(model_test[[i]] <- fit.all.univariate(ln_data[[i]]))
}
# add time series IDs
names_list <- names(ln_data)
names(model_test) <- names_list
# remove time series that cannot be processed by the loglikelihood function
ln_data = ln_data[-which(sapply(model_test, is.null))]
ln_data_meta = ln_data_meta[-which(sapply(model_test, is.null))]
model_test = model_test[-which(sapply(model_test, is.null))]

# save data
#save(model_test, file = "./model_test_uni.Rdata")
#save(ln_data, file = "./ln_data_uni.Rdata")
#save(ln_data_meta, file = "./ln_data_meta_uni.Rdata")

## load model test and data used in article
load("./model_test_uni.Rdata")
load("./ln_data_uni.Rdata")
load("./ln_data_meta_uni.Rdata")


#----------------------------
# ASSESS RELATIVE FIT (AICc)
#----------------------------


# extract AICc values on all results
aicc <- lapply(model_test, function(x) x[(names(x) %in% c("AICc"))])

# check which AICc value is the lowest
aicc <- lapply(aicc, function(x) {
  which.min(as.numeric(unlist(x)))
})

# get percentage
aicc_unlist <- unlist(aicc)
aicc_results <- table(aicc_unlist)
names(aicc_results) <- c("GRW", "URW", "Stasis", "Strict stasis", "Decel", "Accel", "OU",
                   "OU mov. optm. (ancestral state)", "OU mov. optm.")
aicc_results <- as.data.frame(aicc_results)
colnames(aicc_results) <- c("model", "count")
aicc_results$percentage <- (aicc_results$count/sum(aicc_results$count))*100
percent2 <- sum(aicc_results$percentage[4:9])
percent3 <- sum(aicc_results$percentage[5:9])

sink(file = "./results_paleoTS_v0.6.1/AICc_uni_results.txt")
aicc_results
paste("Total count:", length(aicc))
paste("Percentage not URW, GRW or stasis:", percent2)
paste("Percentage not URW, GRW, stasis or strict stasis:", percent3)
sink()


#----------------------------
# MAKE RELATIVE FIT DATA SET
#----------------------------


# filter time series according to best AICc (load data used in article below)
data_aicc <- mapply(c, ln_data, aicc, SIMPLIFY = FALSE) #adds index of lowest AICc as column in ln_data

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

# save data
#save(GRW, URW, stasis,
#strict_stasis, decel,
#accel, OU, OU_mov_opt, 
#OU_mov_opt_anc, file = "aicc_uni_passed.Rdata")

## load data used in article
load("./aicc_uni_passed.Rdata")


#--------------------------------
# ASSESS ABSOLUTE FIT (ADEQUACY) 
#--------------------------------


# test adequacy (load data used in article below)

GRW_adeq <- list()
for(i in 1:length(GRW)){
  try(GRW_adeq[[i]] <- fit3adequacy.trend(GRW[[i]], plot = FALSE))
}
# add time series IDs
names_list <- names(GRW)
names(GRW_adeq) <- names_list
# remove time series that cannot be processed by the loglikelihood function
GRW = GRW[-which(sapply(GRW_adeq, is.null))]
GRW_adeq = GRW_adeq[-which(sapply(GRW_adeq, is.null))]

URW_adeq <- mclapply(URW, fit3adequacy.RW, plot = FALSE)

stasis_adeq <- mclapply(stasis, fit4adequacy.stasis, plot = FALSE) 

strict_stasis_adeq <- mclapply(strict_stasis, fit4adequacy.stasis, plot = FALSE)

decel_adeq <- mclapply(decel, fit3adequacy.decel, plot = FALSE)

accel_adeq <- mclapply(accel, fit3adequacy.RW, plot = FALSE)

# example of how to test adequacy for OU models (this will take time, 
# load tests used in article below)
OU_adeq <- list()
for(i in 1:length(OU)){
  OU_adeq[[i]] <- tryCatch({
    fit3adequacy.OU(OU[[i]], plot = FALSE)
  }, error = function(x) {
    cat("An error occurred:", x$message, "\n")
    return(OU_adeq[[i]] <- NA)
  })
}
# set names
names(OU_adeq) <- names(OU)
# remove time series that cannot be processed by the loglikelihood function
OU = Filter(function(x) length(x) > 1, OU_adeq)
OU_adeq <- Filter(function(x) length(x) > 1, OU_adeq)
# save
save(file = "OU_uni.Rdata", OU)
save(file = "OU_uni_adeq.Rdata", OU_adeq)


## same approach for OU_mov_opt and OU_mov_opt_anc

## load OU tests and data used in article
load("OU_uni_adeq.Rdata")
load("OU_mov_opt_uni_adeq.Rdata")
load("OU_mov_opt_anc_uni_adeq.Rdata")
load("OU_uni.Rdata")
load("OU_mov_opt_uni.Rdata")
load("OU_mov_opt_anc_uni.Rdata")

# get only adequate time series
GRW_adeq_passed <- adequate3tests(GRW_adeq)
URW_adeq_passed <- adequate3tests(URW_adeq)
stasis_adeq_passed <- adequate4tests(stasis_adeq)
strict_stasis_adeq_passed <- adequate4tests(strict_stasis_adeq)
decel_adeq_passed <- adequate3tests(decel_adeq)
accel_adeq_passed <- adequate3tests(accel_adeq)
OU_adeq_passed <- adequate2tests(OU_adeq)
OU_mov_opt_anc_adeq_passed <- adequate2tests(OU_mov_opt_anc_adeq)
OU_mov_opt_adeq_passed <- adequate2tests(OU_mov_opt_adeq)

# save data
#save(GRW_adeq_passed, URW_adeq_passed, stasis_adeq_passed,
#strict_stasis_adeq_passed, decel_adeq_passed,
#accel_adeq_passed, OU_adeq_passed, OU_mov_opt_adeq_passed, OU_mov_opt_anc_adeq_passed, 
#file = "adeq_uni_passed.Rdata")

## load all adequacy data used in article
load("./adeq_uni_passed.Rdata")

# get counts passed
GRW_c <- length(GRW_adeq_passed)
URW_c <- length(URW_adeq_passed)
stasis_c <- length(stasis_adeq_passed)
strict_stasis_c <- length(strict_stasis_adeq_passed)
decel_c <- length(decel_adeq_passed)
accel_c <- length(accel_adeq_passed)
OU_c <- length(OU_adeq_passed)
OU_mov_opt_anc_c <- length(OU_mov_opt_anc_adeq_passed)
OU_mov_opt_c <- length(OU_mov_opt_adeq_passed)

# get percentage passed
GRW_p <- (length(GRW_adeq_passed)/length(GRW_adeq))*100
URW_p <- (length(URW_adeq_passed)/length(URW_adeq))*100
stasis_p <- (length(stasis_adeq_passed)/length(stasis_adeq))*100
strict_stasis_p <- (length(strict_stasis_adeq_passed)/length(strict_stasis_adeq))*100
decel_p <- (length(decel_adeq_passed)/length(decel_adeq))*100
accel_p <- (length(accel_adeq_passed)/length(accel_adeq))*100
OU_p <- (length(OU_adeq_passed)/length(OU_adeq))*100
OU_mov_opt_anc_p <- (length(OU_mov_opt_anc_adeq_passed)/length(OU_mov_opt_anc_adeq))*100
OU_mov_opt_p <- (length(OU_mov_opt_adeq_passed)/length(OU_mov_opt_adeq))*100

# make output table
total_count <- sum(GRW_c,URW_c, stasis_c, strict_stasis_c, decel_c, accel_c, OU_c, OU_mov_opt_anc_c, OU_mov_opt_c)
adeq_table <- as.data.frame(c("GRW", "URW", "stasis", "strict stasis", "decel", "accel", "OU",
                "OU mov. optm. (ancestral state)", "OU mov. optm."))
colnames(adeq_table) <- "model"
adeq_table$count_passed <- c(GRW_c,URW_c, stasis_c, strict_stasis_c, decel_c, accel_c, OU_c, OU_mov_opt_anc_c, OU_mov_opt_c)
adeq_table$percentage_passed <- c(GRW_p,URW_p, stasis_p, strict_stasis_p, decel_p, accel_p, OU_p, OU_mov_opt_anc_p, OU_mov_opt_p)

# get percentage
percent2 <- (sum(adeq_table$count_passed[4:9])/total_count)*100
percent3 <- (sum(adeq_table$count_passed[5:9])/total_count)*100

# write to file
sink(file = "./results_paleoTS_v0.6.1/adequacy_uni_passed.txt")
adeq_table
paste("Total count:", total_count)
paste("Percentage not URW, GRW or stasis:", percent2)
paste("Percentage not URW, GRW, stasis or strict stasis:", percent3)
sink()




