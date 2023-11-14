#########################################
## Adequacy of models of evolution     ##
#########################################

#evoTS GitHub version
#adePEM new models version
#paleoTS verison 0.5.3

rm(list = ls())

library(foreach)
library(iterators)
library(parallel)
library(doParallel)

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

sink(file = "./results/AICc_results.txt")
aicc_results
paste("Total count:", length(ln_data))
paste("Percentage not URW, GRW or stasis:", percent2)
paste("Percentage not URW, GRW, stasis or strict stasis:", percent3)
sink()

###################
## Test adequacy ##
###################

# filter time series according to best AICc
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

# remove problem time series
OU <- OU[names(OU) != 427]
OU <- OU[names(OU) != 428]
OU <- OU[names(OU) != 584]
OU <- OU[names(OU) != 585]

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


# test adequacy
GRW_adeq <- mclapply(GRW, fit3adequacy.trend, plot = FALSE)
URW_adeq <- mclapply(URW, fit3adequacy.RW, plot = FALSE)
stasis_adeq <- mclapply(stasis, fit4adequacy.stasis, plot = FALSE) 
strict_stasis_adeq <- mclapply(strict_stasis, fit4adequacy.stasis, plot = FALSE)
decel_adeq <- mclapply(decel, fit3adequacy.decel, plot = FALSE)
accel_adeq <- mclapply(accel, fit3adequacy.RW, plot = FALSE)
#OU_adeq <- mclapply(OU, fit3adequacy.OU, plot = FALSE)
#save(file = "OU_adeq.Rdata", OU_adeq)
load("OU_adeq.Rdata")
#OU_mov_opt_anc_adeq <- mclapply(OU_mov_opt_anc, fit3adequacy.OU, plot = FALSE)
#save(file = "OU_mov_opt_anc_adeq.Rdata", OU_mov_opt_anc_adeq)
load("OU_mov_opt_anc_adeq.Rdata")
#OU_mov_opt_adeq <- mclapply(OU_mov_opt, fit3adequacy.OU, plot = FALSE)
save(file = "OU_mov_opt_adeq.Rdata", OU_mov_opt_adeq)
load("OU_mov_opt_adeq.Rdata")

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
adeq_table <- as.data.frame(c("GRW", "URW", "stasis", "strict stasis", "decel", "accel", "OU",
                "OU mov. optm. (ancestral state)", "OU mov. optm."))
colnames(adeq_table) <- "model"
adeq_table$count_passed <- c(GRW_c,URW_c, stasis_c, strict_stasis_c, decel_c, accel_c, OU_c, OU_mov_opt_anc_c, OU_mov_opt_c)
adeq_table$percentage_passed <- c(GRW_p,URW_p, stasis_p, strict_stasis_p, decel_p, accel_p, OU_p, OU_mov_opt_anc_p, OU_mov_opt_p)


# write to file
sink(file = "./results/adequacy_passed.txt")
adeq_table
sink()
