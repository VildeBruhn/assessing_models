#########################################
## Evolutionary rates and time scaling ##
#########################################

#paleoTS.v.0.5.3
#evoTS GitHub version

rm(list = ls())

library(doParallel)
library(adePEM)

source("/Users/vildeki/GitHub/assessing_models/assessing_models_functions.R")

# set working directory for database
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

# register it to be used by %dopar%
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

#####################################
## Fit models and find best (AICc) ##
#####################################

# test all possible univariate models from evoTS on timeseries
model_test <- mclapply(ln_data, fit.all.univariate, pool = TRUE)

# extract AICc values from URW on all results
aicc <- lapply(model_test, function(x) x[(names(x) %in% c("AICc"))])

# check which AICc value is the lowest
aicc <- lapply(aicc, function(x) {
  which.min(as.numeric(unlist(x)))
})

# get percentage
aicc_unlist <- unlist(aicc)
counts <- table(aicc_unlist)
percent <- (counts/sum(counts))*100
names(percent) <- c("GRW", "URW", "Stasis", "Strict stasis", "Decel", "Accel", "OU",
                "OU mov. optm. (ancestral state)", "OU mov. optm.")

# write to file
sink(file = "./results/percent_AICc.txt")
percent
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
strict_stasis_adeq <- mclapply(strict_stasis, adeq_stasis, plot = FALSE)
decel_adeq <- mclapply(decel, fit3adequacy.decel, plot = FALSE)
accel_adeq <- mclapply(accel, fit3adequacy.RW, plot = FALSE)
OU_adeq <- mclapply(OU, adeq_OU, plot = FALSE)
OU_mov_opt_anc_adeq <- mclapply(OU_mov_opt_anc, adeq_OU, plot = FALSE)
OU_mov_opt_adeq <- mclapply(OU_mov_opt, adeq_OU, plot = FALSE)

# Get only adequate time series
GRW_adeq_passed <- adequate3tests(GRW_adeq)
URW_adeq_passed <- adequate3tests(URW_adeq)
stasis_adeq_passed <- adequate4tests(stasis_adeq)
strict_stasis_adeq_passed <- adequate4tests(strict_stasis_adeq)
decel_adeq_passed <- adequate3tests(decel_adeq)
accel_adeq_passed <- adequate3tests(accel_adeq)
OU_adeq_passed <- adequate3tests(OU_adeq)
OU_mov_opt_anc_adeq_passed <- adequate3tests(OU_mov_opt_anc_adeq)
OU_mov_opt_adeq_passed <- adequate3tests(OU_mov_opt_adeq)

