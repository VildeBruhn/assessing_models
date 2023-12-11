#########################################
## Adequacy of models of evolution     ##
##          statistics                 ##
#########################################

#evoTS GitHub version
#adePEM new models version
#paleoTS verison 0.5.3

rm(list = ls())

# set working directory
setwd("/Users/vildeki/GitHub/assessing_models/")
source("/Users/vildeki/GitHub/assessing_models/assessing_models_functions.R")

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

# load timeseries best according to aicc
load("./adeq_uni_aicc.Rdata")

# load adequate timeseries 
load("./adeq_uni_passed.Rdata")

# get aicc model info into metadata
ln_data_meta <- model_aicc(ln_data_meta, GRW, "GRW")
ln_data_meta <- model_aicc(ln_data_meta, URW, "URW")
ln_data_meta <- model_aicc(ln_data_meta, stasis, "stasis")
ln_data_meta <- model_aicc(ln_data_meta, strict_stasis, "strict stasis")
ln_data_meta <- model_aicc(ln_data_meta, decel, "decel")
ln_data_meta <- model_aicc(ln_data_meta, accel, "accel")
ln_data_meta <- model_aicc(ln_data_meta, OU, "OU")
ln_data_meta <- model_aicc(ln_data_meta, OU_mov_opt_anc, "OU mov opt anc")
ln_data_meta <- model_aicc(ln_data_meta, OU_mov_opt, "OU mov opt")

# get adequate model info into metadata
ln_data_meta <- model_adeq(ln_data_meta, GRW_adeq_passed, "GRW")
ln_data_meta <- model_adeq(ln_data_meta, URW_adeq_passed, "URW")
ln_data_meta <- model_adeq(ln_data_meta, stasis_adeq_passed, "stasis")
ln_data_meta <- model_adeq(ln_data_meta, strict_stasis_adeq_passed, "strict stasis")
ln_data_meta <- model_adeq(ln_data_meta, decel_adeq_passed, "decel")
ln_data_meta <- model_adeq(ln_data_meta, accel_adeq_passed, "accel")
ln_data_meta <- model_adeq(ln_data_meta, OU_adeq_passed, "OU")
ln_data_meta <- model_adeq(ln_data_meta, OU_mov_opt_anc_adeq_passed, "OU mov opt anc")
ln_data_meta <- model_adeq(ln_data_meta, OU_mov_opt_adeq_passed, "OU mov opt")









