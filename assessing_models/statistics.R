#########################################
## Adequacy of models of evolution     ##
##          statistics                 ##
#########################################

#evoTS version 1.0.3
#adePEM version 1.1.1
#paleoTS verison 0.6.1
library(tidyverse)
library(devtools)
library(wesanderson)
rm(list = ls())

# set working directory
setwd("/Users/vildeki/GitHub/assessing_models/")
source("/Users/vildeki/GitHub/assessing_models/assessing_models_uni_functions.R")



#-----------------
# IMPORT FILES
#-----------------

# import
timeseries <- read_delim("./timeseries/timeseries.txt", col_names = TRUE, delim = "\t")
metadata <- read_delim("./timeseries/metadata.txt", col_names = TRUE, delim = "\t")
envi <- read_delim("./timeseries/Microfossils_envi.txt", col_names = TRUE, delim = "\t")

# join dataframes
df <- left_join(timeseries, metadata, by = c("tsID"))
envi <- select(envi, `Marine environment`, Environment, tsID)
df <- left_join(df, envi, by = c("tsID"))

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
load("./aicc_uni_passed.Rdata")

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

# bind data to dataframe
unit_list <- c("total_N", "steps", "interval_MY", "trait_type", "microfossil",
               "lat", "lon", "sediment", "model_aicc", "model_adequate", "Environment", "Marine environment")
bind_models <- bind(ln_data_meta, unit_list)

#library(nnet)
#model <- multinom(model_aicc ~ `Marine environment`, bind_models)
#summary(model)

# set colors
col_val <- c(wes_palette("Chevalier1"), wes_palette("IsleofDogs1")[6])

# aicc
models_marin_env <- bind_models
models_marin_env <- models_marin_env %>% drop_na(`Marine environment`)
pdf("/Users/vildeki/Downloads/marin_env_aicc.pdf")
ggplot(models_marin_env, aes(model_aicc, fill = `Marine environment`)) + geom_bar() +
  scale_fill_manual(values = col_val) + theme_classic()
dev.off()

models_env <- bind_models
models_env <- models_env %>% drop_na(Environment)
ggplot(models_env, aes(model_aicc, fill = Environment)) + geom_bar() +
  scale_fill_manual(values = wes_palette("Chevalier1")) + theme_classic()

ggplot(bind_models, aes(interval_MY, model_aicc)) + geom_boxplot() + theme_classic()

ggplot(bind_models, aes(total_N, model_aicc)) + geom_boxplot() + theme_classic()

ggplot(bind_models, aes(steps, model_aicc)) + geom_boxplot() + theme_classic()


# adequate
bind_models2 <- bind_models
bind_models2 <- bind_models2 %>% drop_na(model_adequate)

models_marin_env2 <- bind_models2
models_marin_env2 <- models_marin_env2 %>% drop_na(`Marine environment`)
pdf("/Users/vildeki/Downloads/marin_env_adeq.pdf")
ggplot(models_marin_env2, aes(model_adequate, fill = `Marine environment`)) + geom_bar() +
  scale_fill_manual(values = wes_palette("Chevalier1")) + theme_classic()
dev.off()

models_env2 <- bind_models2
models_env2 <- models_env2 %>% drop_na(Environment)
ggplot(models_env2, aes(model_adequate, fill = Environment)) + geom_bar() +
  scale_fill_manual(values = wes_palette("Chevalier1")) + theme_classic()

ggplot(models_env2, aes(model_adequate, fill = trait_type)) + geom_bar()

ggplot(bind_models2, aes(interval_MY, model_adequate)) + geom_boxplot() + theme_classic()

ggplot(bind_models2, aes(total_N, model_adequate)) + geom_boxplot() + theme_classic()

ggplot(bind_models2, aes(steps, model_adequate)) + geom_boxplot() + theme_classic()


#model2 <- multinom(model_adequate ~ `Marine environment`, models_marin_env2)
#summary(model2)
