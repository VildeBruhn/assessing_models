##########################################################################
## Evolutionary rates and time series scaling with shift model included ##
##########################################################################

#paleoTS.v.0.5.3
#evoTS GitHub version

rm(list = ls())

library(parallel)
library(doParallel)
library(adePEM)
library(evoTS)
library(paleoTS)

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

#---------------------------------------------------------------------------------------------
# EXCLUDING THE SHORTEST TIMESERIES (Keep timeseries only if containing more than 10 steps)
#---------------------------------------------------------------------------------------------

metadatalong <- matrix(nrow = 0, ncol = ncol(metadata))
timeserieslong <- matrix(nrow = 0, ncol = ncol(timeseries))
tsIDremoved <- c()

# Filter rows of the metadata file based on step number 
for (i in 1:nrow(metadata)) {                          
  if (metadata[i, "steps"] >= 10) {   
    metadatalong <- rbind(metadatalong, metadata[i, ])  # Add the filtered row to metadatalong
  } else {
    tsIDremoved <- c(tsIDremoved, metadata[i, "tsID"]) # Save the ID of removed timeseries
  }
}

# Filter rows of the timeseries file based on tsID number 
for (j in 1:nrow(timeseries)) {
    tsID <- timeseries[j, "tsID"]
    matching_tsID <- tsID %in% tsIDremoved # Check if tsID exists in tsIDremoved
    if (!matching_tsID) {
      timeserieslong <- rbind(timeserieslong, timeseries[j, ])  # Add the filtered row to timeserieslong if tsID does not exist in tsIDremoved
    }
}

# join dataframes
dflong <- left_join(timeserieslong, metadatalong, by = c("tsID"))

# make list based on ID
dflong <- lapply(split(dflong,dflong$tsID), function(x) as.list(x))

# process data
ln_data_metalong <- dt(dflong, "tsID")
ln_datalong <- lapply(ln_data_metalong, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$N, tt = x$tt, oldest = "first")
})
                        
#####################################################
## Fit models including shift and find best (AICc) ##
#####################################################

# test all possible univariate models from evoTS on every timeseries
model_test_noshift <- mclapply(ln_datalong, fit.all.univariate, pool = TRUE)

# test all possible shift models from evoTS on time series
fit_mode_shift <- function(ln_datalong) {
  models_list <- c("Stasis", "URW", "GRW", "OU")
  store_results <- list()
  k <- 0
  for (i in 1:4) {
    model1 <- models_list[i]  
    for (j in 1:4) {
      model2 <- models_list[j]
        fit_result <- fit.mode.shift(ln_datalong, model1, model2, minb = 5)
        k = k + 1
        store_results[[k]] <- fit_result
    } 
  }
  return(store_results)
}

model_shift_results <- mclapply(ln_datalong, fit_mode_shift)

########################
##  Extract the AICcs ##
########################

### Remove problematic timeseries ###
pblm_TS = c("567","575","576")
keep_TS <- !names(model_shift_results) %in% pblm_TS
model_shift_results_clean = model_shift_results[keep_TS]
keep_TS2 <- !names(model_noshift_results) %in% pblm_TS
model_noshift_results_clean = model_noshift_results[keep_TS2]

#------------------------------------------
# Extract AICcs for the no shift models
#------------------------------------------

aicc_noshift <- lapply(model_noshift_results_clean, function(x) x[(names(x) %in% c("AICc"))])

#------------------------------------------
# Extract AICcs for the shift models
#------------------------------------------

# extract AICc values of shift models on all results
aicc_shift_extraction <- lapply(model_shift_results_clean, function(x) { ### remettre model_shift_result when problem will be solved
  sapply(x, function(result) result$AICc)
})

model_names <- c(
  "Stasis-Stasis","Stasis-URW", "Stasis-GRW", "Stasis-OU", "URW-URW", "URW-GRW",
  "URW-OU","GRW-GRW", "GRW-OU", "OU-OU", "OU-GRW", "OU-URW", "OU-Stasis",
  "GRW-URW","GRW-Stasis", "URW-Stasis"
)

aicc_shift_extraction <- lapply(aicc_shift_extraction, function(x) {
  names(x) <- model_names
  return(x)
})

# create a dataframe with the results
model_names <- lapply(aicc_shift_extraction, function(x) {
  names(x)
})

model_aiccs <- lapply(aicc_shift_extraction, function(x) {
  unname(x)
})

aicc_shift <- Map(function(x, y) {
  data.frame(AICc = unlist(y), row.names = unlist(x))
}, model_names, model_aiccs)

  
###################################
## Get percent of the best AICcs ##
###################################

# combine the AICcs of no shift and shift models
aicc <- list()

for (i in 1:length(aicc_noshift)) {
  aicc[[i]] <- rbind(aicc_noshift[[i]], aicc_shift[[i]])
}

#------------------------------------------
# Find the best AICs for each timeseries
#------------------------------------------

# check which AICc value is the lowest
aicc_min <- lapply(aicc, function(x) {
  which.min(as.numeric(unlist(x)))
})

# get percentage
aicc_unlist <- unlist(aicc_min)
counts <- table(aicc_unlist)

# Create the outcomes for models with 0 iterations (model OU with moving optimum in this case)
modified_counts <- numeric(25)

for (i in 1:25) {
  if (i %in% names(counts)) {
    modified_counts[i] <- counts[[as.character(i)]]
  } else {
    modified_counts[i] <- 0
  }
}

percent <- (modified_counts/sum(modified_counts))*100
names(percent) <- c("GRW", "URW", "Stasis", "Strict stasis", "Decel", "Accel", "OU",
                    "OU mov. optm. (ancestral state)", "OU mov. optm.","Stasis-Stasis", 
                    "Stasis-URW", "Stasis-GRW", "Stasis-OU", "URW-URW", "URW-GRW", "URW-OU",
                    "GRW-GRW", "GRW-OU", "OU-OU", "OU-GRW", "OU-URW", "OU-Stasis", "GRW-URW",
                    "GRW-Stasis", "URW-Stasis")

# write to file
sink(file = "./results/percent_AICc_with_shift.txt")
percent
sink()

  
###################
## Test adequacy ##
###################

# Remove problematic timeseries from ln_datalong 
keep_TS3 <- !names(ln_datalong) %in% pblm_TS
data_aicc = ln_datalong[keep_TS3]

# Add a column with lowest AIC for each time series
for (i in 1:length(data_aicc)) {
  data_aicc[[i]]$Lowest_AICc <- aicc_min[i]
}


#----------------------------------------------------
# Filter the time series according to the best model
#----------------------------------------------------

categories <- c("GRW", "URW", "Stasis", "Strict stasis", "Decel", "Accel", "OU",
                "OU mov. optm. (ancestral state)", "OU mov. optm.", "Stasis_Stasis", 
                "Stasis_URW", "Stasis_GRW", "Stasis_OU", "URW_URW", "URW_GRW", "URW_OU",
                "GRW_GRW", "GRW_OU", "OU_OU", "OU_GRW", "OU_URW", "OU_Stasis", "GRW_URW",
                "GRW_Stasis", "URW_Stasis")

# Store the timeseries in different lists according to the best model fitted
result_list <- list()
for (i in 1:length(categories)) {
  category <- categories[i]
  # Filter data for the current category
  filtered_data <- Filter(function(x) x[[10]] == i, data_aicc)
  filtered_data <- lapply(filtered_data, function(x) { x[[10]] <- NULL; x })
  filtered_data <- lapply(filtered_data, function(x) {
    as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt)
  })
  
  # Store the result in the result_list and in the different best-model lists
  result_list[[category]] <- filtered_data
  assign(paste(category, sep = ""), filtered_data)
}

# remove problem time series
OU <- OU[names(OU) != 427]
OU <- OU[names(OU) != 428]

                          
#------------------------------------
# Splitting the shift models
#------------------------------------
                          
#Load previous version of paleoTS containing function to split PaleoTS objects
install.packages("remotes")
library(remotes)
install_version("paleoTS", version = "0.5-1")
library(paleoTS)

                          
### Stasis-Stasis ###
# Add a column with the best shift point to each timeseries
tsID_Stasis_Stasis = names(Stasis_Stasis)
model_results_Stasis_Stasis <- model_shift_results_clean[tsID_Stasis_Stasis]
for (i in tsID_Stasis_Stasis) {
  Stasis_Stasis[[i]]$Shift_point <- model_shift_results_clean[[i]][[1]]$parameters[["shift1"]]  #1 is for the model Stasis_Stasis, need to be change for the other models
}

#Splitting the model
Stasis_Stasis_subset1 = list()
Stasis_Stasis_subset2 = list()
for (i in 1:length(Stasis_Stasis)) {
  gg <- rep(1:2, c(Stasis_Stasis[[i]]$Shift_point, length(Stasis_Stasis[[i]]$mm) - Stasis_Stasis[[i]]$Shift_point))
  Stasis_Stasis_split = split4punc(Stasis_Stasis[[i]],gg, overlap=TRUE)
  Stasis_Stasis_subset1[[i]] = Stasis_Stasis_split[[1]]
  Stasis_Stasis_subset2[[i]] = Stasis_Stasis_split[[2]]
}


### Stasis-URW ###
# Add a column with the best shift point to each timeseries
tsID_Stasis_URW = names(Stasis_URW)
model_results_Stasis_URW <- model_shift_results_clean[tsID_Stasis_URW]
for (i in tsID_Stasis_URW) {
  Stasis_URW[[i]]$Shift_point <- model_shift_results_clean[[i]][[2]]$parameters[["shift1"]]
}
                          
#Splitting the model
Stasis_URW_subset1 = list()
Stasis_URW_subset2 = list()
for (i in 1:length(Stasis_URW)) {
  gg <- rep(1:2, c(Stasis_URW[[i]]$Shift_point, length(Stasis_URW[[i]]$mm) - Stasis_URW[[i]]$Shift_point))
  Stasis_URW_split = split4punc(Stasis_URW[[i]],gg, overlap=TRUE)
  Stasis_URW_subset1[[i]] = Stasis_URW_split[[1]]
  Stasis_URW_subset2[[i]] = Stasis_URW_split[[2]]
}


### Stasis-GRW ###
# Add a column with the best shift point to each timeseries
tsID_Stasis_GRW = names(Stasis_GRW)
model_results_Stasis_GRW <- model_shift_results_clean[tsID_Stasis_GRW]
for (i in tsID_Stasis_GRW) {
  Stasis_GRW[[i]]$Shift_point <- model_shift_results_clean[[i]][[3]]$parameters[["shift1"]]
}

#Splitting the model
Stasis_GRW_subset1 = list()
Stasis_GRW_subset2 = list()
for (i in 1:length(Stasis_GRW)) {
  gg <- rep(1:2, c(Stasis_GRW[[i]]$Shift_point, length(Stasis_GRW[[i]]$mm) - Stasis_GRW[[i]]$Shift_point))
  Stasis_GRW_split = split4punc(Stasis_GRW[[i]],gg, overlap=TRUE)
  Stasis_GRW_subset1[[i]] = Stasis_GRW_split[[1]]
  Stasis_GRW_subset2[[i]] = Stasis_GRW_split[[2]]
}


### URW-OU ###
# Add a column with the best shift point to each timeseries
tsID_Stasis_OU = names(Stasis_OU)
model_results_Stasis_OU <- model_shift_results_clean[tsID_Stasis_OU]
for (i in tsID_Stasis_OU) {
  Stasis_OU[[i]]$Shift_point <- model_shift_results_clean[[i]][[4]]$parameters[["shift1"]]
}

#Splitting the model
Stasis_OU_subset1 = list()
Stasis_OU_subset2 = list()
for (i in 1:length(Stasis_OU)) {
  gg <- rep(1:2, c(Stasis_OU[[i]]$Shift_point, length(Stasis_OU[[i]]$mm) - Stasis_OU[[i]]$Shift_point))
  Stasis_OU_split = split4punc(Stasis_OU[[i]],gg, overlap=TRUE)
  Stasis_OU_subset1[[i]] = Stasis_OU_split[[1]]
  Stasis_OU_subset2[[i]] = Stasis_OU_split[[2]]
}


### URW-URW ###
# Add a column with the best shift point to each timeseries
tsID_URW_URW = names(URW_URW)
model_results_URW_URW <- model_shift_results_clean[tsID_URW_URW]
for (i in tsID_URW_URW) {
  URW_URW[[i]]$Shift_point <- model_shift_results_clean[[i]][[5]]$parameters[["shift1"]]
}

#Splitting the model
URW_URW_subset1 = list()
URW_URW_subset2 = list()
for (i in 1:length(URW_URW)) {
  gg <- rep(1:2, c(URW_URW[[i]]$Shift_point, length(URW_URW[[i]]$mm) - URW_URW[[i]]$Shift_point))
  URW_URW_split = split4punc(URW_URW[[i]],gg, overlap=TRUE)
  URW_URW_subset1[[i]] = URW_URW_split[[1]]
  URW_URW_subset2[[i]] = URW_URW_split[[2]]
}


### URW-GRW ###
# Add a column with the best shift point to each timeseries
tsID_URW_GRW = names(URW_GRW)
model_results_URW_GRW <- model_shift_results_clean[tsID_URW_GRW]
for (i in tsID_URW_GRW) {
  URW_GRW[[i]]$Shift_point <- model_shift_results_clean[[i]][[6]]$parameters[["shift1"]]
}

#Splitting the model
URW_GRW_subset1 = list()
URW_GRW_subset2 = list()
for (i in 1:length(URW_GRW)) {
  gg <- rep(1:2, c(URW_GRW[[i]]$Shift_point, length(URW_GRW[[i]]$mm) - URW_GRW[[i]]$Shift_point))
  URW_GRW_split = split4punc(URW_GRW[[i]],gg, overlap=TRUE)
  URW_GRW_subset1[[i]] = URW_GRW_split[[1]]
  URW_GRW_subset2[[i]] = URW_GRW_split[[2]]
}


### URW-OU ###
# Add a column with the best shift point to each timeseries
tsID_URW_OU = names(URW_OU)
model_results_URW_OU <- model_shift_results_clean[tsID_URW_OU]
for (i in tsID_URW_OU) {
  URW_OU[[i]]$Shift_point <- model_shift_results_clean[[i]][[7]]$parameters[["shift1"]]
}

#Splitting the model
URW_OU_subset1 = list()
URW_OU_subset2 = list()
for (i in 1:length(URW_OU)) {
  gg <- rep(1:2, c(URW_OU[[i]]$Shift_point, length(URW_OU[[i]]$mm) - URW_OU[[i]]$Shift_point))
  URW_OU_split = split4punc(URW_OU[[i]],gg, overlap=TRUE)
  URW_OU_subset1[[i]] = URW_OU_split[[1]]
  URW_OU_subset2[[i]] = URW_OU_split[[2]]
}


### GRW-GRW ###
# Add a column with the best shift point to each timeseries
tsID_GRW_GRW = names(GRW_GRW)
model_results_GRW_GRW <- model_shift_results_clean[tsID_GRW_GRW]
for (i in tsID_GRW_GRW) {
  GRW_GRW[[i]]$Shift_point <- model_shift_results_clean[[i]][[8]]$parameters[["shift1"]]
}

#Splitting the model
GRW_GRW_subset1 = list()
GRW_GRW_subset2 = list()
for (i in 1:length(GRW_GRW)) {
  gg <- rep(1:2, c(GRW_GRW[[i]]$Shift_point, length(GRW_GRW[[i]]$mm) - GRW_GRW[[i]]$Shift_point))
  GRW_GRW_split = split4punc(GRW_GRW[[i]],gg, overlap=TRUE)
  GRW_GRW_subset1[[i]] = GRW_GRW_split[[1]]
  GRW_GRW_subset2[[i]] = GRW_GRW_split[[2]]
}


### GRW-OU ###
# Add a column with the best shift point to each timeseries
tsID_GRW_OU = names(GRW_OU)
model_results_GRW_OU <- model_shift_results_clean[tsID_GRW_OU]
for (i in tsID_GRW_OU) {
  GRW_OU[[i]]$Shift_point <- model_shift_results_clean[[i]][[9]]$parameters[["shift1"]]
}

#Splitting the model
GRW_OU_subset1 = list()
GRW_OU_subset2 = list()
for (i in 1:length(GRW_OU)) {
  gg <- rep(1:2, c(GRW_OU[[i]]$Shift_point, length(GRW_OU[[i]]$mm) - GRW_OU[[i]]$Shift_point))
  GRW_OU_split = split4punc(GRW_OU[[i]],gg, overlap=TRUE)
  GRW_OU_subset1[[i]] = GRW_OU_split[[1]]
  GRW_OU_subset2[[i]] = GRW_OU_split[[2]]
}


### OU-OU ###
# Add a column with the best shift point to each timeseries
tsID_OU_OU = names(OU_OU)
model_results_OU_OU <- model_shift_results_clean[tsID_OU_OU]
for (i in tsID_OU_OU) {
  OU_OU[[i]]$Shift_point <- model_shift_results_clean[[i]][[10]]$parameters[["shift1"]]
}

#Splitting the model
OU_OU_subset1 = list()
OU_OU_subset2 = list()
for (i in 1:length(OU_OU)) {
  gg <- rep(1:2, c(OU_OU[[i]]$Shift_point, length(OU_OU[[i]]$mm) - OU_OU[[i]]$Shift_point))
  OU_OU_split = split4punc(OU_OU[[i]],gg, overlap=TRUE)
  OU_OU_subset1[[i]] = OU_OU_split[[1]]
  OU_OU_subset2[[i]] = OU_OU_split[[2]]
}


### OU-GRW ###
# Add a column with the best shift point to each timeseries
tsID_OU_GRW = names(OU_GRW)
model_results_OU_GRW <- model_shift_results_clean[tsID_OU_GRW]
for (i in tsID_OU_GRW) {
  OU_GRW[[i]]$Shift_point <- model_shift_results_clean[[i]][[11]]$parameters[["shift1"]]
}

#Splitting the model
OU_GRW_subset1 = list()
OU_GRW_subset2 = list()
for (i in 1:length(OU_GRW)) {
  gg <- rep(1:2, c(OU_GRW[[i]]$Shift_point, length(OU_GRW[[i]]$mm) - OU_GRW[[i]]$Shift_point))
  OU_GRW_split = split4punc(OU_GRW[[i]],gg, overlap=TRUE)
  OU_GRW_subset1[[i]] = OU_GRW_split[[1]]
  OU_GRW_subset2[[i]] = OU_GRW_split[[2]]
}


### OU-URW ###
# Add a column with the best shift point to each timeseries
tsID_OU_URW = names(OU_URW)
model_results_OU_URW <- model_shift_results_clean[tsID_OU_URW]
for (i in tsID_OU_URW) {
  OU_URW[[i]]$Shift_point <- model_shift_results_clean[[i]][[12]]$parameters[["shift1"]]
}

#Splitting the model
OU_URW_subset1 = list()
OU_URW_subset2 = list()
for (i in 1:length(OU_URW)) {
  gg <- rep(1:2, c(OU_URW[[i]]$Shift_point, length(OU_URW[[i]]$mm) - OU_URW[[i]]$Shift_point))
  OU_URW_split = split4punc(OU_URW[[i]],gg, overlap=TRUE)
  OU_URW_subset1[[i]] = OU_URW_split[[1]]
  OU_URW_subset2[[i]] = OU_URW_split[[2]]
}


### OU-Stasis ###
# Add a column with the best shift point to each timeseries
tsID_OU_Stasis = names(OU_Stasis)
model_results_OU_Stasis <- model_shift_results_clean[tsID_OU_Stasis]
for (i in tsID_OU_Stasis) {
  OU_Stasis[[i]]$Shift_point <- model_shift_results_clean[[i]][[13]]$parameters[["shift1"]]
}

#Splitting the model
OU_Stasis_subset1 = list()
OU_Stasis_subset2 = list()
for (i in 1:length(OU_Stasis)) {
  gg <- rep(1:2, c(OU_Stasis[[i]]$Shift_point, length(OU_Stasis[[i]]$mm) - OU_Stasis[[i]]$Shift_point))
  OU_Stasis_split = split4punc(OU_Stasis[[i]],gg, overlap=TRUE)
  OU_Stasis_subset1[[i]] = OU_Stasis_split[[1]]
  OU_Stasis_subset2[[i]] = OU_Stasis_split[[2]]
}


### GRW-URW ###
# Add a column with the best shift point to each timeseries
tsID_GRW_URW = names(GRW_URW)
model_results_GRW_URW <- model_shift_results_clean[tsID_GRW_URW]
for (i in tsID_GRW_URW) {
  GRW_URW[[i]]$Shift_point <- model_shift_results_clean[[i]][[14]]$parameters[["shift1"]]
}

#Splitting the model
GRW_URW_subset1 = list()
GRW_URW_subset2 = list()
for (i in 1:length(GRW_URW)) {
  gg <- rep(1:2, c(GRW_URW[[i]]$Shift_point, length(GRW_URW[[i]]$mm) - GRW_URW[[i]]$Shift_point))
  GRW_URW_split = split4punc(GRW_URW[[i]],gg, overlap=TRUE)
  GRW_URW_subset1[[i]] = GRW_URW_split[[1]]
  GRW_URW_subset2[[i]] = GRW_URW_split[[2]]
}


### GRW-Stasis ###
# Add a column with the best shift point to each timeseries
tsID_GRW_Stasis = names(GRW_Stasis)
model_results_GRW_Stasis <- model_shift_results_clean[tsID_GRW_Stasis]
for (i in tsID_GRW_Stasis) {
  GRW_Stasis[[i]]$Shift_point <- model_shift_results_clean[[i]][[15]]$parameters[["shift1"]]
}

#Splitting the model
GRW_Stasis_subset1 = list()
GRW_Stasis_subset2 = list()
for (i in 1:length(GRW_Stasis)) {
  gg <- rep(1:2, c(GRW_Stasis[[i]]$Shift_point, length(GRW_Stasis[[i]]$mm) - GRW_Stasis[[i]]$Shift_point))
  GRW_Stasis_split = split4punc(GRW_Stasis[[i]],gg, overlap=TRUE)
  GRW_Stasis_subset1[[i]] = GRW_Stasis_split[[1]]
  GRW_Stasis_subset2[[i]] = GRW_Stasis_split[[2]]
}


### URW-Stasis ###
# Add a column with the best shift point to each timeseries
tsID_URW_Stasis = names(URW_Stasis)
model_results_URW_Stasis <- model_shift_results_clean[tsID_URW_Stasis]
for (i in tsID_URW_Stasis) {
  URW_Stasis[[i]]$Shift_point <- model_shift_results_clean[[i]][[16]]$parameters[["shift1"]]
}

#Splitting the model
URW_Stasis_subset1 = list()
URW_Stasis_subset2 = list()
for (i in 1:length(URW_Stasis)) {
  gg <- rep(1:2, c(URW_Stasis[[i]]$Shift_point, length(URW_Stasis[[i]]$mm) - URW_Stasis[[i]]$Shift_point))
  URW_Stasis_split = split4punc(URW_Stasis[[i]],gg, overlap=TRUE)
  URW_Stasis_subset1[[i]] = URW_Stasis_split[[1]]
  URW_Stasis_subset2[[i]] = URW_Stasis_split[[2]]
}


#------------------------------------
# Testing the adequacy of the models
#------------------------------------

# test adequacy
GRW_adeq <- mclapply(GRW, fit3adequacy.trend, plot = FALSE)
URW_adeq <- mclapply(URW, fit3adequacy.RW, plot = FALSE)
stasis_adeq <- mclapply(stasis, fit4adequacy.stasis, plot = FALSE) 
strict_stasis_adeq <- mclapply(strict_stasis, fit4adequacy.stasis, plot = FALSE)
decel_adeq <- mclapply(decel, fit3adequacy.decel, plot = FALSE)
accel_adeq <- mclapply(accel, fit3adequacy.RW, plot = FALSE)
OU_adeq <- mclapply(OU, fit3adequacy.OU, plot = FALSE)
OU_mov_opt_anc_adeq <- mclapply(OU_mov_opt_anc, fit3adequacy.OU, plot = FALSE)
OU_mov_opt_adeq <- mclapply(OU_mov_opt, fit3adequacy.OU, plot = FALSE)

Stasis_Stasis_subset1_adeq <- mclapply(Stasis_Stasis_subset1, fit4adequacy.stasis, plot = FALSE)
Stasis_Stasis_subset2_adeq <- mclapply(Stasis_Stasis_subset2, fit4adequacy.stasis, plot = FALSE)
names(Stasis_Stasis_subset1_adeq) = tsID_Stasis_Stasis
names(Stasis_Stasis_subset2_adeq) = tsID_Stasis_Stasis

Stasis_URW_subset1_adeq <- mclapply(Stasis_URW_subset1, fit4adequacy.stasis, plot = FALSE)
Stasis_URW_subset2_adeq <- mclapply(Stasis_URW_subset2, fit3adequacy.RW, plot = FALSE)
names(Stasis_URW_subset1_adeq) = tsID_Stasis_URW
names(Stasis_URW_subset2_adeq) = tsID_Stasis_URW

Stasis_GRW_subset1_adeq <- mclapply(Stasis_GRW_subset1, fit4adequacy.stasis, plot = FALSE)
Stasis_GRW_subset2_adeq <- mclapply(Stasis_GRW_subset2, fit3adequacy.trend, plot = FALSE)
names(Stasis_GRW_subset1_adeq) = tsID_Stasis_GRW
names(Stasis_GRW_subset2_adeq) = tsID_Stasis_GRW

Stasis_OU_subset1_adeq <- mclapply(Stasis_OU_subset1, fit4adequacy.stasis, plot = FALSE)
Stasis_OU_subset2_adeq <- mclapply(Stasis_OU_subset2, fit3adequacy.OU, plot = FALSE)
names(Stasis_OU_subset1_adeq) = tsID_Stasis_OU
names(Stasis_OU_subset2_adeq) = tsID_Stasis_OU

URW_URW_subset1_adeq <- mclapply(URW_URW_subset1, fit3adequacy.RW, plot = FALSE)
URW_URW_subset2_adeq <- mclapply(URW_URW_subset2, fit3adequacy.RW, plot = FALSE)
names(URW_URW_subset1_adeq) = tsID_URW_URW
names(URW_URW_subset2_adeq) = tsID_URW_URW

URW_GRW_subset1_adeq <- mclapply(URW_GRW_subset1, fit3adequacy.RW, plot = FALSE)
URW_GRW_subset2_adeq <- mclapply(URW_GRW_subset2, fit3adequacy.trend, plot = FALSE)
names(URW_GRW_subset1_adeq) = tsID_URW_GRW
names(URW_GRW_subset2_adeq) = tsID_URW_GRW

URW_OU_subset1_adeq <- mclapply(URW_OU_subset1, fit3adequacy.RW, plot = FALSE)
URW_OU_subset2_adeq <- mclapply(URW_OU_subset2, fit3adequacy.OU, plot = FALSE)
names(URW_OU_subset1_adeq) = tsID_URW_OU
names(URW_OU_subset2_adeq) = tsID_URW_OU

GRW_GRW_subset1_adeq <- mclapply(GRW_GRW_subset1, fit3adequacy.trend, plot = FALSE)
GRW_GRW_subset2_adeq <- mclapply(GRW_GRW_subset2, fit3adequacy.trend, plot = FALSE)
names(GRW_GRW_subset1_adeq) = tsID_GRW_GRW
names(GRW_GRW_subset2_adeq) = tsID_GRW_GRW

GRW_OU_subset1_adeq <- mclapply(GRW_OU_subset1, fit3adequacy.trend, plot = FALSE)
GRW_OU_subset2_adeq <- mclapply(GRW_OU_subset2, fit3adequacy.OU, plot = FALSE)
names(GRW_OU_subset1_adeq) = tsID_GRW_OU
names(GRW_OU_subset2_adeq) = tsID_GRW_OU

OU_OU_subset1_adeq <- mclapply(OU_OU_subset1, fit3adequacy.OU, plot = FALSE)
OU_OU_subset2_adeq <- mclapply(OU_OU_subset2, fit3adequacy.OU, plot = FALSE)
names(OU_OU_subset1_adeq) = tsID_OU_OU
names(OU_OU_subset2_adeq) = tsID_OU_OU

OU_GRW_subset1_adeq <- mclapply(OU_GRW_subset1, fit3adequacy.OU, plot = FALSE)
OU_GRW_subset2_adeq <- mclapply(OU_GRW_subset2, fit3adequacy.trend, plot = FALSE)
names(OU_GRW_subset1_adeq) = tsID_OU_GRW
names(OU_GRW_subset2_adeq) = tsID_OU_GRW

OU_URW_subset1_adeq <- mclapply(OU_URW_subset1, fit3adequacy.OU, plot = FALSE)
OU_URW_subset2_adeq <- mclapply(OU_URW_subset2, fit3adequacy.RW, plot = FALSE)
names(OU_URW_subset1_adeq) = tsID_OU_URW
names(OU_URW_subset2_adeq) = tsID_OU_URW

OU_Stasis_subset1_adeq <- mclapply(OU_Stasis_subset1, fit3adequacy.OU, plot = FALSE)
OU_Stasis_subset2_adeq <- mclapply(OU_Stasis_subset2, fit4adequacy.stasis, plot = FALSE)
names(OU_Stasis_subset1_adeq) = tsID_OU_Stasis
names(OU_Stasis_subset2_adeq) = tsID_OU_Stasis

GRW_URW_subset1_adeq <- mclapply(GRW_URW_subset1, fit3adequacy.trend, plot = FALSE)
GRW_URW_subset2_adeq <- mclapply(GRW_URW_subset2, fit3adequacy.RW, plot = FALSE)
names(GRW_URW_subset1_adeq) = tsID_GRW_URW
names(GRW_URW_subset2_adeq) = tsID_GRW_URW

GRW_Stasis_subset1_adeq <- mclapply(GRW_Stasis_subset1, fit3adequacy.trend, plot = FALSE)
GRW_Stasis_subset2_adeq <- mclapply(GRW_Stasis_subset2, fit4adequacy.stasis, plot = FALSE)
names(GRW_Stasis_subset1_adeq) = tsID_GRW_Stasis
names(GRW_Stasis_subset2_adeq) = tsID_GRW_Stasis

URW_Stasis_subset1_adeq <- mclapply(URW_Stasis_subset1, fit3adequacy.RW, plot = FALSE)
URW_Stasis_subset2_adeq <- mclapply(URW_Stasis_subset2, fit4adequacy.stasis, plot = FALSE)
names(URW_Stasis_subset1_adeq) = tsID_URW_Stasis
names(URW_Stasis_subset2_adeq) = tsID_URW_Stasis

                          
# get only adequate time series
GRW_adeq_passed <- adequate3tests(GRW_adeq)
URW_adeq_passed <- adequate3tests(URW_adeq)
stasis_adeq_passed <- adequate4tests(stasis_adeq)
strict_stasis_adeq_passed <- adequate4tests(strict_stasis_adeq)
decel_adeq_passed <- adequate3tests(decel_adeq)
accel_adeq_passed <- adequate3tests(accel_adeq)
OU_adeq_passed <- adequate3tests(OU_adeq)
OU_mov_opt_anc_adeq_passed <- adequate3tests(OU_mov_opt_anc_adeq)
OU_mov_opt_adeq_passed <- adequate3tests(OU_mov_opt_adeq)

Stasis_Stasis_subset1_adeq_passed <- adequate4tests(Stasis_Stasis_subset1_adeq)
Stasis_Stasis_subset2_adeq_passed <- adequate4tests(Stasis_Stasis_subset2_adeq)

Stasis_URW_subset1_adeq_passed <- adequate4tests(Stasis_URW_subset1_adeq)
Stasis_URW_subset2_adeq_passed <- adequate3tests(Stasis_URW_subset2_adeq)

Stasis_GRW_subset1_adeq_passed <- adequate4tests(Stasis_GRW_subset1_adeq)
Stasis_GRW_subset2_adeq_passed <- adequate3tests(Stasis_GRW_subset2_adeq)

Stasis_OU_subset1_adeq_passed <- adequate4tests(Stasis_OU_subset1_adeq)
Stasis_OU_subset2_adeq_passed <- adequate3tests(Stasis_OU_subset2_adeq)

URW_URW_subset1_adeq_passed <- adequate3tests(URW_URW_subset1_adeq)
URW_URW_subset2_adeq_passed <- adequate3tests(URW_URW_subset2_adeq)

URW_GRW_subset1_adeq_passed <- adequate3tests(URW_GRW_subset1_adeq)
URW_GRW_subset2_adeq_passed <- adequate3tests(URW_GRW_subset2_adeq)

URW_OU_subset1_adeq_passed <- adequate3tests(URW_OU_subset1_adeq)
URW_OU_subset2_adeq_passed <- adequate3tests(URW_OU_subset2_adeq)

GRW_GRW_subset1_adeq_passed <- adequate3tests(GRW_GRW_subset1_adeq)
GRW_GRW_subset2_adeq_passed <- adequate3tests(GRW_GRW_subset2_adeq)

GRW_OU_subset1_adeq_passed <- adequate3tests(GRW_OU_subset1_adeq)
GRW_OU_subset2_adeq_passed <- adequate3tests(GRW_OU_subset2_adeq)

OU_OU_subset1_adeq_passed <- adequate3tests(OU_OU_subset1_adeq)
OU_OU_subset2_adeq_passed <- adequate3tests(OU_OU_subset2_adeq)

OU_GRW_subset1_adeq_passed <- adequate3tests(OU_GRW_subset1_adeq)
OU_GRW_subset2_adeq_passed <- adequate3tests(OU_GRW_subset2_adeq)

OU_URW_subset1_adeq_passed <- adequate3tests(OU_URW_subset1_adeq)
OU_URW_subset2_adeq_passed <- adequate3tests(OU_URW_subset2_adeq)

OU_Stasis_subset1_adeq_passed <- adequate3tests(OU_Stasis_subset1_adeq)
OU_Stasis_subset2_adeq_passed <- adequate4tests(OU_Stasis_subset2_adeq)

GRW_URW_subset1_adeq_passed <- adequate3tests(GRW_URW_subset1_adeq)
GRW_URW_subset2_adeq_passed <- adequate3tests(GRW_URW_subset2_adeq)

GRW_Stasis_subset1_adeq_passed <- adequate3tests(GRW_Stasis_subset1_adeq)
GRW_Stasis_subset2_adeq_passed <- adequate4tests(GRW_Stasis_subset2_adeq)

URW_Stasis_subset1_adeq_passed <- adequate3tests(URW_Stasis_subset1_adeq)
URW_Stasis_subset2_adeq_passed <- adequate4tests(URW_Stasis_subset2_adeq)

                          
# merge split timeseries if the two subsets passed the adequacy tests

Stasis_Stasis_adeq_passed <- Stasis_Stasis_subset1_adeq_passed[intersect(names(Stasis_Stasis_subset1_adeq_passed), names(Stasis_Stasis_subset2_adeq_passed))]
Stasis_URW_adeq_passed <- Stasis_URW_subset1_adeq_passed[intersect(names(Stasis_URW_subset1_adeq_passed), names(Stasis_URW_subset2_adeq_passed))]
Stasis_GRW_adeq_passed <- Stasis_GRW_subset1_adeq_passed[intersect(names(Stasis_GRW_subset1_adeq_passed), names(Stasis_GRW_subset2_adeq_passed))]
Stasis_OU_adeq_passed <- Stasis_OU_subset1_adeq_passed[intersect(names(Stasis_OU_subset1_adeq_passed), names(Stasis_OU_subset2_adeq_passed))]
URW_URW_adeq_passed <- URW_URW_subset1_adeq_passed[intersect(names(URW_URW_subset1_adeq_passed), names(URW_URW_subset2_adeq_passed))]
URW_GRW_adeq_passed <- URW_GRW_subset1_adeq_passed[intersect(names(URW_GRW_subset1_adeq_passed), names(URW_GRW_subset2_adeq_passed))]
URW_OU_adeq_passed <- URW_OU_subset1_adeq_passed[intersect(names(URW_OU_subset1_adeq_passed), names(URW_OU_subset2_adeq_passed))]
GRW_GRW_adeq_passed <- GRW_GRW_subset1_adeq_passed[intersect(names(GRW_GRW_subset1_adeq_passed), names(GRW_GRW_subset2_adeq_passed))]
GRW_OU_adeq_passed <- GRW_OU_subset1_adeq_passed[intersect(names(GRW_OU_subset1_adeq_passed), names(GRW_OU_subset2_adeq_passed))]
OU_OU_adeq_passed <- OU_OU_subset1_adeq_passed[intersect(names(OU_OU_subset1_adeq_passed), names(OU_OU_subset2_adeq_passed))]
OU_GRW_adeq_passed <- OU_GRW_subset1_adeq_passed[intersect(names(OU_GRW_subset1_adeq_passed), names(OU_GRW_subset2_adeq_passed))]
OU_URW_adeq_passed <- OU_URW_subset1_adeq_passed[intersect(names(OU_URW_subset1_adeq_passed), names(OU_URW_subset2_adeq_passed))]
OU_Stasis_adeq_passed <- OU_Stasis_subset1_adeq_passed[intersect(names(OU_Stasis_subset1_adeq_passed), names(OU_Stasis_subset2_adeq_passed))]
GRW_URW_adeq_passed <- GRW_URW_subset1_adeq_passed[intersect(names(GRW_URW_subset1_adeq_passed), names(GRW_URW_subset2_adeq_passed))]
GRW_Stasis_adeq_passed <- GRW_Stasis_subset1_adeq_passed[intersect(names(GRW_Stasis_subset1_adeq_passed), names(GRW_Stasis_subset2_adeq_passed))]
URW_Stasis_adeq_passed <- URW_Stasis_subset1_adeq_passed[intersect(names(URW_Stasis_subset1_adeq_passed), names(URW_Stasis_subset2_adeq_passed))]            

# get percentage passed


                          
