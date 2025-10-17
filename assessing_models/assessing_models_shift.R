##########################################################################
## Evolutionary rates and time series scaling with shift model included ##
##########################################################################

#R version 4.2.1
#evoTS version 1.0.3
#adePEM new models version
#paleoTS version 0.6.2

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


source("C:/Users/marionth/OneDrive - Universitetet i Oslo/Skrivebord/PhD/Github/assessing_models_evolution/assessing_models/assessing_models_uni_functions.R")


# set working directory
setwd("C:/Users/marionth/OneDrive - Universitetet i Oslo/Skrivebord/PhD/Github/assessing_models_evolution/assessing_models")

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
# EXCLUDING THE SHORTEST TIMESERIES (Keep timeseries only if containing 20 or more steps)
#---------------------------------------------------------------------------------------------

metadatalong <- matrix(nrow = 0, ncol = ncol(metadata))
timeserieslong <- matrix(nrow = 0, ncol = ncol(timeseries))
tsIDremoved <- c()

# Filter rows of the metadata file based on step number 
for (i in 1:nrow(metadata)) {                          
  if (metadata[i, "steps"] >= 20) {   
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
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$N, tt = x$tt, oldest = "first")   ###NEED TO ADD TIME RESCALING
})

                        
#####################################################
## Fit models including shift and find best (AICc) ##
#####################################################

# test all possible univariate models from evoTS on every timeseries
# for paleoTS v0.5.3
# model_noshift_results <- mclapply(ln_datalong, fit.all.univariate) 

# for paleoTS v0.6.1
model_noshift_results <- list()
for(i in 1:length(ln_datalong)){
  try(model_noshift_results[[i]] <- fit.all.univariate(ln_datalong[[i]]))
}

# add time series IDs
names_list <- names(ln_datalong)
names(model_noshift_results) <- names_list

                 
# test all possible shift models from evoTS on time series (To do on the HPC)
fit_mode_shift <- function(ln_datalong) {
  models_list <- c("Stasis", "URW", "GRW", "OU")
  store_results <- list()
  k <- 0
  for (i in 1:4) {
    model1 <- models_list[i]  
    for (j in 1:4) {
      model2 <- models_list[j]
        fit_result <- fit.mode.shift(ln_datalong, model1, model2, minb = 10)
        k = k + 1
        store_results[[k]] <- fit_result
    } 
  }
  return(store_results)
}

#model_shift_results <- mclapply(ln_datalong, fit_mode_shift)

#save(model_noshift_results, model_shift_results, file = "./results_paleoTS_v0.6.1/Results_fit_shiftmodels.RData")
load("./results_paleoTS_v0.6.1/Results_fit_shiftmodels.RData") # Loading results from the HPC


# remove time series that cannot be processed by the loglikelihood function 
ln_datalong2 = ln_datalong
ln_datalong2 = ln_datalong2[-which(sapply(model_shift_results, is.null))]
model_noshift_results = model_noshift_results[-which(sapply(model_noshift_results, is.null))]

# remove time series that can not be processed by the loglikelihood function from the results with a shift
model_shift_results_names <- names(model_shift_results)
model_noshift_results_names <- names(model_noshift_results)
removed_TS <- setdiff(model_shift_results_names, model_noshift_results_names)
model_shift_results[removed_TS] <- NULL

# remove time series that can not be processed by the loglikelihood function from ln_datalong
ln_datalong[removed_TS] <- NULL


########################
##  Extract the AICcs ##
########################

### Remove problematic timeseries ###
pblm_TS <- c("584","585","75","427","428","574","368") #issues during the adequacy tests
keep_TS <- !names(model_shift_results) %in% pblm_TS
model_shift_results <- model_shift_results[keep_TS]
keep_TS2 <- !names(model_noshift_results) %in% pblm_TS
model_noshift_results <- model_noshift_results[keep_TS2]

#------------------------------------------
# Extract AICcs for the no shift models
#------------------------------------------

aicc_noshift <- lapply(model_noshift_results, function(x) x[(names(x) %in% c("AICc"))]) #USE model_noshift_results_clean if those time series are still problematic

#------------------------------------------
# Extract AICcs for the shift models
#------------------------------------------

# extract AICc values of shift models on all results
aicc_shift_extraction <- lapply(model_shift_results, function(x) { ### remettre model_shift_result when problem will be solved
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

aicc <- list()
for (i in 1:length(aicc_noshift)) {
  aicc_name <- names(aicc_noshift)[i]  # Extract the name of the current sublist in aicc_noshift
  aicc[[aicc_name]] <- rbind(aicc_noshift[[i]], aicc_shift[[i]])
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
aicc_results <- table(aicc_unlist)

# Create the outcomes of models with 0 time series 
aicc_results_complete <- numeric(25)

for (i in 1:25) {
  if (i %in% names(aicc_results)) {
    aicc_results_complete[i] <- aicc_results[[as.character(i)]]
  } else {
    aicc_results_complete[i] <- 0
  }
}

names(aicc_results_complete) <- c("GRW", "URW", "Stasis", "Strict stasis", "Decel", "Accel", "OU",
                                  "OU mov. optm. (ancestral state)", "OU mov. optm.","Stasis-Stasis", 
                                  "Stasis-URW", "Stasis-GRW", "Stasis-OU", "URW-URW", "URW-GRW", "URW-OU",
                                  "GRW-GRW", "GRW-OU", "OU-OU", "OU-GRW", "OU-URW", "OU-Stasis", "GRW-URW",
                                  "GRW-Stasis", "URW-Stasis")

aicc_results_complete <- data.frame(model = names(aicc_results_complete), count = unname(aicc_results_complete))
aicc_results_complete$percentage <- (aicc_results_complete$count/sum(aicc_results_complete$count))*100

percent2 <- sum(aicc_results_complete$percent[4:25])
percent3 <- sum(aicc_results_complete$percent[5:25])
percent4 <- sum(aicc_results_complete$percent[10:25])

# write to file
sink(file = "./results_paleoTS_v0.6.1/AICc_results_with_shift.txt")
aicc_results_complete
paste("Total number of time series investigated:", length(aicc))
paste("Percentage of time series not described by URW, GRW or stasis:", percent2)
paste("Percentage of time series not described by URW, GRW, stasis or strict stasis:", percent3)
paste("Percentage of time series described by models with shift:", percent4)
sink()


#--------------------------------------------------------------------------------------------------
# Time series described by more than one model (AICcs with a difference inferior to 2 units)
#--------------------------------------------------------------------------------------------------

aicc_filtered = list()
aicc_filtered_names = list()
threshold = 2

for (i in 1:length(aicc)) {
  for (j in 1:nrow(aicc[[1]][])) {
    if (aicc[[i]][j,] != aicc[[i]][aicc_min[[i]],]) {
      if (any(abs(aicc[[i]][j,] - aicc[[i]][aicc_min[[i]],]) <= threshold)) {
        aicc_filtered = c(aicc_filtered, list(aicc[[i]]))
        aicc_filtered_names = c(aicc_filtered_names, names(aicc)[i])
        break  # Exit the inner loop once the threshold is met
      } 
    }
  }
}

names(aicc_filtered) <- aicc_filtered_names

sink(file = "./results_paleoTS_v0.6.1/AICc_filter_with_shift.txt")
paste("Total number of time series investigated:", length(aicc))
paste("Total number of time series described by more than one model:", length(aicc_filtered))
paste("Percentage of time series filtered:", length(aicc_filtered)/length(aicc)*100)
sink()

  
###################
## Test adequacy ##
###################


data_aicc = list()

### Remove problematic and too short time series from ln_datalong ###
keep_TS3 <- !names(ln_datalong) %in% pblm_TS
data_aicc <- ln_datalong[keep_TS3]


# Add a column with lowest AIC for each time series
for (i in 1:length(data_aicc)) {
  data_aicc[[i]]$Lowest_AICc <- aicc_min[[i]]
}

#----------------------------------------------------
# Filter the time series according to the best model
#----------------------------------------------------


categories <- c("GRW", "URW", "Stasis", "Strict_stasis", "Decel", "Accel", "OU",
                "OU_mov_opt_anc", "OU_mov_opt", "Stasis_Stasis", 
                "Stasis_URW", "Stasis_GRW", "Stasis_OU", "URW_URW", "URW_GRW", "URW_OU",
                "GRW_GRW", "GRW_OU", "OU_OU", "OU_GRW", "OU_URW", "OU_Stasis", "GRW_URW",
                "GRW_Stasis", "URW_Stasis")

# Create a list to store the results
result_list <- list()

for (i in 1:length(categories)) {
  category <- categories[i]
  
  # Filter data for the current category
  filtered_data <- Filter(function(x) x[[10]] == i, data_aicc)
  filtered_data <- lapply(filtered_data, function(x) { x[[10]] <- NULL; x })
  filtered_data <- lapply(filtered_data, function(x) {
    as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt)
  })
  
  # Store the result in the result_list
  result_list[[category]] <- filtered_data
  assign(paste(category, sep = ""), filtered_data)
}

                          
#-----------------------------------------------------------------------
# Splitting the timeseries best described by models with a shift models
#-----------------------------------------------------------------------
                          
### Stasis-Stasis ###
# Add a column with the best shift point to each timeseries
tsID_Stasis_Stasis = names(Stasis_Stasis)
model_results_Stasis_Stasis <- model_shift_results[tsID_Stasis_Stasis]

for (i in tsID_Stasis_Stasis) {
  Stasis_Stasis[[i]]$Shift_point <- model_shift_results[[i]][[1]]$parameters[["shift1"]]  #1 is for the model Stasis_Stasis, need to be changed for the other models
}

#Splitting the model
Stasis_Stasis_subset1 = list()
Stasis_Stasis_subset2 = list()

if (length(Stasis_Stasis) > 0) {
  for (i in 1:length(Stasis_Stasis)) {
    gg <- rep(1:2, c(Stasis_Stasis[[i]]$Shift_point, length(Stasis_Stasis[[i]]$mm) - Stasis_Stasis[[i]]$Shift_point))
    Stasis_Stasis_split = paleoTS:::split4punc(Stasis_Stasis[[i]],gg, overlap=TRUE)
    Stasis_Stasis_subset1[[i]] = Stasis_Stasis_split[[1]]
    Stasis_Stasis_subset2[[i]] = Stasis_Stasis_split[[2]]
  }
}


### Stasis-URW ###
# Add a column with the best shift point to each timeseries
tsID_Stasis_URW = names(Stasis_URW)
model_results_Stasis_URW <- model_shift_results[tsID_Stasis_URW]

for (i in tsID_Stasis_URW) {
  Stasis_URW[[i]]$Shift_point <- model_shift_results[[i]][[2]]$parameters[["shift1"]]
}

#Splitting the model
Stasis_URW_subset1 = list()
Stasis_URW_subset2 = list()

if (length(Stasis_URW) > 0) {
  for (i in 1:length(Stasis_URW)) {
    gg <- rep(1:2, c(Stasis_URW[[i]]$Shift_point, length(Stasis_URW[[i]]$mm) - Stasis_URW[[i]]$Shift_point))
   Stasis_URW_split = paleoTS:::split4punc(Stasis_URW[[i]],gg, overlap=TRUE)
   Stasis_URW_subset1[[i]] = Stasis_URW_split[[1]]
   Stasis_URW_subset2[[i]] = Stasis_URW_split[[2]]
  }
}


### Stasis-GRW ###
# Add a column with the best shift point to each timeseries
tsID_Stasis_GRW = names(Stasis_GRW)
model_results_Stasis_GRW <- model_shift_results[tsID_Stasis_GRW]

for (i in tsID_Stasis_GRW) {
  Stasis_GRW[[i]]$Shift_point <- model_shift_results[[i]][[3]]$parameters[["shift1"]]
}

#Splitting the model
Stasis_GRW_subset1 = list()
Stasis_GRW_subset2 = list()

if (length(Stasis_GRW) > 0) {
  for (i in 1:length(Stasis_GRW)) {
    gg <- rep(1:2, c(Stasis_GRW[[i]]$Shift_point, length(Stasis_GRW[[i]]$mm) - Stasis_GRW[[i]]$Shift_point))
    Stasis_GRW_split = paleoTS:::split4punc(Stasis_GRW[[i]],gg, overlap=TRUE)
    Stasis_GRW_subset1[[i]] = Stasis_GRW_split[[1]]
    Stasis_GRW_subset2[[i]] = Stasis_GRW_split[[2]]
  }
}


### URW-OU ###
# Add a column with the best shift point to each timeseries
tsID_Stasis_OU = names(Stasis_OU)
model_results_Stasis_OU <- model_shift_results[tsID_Stasis_OU]

for (i in tsID_Stasis_OU) {
  Stasis_OU[[i]]$Shift_point <- model_shift_results[[i]][[4]]$parameters[["shift1"]]
}

#Splitting the model
Stasis_OU_subset1 = list()
Stasis_OU_subset2 = list()

if (length(Stasis_OU) > 0) {
  for (i in 1:length(Stasis_OU)) {
    gg <- rep(1:2, c(Stasis_OU[[i]]$Shift_point, length(Stasis_OU[[i]]$mm) - Stasis_OU[[i]]$Shift_point))
    Stasis_OU_split = paleoTS:::split4punc(Stasis_OU[[i]],gg, overlap=TRUE)
    Stasis_OU_subset1[[i]] = Stasis_OU_split[[1]]
    Stasis_OU_subset2[[i]] = Stasis_OU_split[[2]]
  }
}


### URW-URW ###
# Add a column with the best shift point to each timeseries
tsID_URW_URW = names(URW_URW)
model_results_URW_URW <- model_shift_results[tsID_URW_URW]

for (i in tsID_URW_URW) {
  URW_URW[[i]]$Shift_point <- model_shift_results[[i]][[5]]$parameters[["shift1"]]
}

#Splitting the model
URW_URW_subset1 = list()
URW_URW_subset2 = list()

if (length(URW_URW) > 0) {
  for (i in 1:length(URW_URW)) {
    gg <- rep(1:2, c(URW_URW[[i]]$Shift_point, length(URW_URW[[i]]$mm) - URW_URW[[i]]$Shift_point))
    URW_URW_split = paleoTS:::split4punc(URW_URW[[i]],gg, overlap=TRUE)
    URW_URW_subset1[[i]] = URW_URW_split[[1]]
    URW_URW_subset2[[i]] = URW_URW_split[[2]]
  }
}


### URW-GRW ###
# Add a column with the best shift point to each timeseries
tsID_URW_GRW = names(URW_GRW)
model_results_URW_GRW <- model_shift_results[tsID_URW_GRW]

for (i in tsID_URW_GRW) {
  URW_GRW[[i]]$Shift_point <- model_shift_results[[i]][[6]]$parameters[["shift1"]]
}

#Splitting the model
URW_GRW_subset1 = list()
URW_GRW_subset2 = list()

if (length(URW_GRW) > 0) {
  for (i in 1:length(URW_GRW)) {
    gg <- rep(1:2, c(URW_GRW[[i]]$Shift_point, length(URW_GRW[[i]]$mm) - URW_GRW[[i]]$Shift_point))
    URW_GRW_split = paleoTS:::split4punc(URW_GRW[[i]],gg, overlap=TRUE)
    URW_GRW_subset1[[i]] = URW_GRW_split[[1]]
    URW_GRW_subset2[[i]] = URW_GRW_split[[2]]
  }
}


### URW-OU ###
# Add a column with the best shift point to each timeseries
tsID_URW_OU = names(URW_OU)
model_results_URW_OU <- model_shift_results[tsID_URW_OU]

for (i in tsID_URW_OU) {
  URW_OU[[i]]$Shift_point <- model_shift_results[[i]][[7]]$parameters[["shift1"]]
}

#Splitting the model
URW_OU_subset1 = list()
URW_OU_subset2 = list()

if (length(URW_OU) > 0) {
  for (i in 1:length(URW_OU)) {
   gg <- rep(1:2, c(URW_OU[[i]]$Shift_point, length(URW_OU[[i]]$mm) - URW_OU[[i]]$Shift_point))
   URW_OU_split = paleoTS:::split4punc(URW_OU[[i]],gg, overlap=TRUE)
   URW_OU_subset1[[i]] = URW_OU_split[[1]]
   URW_OU_subset2[[i]] = URW_OU_split[[2]]
  }
}


### GRW-GRW ###
# Add a column with the best shift point to each timeseries
tsID_GRW_GRW = names(GRW_GRW)
model_results_GRW_GRW <- model_shift_results[tsID_GRW_GRW]

for (i in tsID_GRW_GRW) {
  GRW_GRW[[i]]$Shift_point <- model_shift_results[[i]][[8]]$parameters[["shift1"]]
}

#Splitting the model
GRW_GRW_subset1 = list()
GRW_GRW_subset2 = list()

if (length(GRW_GRW) > 0) {
  for (i in 1:length(GRW_GRW)) {
   gg <- rep(1:2, c(GRW_GRW[[i]]$Shift_point, length(GRW_GRW[[i]]$mm) - GRW_GRW[[i]]$Shift_point))
   GRW_GRW_split = paleoTS:::split4punc(GRW_GRW[[i]],gg, overlap=TRUE)
   GRW_GRW_subset1[[i]] = GRW_GRW_split[[1]]
   GRW_GRW_subset2[[i]] = GRW_GRW_split[[2]]
  }
}


### GRW-OU ###
# Add a column with the best shift point to each timeseries
tsID_GRW_OU = names(GRW_OU)
model_results_GRW_OU <- model_shift_results[tsID_GRW_OU]

for (i in tsID_GRW_OU) {
  GRW_OU[[i]]$Shift_point <- model_shift_results[[i]][[9]]$parameters[["shift1"]]
}

#Splitting the model
GRW_OU_subset1 = list()
GRW_OU_subset2 = list()

if (length(GRW_OU) > 0) {
  for (i in 1:length(GRW_OU)) {
    gg <- rep(1:2, c(GRW_OU[[i]]$Shift_point, length(GRW_OU[[i]]$mm) - GRW_OU[[i]]$Shift_point))
    GRW_OU_split = paleoTS:::split4punc(GRW_OU[[i]],gg, overlap=TRUE)
    GRW_OU_subset1[[i]] = GRW_OU_split[[1]]
    GRW_OU_subset2[[i]] = GRW_OU_split[[2]]
  }
}


### OU-OU ###
# Add a column with the best shift point to each timeseries
tsID_OU_OU = names(OU_OU)
model_results_OU_OU <- model_shift_results[tsID_OU_OU]

for (i in tsID_OU_OU) {
  OU_OU[[i]]$Shift_point <- model_shift_results[[i]][[10]]$parameters[["shift1"]]
}

#Splitting the model
OU_OU_subset1 = list()
OU_OU_subset2 = list()

if (length(OU_OU) > 0) {
  for (i in 1:length(OU_OU)) {
   gg <- rep(1:2, c(OU_OU[[i]]$Shift_point, length(OU_OU[[i]]$mm) - OU_OU[[i]]$Shift_point))
   OU_OU_split = paleoTS:::split4punc(OU_OU[[i]],gg, overlap=TRUE)
   OU_OU_subset1[[i]] = OU_OU_split[[1]]
   OU_OU_subset2[[i]] = OU_OU_split[[2]]
  }
}

### OU-GRW ###
# Add a column with the best shift point to each timeseries
tsID_OU_GRW = names(OU_GRW)
model_results_OU_GRW <- model_shift_results[tsID_OU_GRW]

for (i in tsID_OU_GRW) {
  OU_GRW[[i]]$Shift_point <- model_shift_results[[i]][[11]]$parameters[["shift1"]]
}

#Splitting the model
OU_GRW_subset1 = list()
OU_GRW_subset2 = list()

if (length(OU_GRW) > 0) {
  for (i in 1:length(OU_GRW)) {
   gg <- rep(1:2, c(OU_GRW[[i]]$Shift_point, length(OU_GRW[[i]]$mm) - OU_GRW[[i]]$Shift_point))
   OU_GRW_split = paleoTS:::split4punc(OU_GRW[[i]],gg, overlap=TRUE)
   OU_GRW_subset1[[i]] = OU_GRW_split[[1]]
   OU_GRW_subset2[[i]] = OU_GRW_split[[2]]
  }
}


### OU-URW ###
# Add a column with the best shift point to each timeseries
tsID_OU_URW = names(OU_URW)
model_results_OU_URW <- model_shift_results[tsID_OU_URW]

for (i in tsID_OU_URW) {
  OU_URW[[i]]$Shift_point <- model_shift_results[[i]][[12]]$parameters[["shift1"]]
}

#Splitting the model
OU_URW_subset1 = list()
OU_URW_subset2 = list()

if (length(OU_URW) > 0) {
  for (i in 1:length(OU_URW)) {
   gg <- rep(1:2, c(OU_URW[[i]]$Shift_point, length(OU_URW[[i]]$mm) - OU_URW[[i]]$Shift_point))
   OU_URW_split = paleoTS:::split4punc(OU_URW[[i]],gg, overlap=TRUE)
   OU_URW_subset1[[i]] = OU_URW_split[[1]]
   OU_URW_subset2[[i]] = OU_URW_split[[2]]
  }
}


### OU-Stasis ###
# Add a column with the best shift point to each timeseries
tsID_OU_Stasis = names(OU_Stasis)
model_results_OU_Stasis <- model_shift_results[tsID_OU_Stasis]

for (i in tsID_OU_Stasis) {
  OU_Stasis[[i]]$Shift_point <- model_shift_results[[i]][[13]]$parameters[["shift1"]]
}

#Splitting the model
OU_Stasis_subset1 = list()
OU_Stasis_subset2 = list()

if (length(OU_Stasis) > 0) {
  for (i in 1:length(OU_Stasis)) {
    gg <- rep(1:2, c(OU_Stasis[[i]]$Shift_point, length(OU_Stasis[[i]]$mm) - OU_Stasis[[i]]$Shift_point))
    OU_Stasis_split = paleoTS:::split4punc(OU_Stasis[[i]],gg, overlap=TRUE)
    OU_Stasis_subset1[[i]] = OU_Stasis_split[[1]]
    OU_Stasis_subset2[[i]] = OU_Stasis_split[[2]]
  }
}


### GRW-URW ###
# Add a column with the best shift point to each timeseries
tsID_GRW_URW = names(GRW_URW)
model_results_GRW_URW <- model_shift_results[tsID_GRW_URW]

for (i in tsID_GRW_URW) {
  GRW_URW[[i]]$Shift_point <- model_shift_results[[i]][[14]]$parameters[["shift1"]]
}

#Splitting the model
GRW_URW_subset1 = list()
GRW_URW_subset2 = list()

if (length(GRW_URW) > 0) {
  for (i in 1:length(GRW_URW)) {
   gg <- rep(1:2, c(GRW_URW[[i]]$Shift_point, length(GRW_URW[[i]]$mm) - GRW_URW[[i]]$Shift_point))
   GRW_URW_split = paleoTS:::split4punc(GRW_URW[[i]],gg, overlap=TRUE)
   GRW_URW_subset1[[i]] = GRW_URW_split[[1]]
   GRW_URW_subset2[[i]] = GRW_URW_split[[2]]
  }
}


### GRW-Stasis ###
# Add a column with the best shift point to each timeseries
tsID_GRW_Stasis = names(GRW_Stasis)
model_results_GRW_Stasis <- model_shift_results[tsID_GRW_Stasis]

for (i in tsID_GRW_Stasis) {
  GRW_Stasis[[i]]$Shift_point <- model_shift_results[[i]][[15]]$parameters[["shift1"]]
}

#Splitting the model
GRW_Stasis_subset1 = list()
GRW_Stasis_subset2 = list()

if (length(GRW_Stasis) > 0) {
  for (i in 1:length(GRW_Stasis)) {
    gg <- rep(1:2, c(GRW_Stasis[[i]]$Shift_point, length(GRW_Stasis[[i]]$mm) - GRW_Stasis[[i]]$Shift_point))
    GRW_Stasis_split = paleoTS:::split4punc(GRW_Stasis[[i]],gg, overlap=TRUE)
    GRW_Stasis_subset1[[i]] = GRW_Stasis_split[[1]]
    GRW_Stasis_subset2[[i]] = GRW_Stasis_split[[2]]
  }
}


### URW-Stasis ###
# Add a column with the best shift point to each timeseries
tsID_URW_Stasis = names(URW_Stasis)
model_results_URW_Stasis <- model_shift_results[tsID_URW_Stasis]

for (i in tsID_URW_Stasis) {
  URW_Stasis[[i]]$Shift_point <- model_shift_results[[i]][[16]]$parameters[["shift1"]]
}

#Splitting the model
URW_Stasis_subset1 = list()
URW_Stasis_subset2 = list()

if (length(URW_Stasis) > 0) {
  for (i in 1:length(URW_Stasis)) {
    gg <- rep(1:2, c(URW_Stasis[[i]]$Shift_point, length(URW_Stasis[[i]]$mm) - URW_Stasis[[i]]$Shift_point))
    URW_Stasis_split = paleoTS:::split4punc(URW_Stasis[[i]],gg, overlap=TRUE)
    URW_Stasis_subset1[[i]] = URW_Stasis_split[[1]]
    URW_Stasis_subset2[[i]] = URW_Stasis_split[[2]]
  }
}

# Saving OU subsets to work on them in another script (the adequacy of OU models is not working in parallel so need to be implemented in a loop)
# see assessing_models_shift_OU_adequacy.R for the code
save(OU, OU_mov_opt_anc, OU_mov_opt, Stasis_OU_subset2, URW_OU_subset2, GRW_OU_subset2,
     OU_OU_subset1, OU_OU_subset2, OU_GRW_subset1, OU_URW_subset1, OU_Stasis_subset1,
     file = "./results_paleoTS_v0.6.1/Results_OUsubsets_shiftmodels.RData")

#------------------------------------
# Testing the adequacy of the models
#------------------------------------

# Loading results of the OU adequacy (the OU models are not working in parallel)
load("./results_paleoTS_v0.6.1/Results_OUadeq_shiftmodels.RData")

# test adequacy
GRW_adeq <- mclapply(GRW, fit3adequacy.trend, plot = FALSE)
URW_adeq <- mclapply(URW, fit3adequacy.RW, plot = FALSE)
stasis_adeq <- mclapply(Stasis, fit4adequacy.stasis, plot = FALSE) 
strict_stasis_adeq <- mclapply(Strict_stasis, fit4adequacy.stasis, plot = FALSE)
decel_adeq <- mclapply(Decel, fit3adequacy.decel, plot = FALSE)
accel_adeq <- mclapply(Accel, fit3adequacy.RW, plot = FALSE)
#OU_adeq <- mclapply(OU, fit3adequacy.OU, plot = FALSE)
#OU_mov_opt_anc_adeq <- mclapply(OU_mov_opt_anc, fit3adequacy.OU, plot = FALSE)
#OU_mov_opt_adeq <- mclapply(OU_mov_opt, fit3adequacy.OU, plot = FALSE)

tsID_OU = names(OU)
names(OU_adeq) = tsID_OU

tsID_OU_mov_opt_anc = names(OU_mov_opt_anc)
names(OU_mov_opt_anc_adeq) = tsID_OU_mov_opt_anc

tsID_OU_mov_opt = names(OU_mov_opt)
names(OU_mov_opt) = tsID_OU_mov_opt

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
#Stasis_OU_subset2_adeq <- mclapply(Stasis_OU_subset2, fit3adequacy.OU, plot = FALSE)
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
#URW_OU_subset2_adeq <- mclapply(URW_OU_subset2, fit3adequacy.OU, plot = FALSE)
names(URW_OU_subset1_adeq) = tsID_URW_OU
names(URW_OU_subset2_adeq) = tsID_URW_OU

GRW_GRW_subset1_adeq <- mclapply(GRW_GRW_subset1, fit3adequacy.trend, plot = FALSE)
GRW_GRW_subset2_adeq <- mclapply(GRW_GRW_subset2, fit3adequacy.trend, plot = FALSE)
names(GRW_GRW_subset1_adeq) = tsID_GRW_GRW
names(GRW_GRW_subset2_adeq) = tsID_GRW_GRW

GRW_OU_subset1_adeq <- mclapply(GRW_OU_subset1, fit3adequacy.trend, plot = FALSE)
#GRW_OU_subset2_adeq <- mclapply(GRW_OU_subset2, fit3adequacy.OU, plot = FALSE)
names(GRW_OU_subset1_adeq) = tsID_GRW_OU
names(GRW_OU_subset2_adeq) = tsID_GRW_OU

#OU_OU_subset1_adeq <- mclapply(OU_OU_subset1, fit3adequacy.OU, plot = FALSE)
#OU_OU_subset2_adeq <- mclapply(OU_OU_subset2, fit3adequacy.OU, plot = FALSE)
names(OU_OU_subset1_adeq) = tsID_OU_OU 
names(OU_OU_subset2_adeq) = tsID_OU_OU

#OU_GRW_subset1_adeq <- mclapply(OU_GRW_subset1, fit3adequacy.OU, plot = FALSE)
OU_GRW_subset2_adeq <- mclapply(OU_GRW_subset2, fit3adequacy.trend, plot = FALSE)
names(OU_GRW_subset1_adeq) = tsID_OU_GRW
names(OU_GRW_subset2_adeq) = tsID_OU_GRW

#OU_URW_subset1_adeq <- mclapply(OU_URW_subset1, fit3adequacy.OU, plot = FALSE)
OU_URW_subset2_adeq <- mclapply(OU_URW_subset2, fit3adequacy.RW, plot = FALSE)
names(OU_URW_subset1_adeq) = tsID_OU_URW
names(OU_URW_subset2_adeq) = tsID_OU_URW

#OU_Stasis_subset1_adeq <- mclapply(OU_Stasis_subset1, fit3adequacy.OU, plot = FALSE)
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

                          
# get adequacy results for only adequate time series
GRW_adeq_passed <- adequate3tests(GRW_adeq)
URW_adeq_passed <- adequate3tests(URW_adeq)
stasis_adeq_passed <- adequate4tests(stasis_adeq)
strict_stasis_adeq_passed <- adequate4tests(strict_stasis_adeq)
decel_adeq_passed <- adequate3tests(decel_adeq)
accel_adeq_passed <- adequate3tests(accel_adeq)
OU_adeq_passed <- adequate2tests(OU_adeq)
OU_mov_opt_anc_adeq_passed <- adequate2tests(OU_mov_opt_anc_adeq)
OU_mov_opt_adeq_passed <- adequate2tests(OU_mov_opt_adeq)

Stasis_Stasis_subset1_adeq_passed <- adequate4tests(Stasis_Stasis_subset1_adeq)
Stasis_Stasis_subset2_adeq_passed <- adequate4tests(Stasis_Stasis_subset2_adeq)

Stasis_URW_subset1_adeq_passed <- adequate4tests(Stasis_URW_subset1_adeq)
Stasis_URW_subset2_adeq_passed <- adequate3tests(Stasis_URW_subset2_adeq)

Stasis_GRW_subset1_adeq_passed <- adequate4tests(Stasis_GRW_subset1_adeq)
Stasis_GRW_subset2_adeq_passed <- adequate3tests(Stasis_GRW_subset2_adeq)

Stasis_OU_subset1_adeq_passed <- adequate4tests(Stasis_OU_subset1_adeq)
Stasis_OU_subset2_adeq_passed <- adequate2tests(Stasis_OU_subset2_adeq)

URW_URW_subset1_adeq_passed <- adequate3tests(URW_URW_subset1_adeq)
URW_URW_subset2_adeq_passed <- adequate3tests(URW_URW_subset2_adeq)

URW_GRW_subset1_adeq_passed <- adequate3tests(URW_GRW_subset1_adeq)
URW_GRW_subset2_adeq_passed <- adequate3tests(URW_GRW_subset2_adeq)

URW_OU_subset1_adeq_passed <- adequate3tests(URW_OU_subset1_adeq)
URW_OU_subset2_adeq_passed <- adequate2tests(URW_OU_subset2_adeq)

GRW_GRW_subset1_adeq_passed <- adequate3tests(GRW_GRW_subset1_adeq)
GRW_GRW_subset2_adeq_passed <- adequate3tests(GRW_GRW_subset2_adeq)

GRW_OU_subset1_adeq_passed <- adequate3tests(GRW_OU_subset1_adeq)
GRW_OU_subset2_adeq_passed <- adequate2tests(GRW_OU_subset2_adeq)

OU_OU_subset1_adeq_passed <- adequate2tests(OU_OU_subset1_adeq)
OU_OU_subset2_adeq_passed <- adequate2tests(OU_OU_subset2_adeq)

OU_GRW_subset1_adeq_passed <- adequate2tests(OU_GRW_subset1_adeq)
OU_GRW_subset2_adeq_passed <- adequate3tests(OU_GRW_subset2_adeq)

OU_URW_subset1_adeq_passed <- adequate2tests(OU_URW_subset1_adeq)
OU_URW_subset2_adeq_passed <- adequate3tests(OU_URW_subset2_adeq)

OU_Stasis_subset1_adeq_passed <- adequate2tests(OU_Stasis_subset1_adeq)
OU_Stasis_subset2_adeq_passed <- adequate4tests(OU_Stasis_subset2_adeq)

GRW_URW_subset1_adeq_passed <- adequate3tests(GRW_URW_subset1_adeq)
GRW_URW_subset2_adeq_passed <- adequate3tests(GRW_URW_subset2_adeq)

GRW_Stasis_subset1_adeq_passed <- adequate3tests(GRW_Stasis_subset1_adeq)
GRW_Stasis_subset2_adeq_passed <- adequate4tests(GRW_Stasis_subset2_adeq)

URW_Stasis_subset1_adeq_passed <- adequate3tests(URW_Stasis_subset1_adeq)
URW_Stasis_subset2_adeq_passed <- adequate4tests(URW_Stasis_subset2_adeq)

                          
# merge split adequacy results of the time series if the two subsets passed the adequacy tests
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

# Get the list of the time series which passed the adequacy tests
GRW_list_adequate <- names(GRW_adeq_passed)
GRW_adeq_list <- GRW[names(GRW) %in% GRW_list_adequate]
URW_list_adequate <- names(URW_adeq_passed)
URW_adeq_list <- URW[names(URW) %in% URW_list_adequate]
stasis_list_adequate <- names(stasis_adeq_passed)
stasis_adeq_list <- Stasis[names(Stasis) %in% stasis_list_adequate]
strict_stasis_list_adequate <- names(strict_stasis_adeq_passed)
strict_stasis_adeq_list <- Strict_stasis[names(Strict_stasis) %in% strict_stasis_list_adequate]
decel_list_adequate <- names(decel_adeq_passed)
decel_adeq_list <- Decel[names(Decel) %in% decel_list_adequate]
accel_list_adequate <- names(accel_adeq_passed)
accel_adeq_list <- Accel[names(Accel) %in% accel_list_adequate]
OU_list_adequate <- names(OU_adeq_passed)
OU_adeq_list <- OU[names(OU) %in% OU_list_adequate]
OU_mov_opt_anc_list_adequate <- names(OU_mov_opt_anc_adeq_passed)
OU_mov_opt_anc_adeq_list <- OU_mov_opt_anc[names(OU_mov_opt_anc) %in% OU_mov_opt_anc_list_adequate]
OU_mov_opt_list_adequate <- names(OU_mov_opt_adeq_passed)
OU_mov_opt_adeq_list <- OU_mov_opt[names(OU_mov_opt) %in% OU_mov_opt_list_adequate]
Stasis_Stasis_list_adequate <- names(Stasis_Stasis_adeq_passed)
Stasis_Stasis_adeq_list <- Stasis_Stasis[names(Stasis_Stasis) %in% Stasis_Stasis_list_adequate]
Stasis_URW_list_adequate <- names(Stasis_URW_adeq_passed)
Stasis_URW_adeq_list <- Stasis_URW[names(Stasis_URW) %in% Stasis_URW_list_adequate]
Stasis_GRW_list_adequate <- names(Stasis_GRW_adeq_passed)
Stasis_GRW_adeq_list <- Stasis_GRW[names(Stasis_GRW) %in% Stasis_GRW_list_adequate]
Stasis_OU_list_adequate <- names(Stasis_OU_adeq_passed)
Stasis_OU_adeq_list <- Stasis_OU[names(Stasis_OU) %in% Stasis_OU_list_adequate]
URW_URW_list_adequate <- names(URW_URW_adeq_passed)
URW_URW_adeq_list <- URW_URW[names(URW_URW) %in% URW_URW_list_adequate]
URW_GRW_list_adequate <- names(URW_GRW_adeq_passed)
URW_GRW_adeq_list <- URW_GRW[names(URW_GRW) %in% URW_GRW_list_adequate]
URW_OU_list_adequate <- names(URW_OU_adeq_passed)
URW_OU_adeq_list <- URW_OU[names(URW_OU) %in% URW_OU_list_adequate]
GRW_GRW_list_adequate <- names(GRW_GRW_adeq_passed)
GRW_GRW_adeq_list <- GRW_GRW[names(GRW_GRW) %in% GRW_GRW_list_adequate]
GRW_OU_list_adequate <- names(GRW_OU_adeq_passed)
GRW_OU_adeq_list <- GRW_OU[names(GRW_OU) %in% GRW_OU_list_adequate]
OU_OU_list_adequate <- names(OU_OU_adeq_passed)
OU_OU_adeq_list <- OU_OU[names(OU_OU) %in% OU_OU_list_adequate]
OU_GRW_list_adequate <- names(OU_GRW_adeq_passed)
OU_GRW_adeq_list <- OU_GRW[names(OU_GRW) %in% OU_GRW_list_adequate]
OU_URW_list_adequate <- names(OU_URW_adeq_passed)
OU_URW_adeq_list <- OU_URW[names(OU_URW) %in% OU_URW_list_adequate]
OU_Stasis_list_adequate <- names(OU_Stasis_adeq_passed)
OU_Stasis_adeq_list <- OU_Stasis[names(OU_Stasis) %in% OU_Stasis_list_adequate]
GRW_URW_list_adequate <- names(GRW_URW_adeq_passed)
GRW_URW_adeq_list <- GRW_URW[names(GRW_URW) %in% GRW_URW_list_adequate]
GRW_Stasis_list_adequate <- names(GRW_Stasis_adeq_passed)
GRW_Stasis_adeq_list <- GRW_Stasis[names(GRW_Stasis) %in% GRW_Stasis_list_adequate]
URW_Stasis_list_adequate <- names(URW_Stasis_adeq_passed)
URW_Stasis_adeq_list <- URW_Stasis[names(URW_Stasis) %in% URW_Stasis_list_adequate]

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
Stasis_Stasis_c <- length(Stasis_Stasis_adeq_passed)
Stasis_URW_c <- length(Stasis_URW_adeq_passed)
Stasis_GRW_c <- length(Stasis_GRW_adeq_passed)
Stasis_OU_c <- length(Stasis_OU_adeq_passed)
URW_URW_c <- length(URW_URW_adeq_passed)
URW_GRW_c <- length(URW_GRW_adeq_passed)
URW_OU_c <- length(URW_OU_adeq_passed)
GRW_GRW_c <- length(GRW_GRW_adeq_passed)
GRW_OU_c <- length(GRW_OU_adeq_passed)
OU_OU_c <- length(OU_OU_adeq_passed)
OU_GRW_c <- length(OU_GRW_adeq_passed)
OU_URW_c <- length(OU_URW_adeq_passed)
OU_Stasis_c <- length(OU_Stasis_adeq_passed) 
GRW_URW_c <- length(GRW_URW_adeq_passed)
GRW_Stasis_c <- length(GRW_Stasis_adeq_passed)
URW_Stasis_c <- length(URW_Stasis_adeq_passed)

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
Stasis_Stasis_p <- (length(Stasis_Stasis_adeq_passed)/length(Stasis_Stasis_subset1_adeq))*100
Stasis_URW_p <- (length(Stasis_URW_adeq_passed)/length(Stasis_URW_subset1_adeq))*100
Stasis_GRW_p <- (length(Stasis_GRW_adeq_passed)/length(Stasis_GRW_subset1_adeq))*100
Stasis_OU_p <- (length(Stasis_OU_adeq_passed)/length(Stasis_OU_subset1_adeq))*100
URW_URW_p <- (length(URW_URW_adeq_passed)/length(URW_URW_subset1_adeq))*100
URW_GRW_p <- (length(URW_GRW_adeq_passed)/length(URW_GRW_subset1_adeq))*100
URW_OU_p <- (length(URW_OU_adeq_passed)/length(URW_OU_subset1_adeq))*100
GRW_GRW_p <- (length(GRW_GRW_adeq_passed)/length(GRW_GRW_subset1_adeq))*100
GRW_OU_p <- (length(GRW_OU_adeq_passed)/length(GRW_OU_subset1_adeq))*100
OU_OU_p <- (length(OU_OU_adeq_passed)/length(OU_OU_subset1_adeq))*100
OU_GRW_p <- (length(OU_GRW_adeq_passed)/length(OU_GRW_subset1_adeq))*100
OU_URW_p <- (length(OU_URW_adeq_passed)/length(OU_URW_subset1_adeq))*100
OU_Stasis_p <- (length(OU_Stasis_adeq_passed)/length(OU_Stasis_subset1_adeq))*100
GRW_URW_p <- (length(GRW_URW_adeq_passed)/length(GRW_URW_subset1_adeq))*100
GRW_Stasis_p <- (length(GRW_Stasis_adeq_passed)/length(GRW_Stasis_subset1_adeq))*100
URW_Stasis_p <- (length(URW_Stasis_adeq_passed)/length(URW_Stasis_subset1_adeq))*100

# make output table
adeq_table <- data.frame(
  model = c("GRW", "URW", "Stasis", "Strict Stasis", "Decel", "Accel", "OU",
            "OU Mov. Optm. (Ancestral State)", "OU Mov. Optm.", "Stasis-Stasis", 
            "Stasis-URW", "Stasis-GRW", "Stasis-OU", "URW-URW", "URW-GRW", "URW-OU",
            "GRW-GRW", "GRW-OU", "OU-OU", "OU-GRW", "OU-URW", "OU-Stasis", "GRW-URW",
            "GRW-Stasis", "URW-Stasis"),
  
  count_passed = c(GRW_c, URW_c, stasis_c, strict_stasis_c, decel_c, accel_c, OU_c, 
                   OU_mov_opt_anc_c, OU_mov_opt_c, Stasis_Stasis_c, Stasis_URW_c, Stasis_GRW_c, 
                   Stasis_OU_c, URW_URW_c, URW_GRW_c, URW_OU_c, GRW_GRW_c, GRW_OU_c, 
                   OU_OU_c, OU_GRW_c, OU_URW_c, OU_Stasis_c, GRW_URW_c, GRW_Stasis_c, 
                   URW_Stasis_c),
  
  percentage_passed = c(GRW_p, URW_p, stasis_p, strict_stasis_p, decel_p, accel_p, OU_p, 
                        OU_mov_opt_anc_p, OU_mov_opt_p, Stasis_Stasis_p, Stasis_URW_p, 
                        Stasis_GRW_p, Stasis_OU_p, URW_URW_p, URW_GRW_p, URW_OU_p, GRW_GRW_p, 
                        GRW_OU_p, OU_OU_p, OU_GRW_p, OU_URW_p, OU_Stasis_p, GRW_URW_p, 
                        GRW_Stasis_p, URW_Stasis_p)
)

Total_adeq_passed = sum(GRW_c, URW_c, stasis_c, strict_stasis_c, decel_c, accel_c, OU_c, 
    OU_mov_opt_anc_c, OU_mov_opt_c, Stasis_Stasis_c, Stasis_URW_c, Stasis_GRW_c, 
    Stasis_OU_c, URW_URW_c, URW_GRW_c, URW_OU_c, GRW_GRW_c, GRW_OU_c, 
    OU_OU_c, OU_GRW_c, OU_URW_c, OU_Stasis_c, GRW_URW_c, GRW_Stasis_c, 
    URW_Stasis_c)

# write to file
sink(file = "./results_paleoTS_v0.6.1/adequacy_passed_with_shift.txt")
adeq_table
paste("Total number of time series investigated:", length(model_noshift_results))
paste("Total number of time series which passed adequacy tests:", Total_adeq_passed)
paste("Percentage of time series which passed adequacy tests:", (Total_adeq_passed*100)/length(model_noshift_results))
sink()

# Save all the results
save.image(file='./results_paleoTS_v0.6.1/Results_fit_adequacy_shiftmodels.RData')


#----------------------------------------------------
# Table of global results for the models with shift
#----------------------------------------------------

# get the number of parameters for each model
K_noshift <- model_noshift_results[[1]]$K

model_shift_results_subsample <- model_shift_results[[1]]
model_shift_results_subsample2 <- list()

for (sublist in model_shift_results_subsample) { # remove the mirror models
  if (sublist$modelName %in% c("Strict Stasis", "Decel", "Accel", "OU",
                               "OU Mov. Optm. (Ancestral State)", "OU Mov. Optm.", "Stasis-Stasis", 
                               "Stasis-URW", "Stasis-GRW", "Stasis-OU", "URW-URW", "URW-GRW", "URW-OU",
                               "GRW-GRW", "GRW-OU", "OU-OU")) {
    model_shift_results_subsample2 <- append(model_shift_results_subsample2, list(sublist))
  }
}

K_shift <- sapply(model_shift_results_subsample2, function(item) item$K)

# Relative count - number of time series that are fitted the best by each model before adequacy
GRW_rc = length(GRW)
URW_rc = length(URW)
stasis_rc = length(Stasis)
strict_stasis_rc = length(Strict_stasis)
decel_rc = length(Decel)
accel_rc = length(Accel)
OU_rc = length(OU)
OU_mov_opt_anc_rc = length(OU_mov_opt_anc)
OU_mov_opt_rc = length(OU_mov_opt)
Stasis_Stasis_rc = length(Stasis_Stasis)
Stasis_URW_rc = length(Stasis_URW) + length(URW_Stasis)
Stasis_GRW_rc = length(Stasis_GRW) + length(GRW_Stasis)
Stasis_OU_rc = length(Stasis_OU) + length(OU_Stasis)
URW_URW_rc = length(URW_URW)
URW_GRW_rc = length(URW_GRW) + length(URW_GRW)
URW_OU_rc = length(URW_OU) + length(OU_URW)
GRW_GRW_rc = length(GRW_GRW)
GRW_OU_rc = length(GRW_OU) + length(OU_GRW)
OU_OU_rc = length(OU_OU) 

# Relative percentage - percentage of time series that are fitted the best by each model over the total number of time series before adequacy
GRW_rp = GRW_rc/length(model_noshift_results)*100
URW_rp = URW_rc/length(model_noshift_results)*100
stasis_rp = stasis_rc/length(model_noshift_results)*100
strict_stasis_rp = strict_stasis_rc/length(model_noshift_results)*100
decel_rp = decel_rc/length(model_noshift_results)*100
accel_rp = accel_rc/length(model_noshift_results)*100
OU_rp = OU_rc/length(model_noshift_results)*100
OU_mov_opt_anc_rp = OU_mov_opt_anc_rc/length(model_noshift_results)*100
OU_mov_opt_rp = OU_mov_opt_rc/length(model_noshift_results)*100
Stasis_Stasis_rp = Stasis_Stasis_rc/length(model_noshift_results)*100
Stasis_URW_rp = Stasis_URW_rc/length(model_noshift_results)*100
Stasis_GRW_rp = Stasis_GRW_rc/length(model_noshift_results)*100
Stasis_OU_rp = Stasis_OU_rc/length(model_noshift_results)*100
URW_URW_rp = URW_URW_rc/length(model_noshift_results)*100
URW_GRW_rp = URW_GRW_rc/length(model_noshift_results)*100
URW_OU_rp = URW_OU_rc/length(model_noshift_results)*100
GRW_GRW_rp = GRW_GRW_rc/length(model_noshift_results)*100
GRW_OU_rp = GRW_OU_rc/length(model_noshift_results)*100
OU_OU_rp = OU_OU_rc/length(model_noshift_results)*100

# compile results for the mirror models (stasis + URW = URW + stasis) for the percentage of time series which passed the adequacy test (over number of time series best described by this model)
Stasis_URW_p2 = (Stasis_URW_c + URW_Stasis_c)/(Stasis_URW_rc)*100
Stasis_GRW_p2 = (Stasis_GRW_c + GRW_Stasis_c)/(Stasis_GRW_rc)*100
Stasis_OU_p2 = (Stasis_OU_c + OU_Stasis_c)/(Stasis_OU_rc)*100
URW_GRW_p2 = (URW_GRW_c + GRW_URW_c)/(URW_GRW_rc)*100
URW_OU_p2 = (URW_OU_c + OU_URW_c)/(URW_OU_rc)*100
GRW_OU_p2 = (GRW_OU_c + OU_GRW_c)/(GRW_OU_rc)*100

# compile results for the mirror models (stasis + URW = URW + stasis) for the count of time series which passed the adequacy test
Stasis_URW_c2 = Stasis_URW_c + URW_Stasis_c
Stasis_GRW_c2 = Stasis_GRW_c + GRW_Stasis_c
Stasis_OU_c2 = Stasis_OU_c + OU_Stasis_c
URW_GRW_c2 = URW_GRW_c + GRW_URW_c
URW_OU_c2 = URW_OU_c + OU_URW_c
GRW_OU_c2 = GRW_OU_c + OU_GRW_c

# get absolute percentage of time series which passed the adequacy test over total of the time series 
GRW_ap <- (GRW_c/sum(aicc_results_complete$count))*100
URW_ap <- (URW_c/sum(aicc_results_complete$count))*100
stasis_ap <- (stasis_c/sum(aicc_results_complete$count))*100
strict_stasis_ap <- (strict_stasis_c/sum(aicc_results_complete$count))*100
decel_ap <- (decel_c/sum(aicc_results_complete$count))*100
accel_ap <- (accel_c/sum(aicc_results_complete$count))*100
OU_ap <- (OU_c/sum(aicc_results_complete$count))*100
OU_mov_opt_anc_ap <- (OU_mov_opt_anc_c/sum(aicc_results_complete$count))*100
OU_mov_opt_ap <- (OU_mov_opt_c/sum(aicc_results_complete$count))*100
Stasis_Stasis_ap <- (Stasis_Stasis_c/sum(aicc_results_complete$count))*100
Stasis_URW_ap <- (Stasis_URW_c2/sum(aicc_results_complete$count))*100
Stasis_GRW_ap <- (Stasis_GRW_c2/sum(aicc_results_complete$count))*100
Stasis_OU_ap <- (Stasis_OU_c2/sum(aicc_results_complete$count))*100
URW_URW_ap <- (URW_URW_c/sum(aicc_results_complete$count))*100
URW_GRW_ap <- (URW_GRW_c2/sum(aicc_results_complete$count))*100
URW_OU_ap <- (URW_OU_c2/sum(aicc_results_complete$count))*100
GRW_GRW_ap <- (GRW_GRW_c/sum(aicc_results_complete$count))*100
GRW_OU_ap <- (GRW_OU_c2/sum(aicc_results_complete$count))*100
OU_OU_ap <- (OU_OU_c/sum(aicc_results_complete$count))*100


# get average interval time for each model 
GRW_it <- list()
if (length(GRW_adeq_list) > 0) {
  for (i in 1:length(GRW_adeq_list)) {
    GRW_it[[i]] <- GRW_adeq_list[[i]]$tt[length(GRW_adeq_list[[i]]$tt)]
  }
}
GRW_it_values <- if (length(GRW_it) > 0) sapply(GRW_it, function(x) x[1]) else numeric(0)
GRW_it_mean <- if (length(GRW_it_values) > 0) mean(GRW_it_values) else NA
GRW_it_sd <- if (length(GRW_it_values) > 0) sd(GRW_it_values) else NA

URW_it <- list()
if (length(URW_adeq_list) > 0) {
  for (i in 1:length(URW_adeq_list)) {
    URW_it[[i]] <- URW_adeq_list[[i]]$tt[length(URW_adeq_list[[i]]$tt)]
  }
}
URW_it_values <- if (length(URW_it) > 0) sapply(URW_it, function(x) x[1]) else numeric(0)
URW_it_mean <- if (length(URW_it_values) > 0) mean(URW_it_values) else NA
URW_it_sd <- if (length(URW_it_values) > 0) sd(URW_it_values) else NA

Stasis_it <- list()
if (length(stasis_adeq_list) > 0) {
  for (i in 1:length(stasis_adeq_list)) {
    Stasis_it[[i]] <- stasis_adeq_list[[i]]$tt[length(stasis_adeq_list[[i]]$tt)]
  }
}
Stasis_it_values <- if (length(Stasis_it) > 0) sapply(Stasis_it, function(x) x[1]) else numeric(0)
Stasis_it_mean <- if (length(Stasis_it_values) > 0) mean(Stasis_it_values) else NA
Stasis_it_sd <- if (length(Stasis_it_values) > 0) sd(Stasis_it_values) else NA

Strict_stasis_it <- list()
if (length(strict_stasis_adeq_list) > 0) {
  for (i in 1:length(strict_stasis_adeq_list)) {
    Strict_stasis_it[[i]] <- strict_stasis_adeq_list[[i]]$tt[length(strict_stasis_adeq_list[[i]]$tt)]
  }
}
Strict_stasis_it_values <- if (length(Strict_stasis_it) > 0) sapply(Strict_stasis_it, function(x) x[1]) else numeric(0)
Strict_stasis_it_mean <- if (length(Strict_stasis_it_values) > 0) mean(Strict_stasis_it_values) else NA
Strict_stasis_it_sd <- if (length(Strict_stasis_it_values) > 0) sd(Strict_stasis_it_values) else NA

Decel_it <- list()
if (length(decel_adeq_list) > 0) {
  for (i in 1:length(decel_adeq_list)) {
    Decel_it[[i]] <- decel_adeq_list[[i]]$tt[length(decel_adeq_list[[i]]$tt)]
  }
}
Decel_it_values <- if (length(Decel_it) > 0) sapply(Decel_it, function(x) x[1]) else numeric(0)
Decel_it_mean <- if (length(Decel_it_values) > 0) mean(Decel_it_values) else NA
Decel_it_sd <- if (length(Decel_it_values) > 0) sd(Decel_it_values) else NA

Accel_it <- list()
if (length(accel_adeq_list) > 0) {
  for (i in 1:length(accel_adeq_list)) {
    Accel_it[[i]] <- accel_adeq_list[[i]]$tt[length(accel_adeq_list[[i]]$tt)]
  }
}
Accel_it_values <- if (length(Accel_it) > 0) sapply(Accel_it, function(x) x[1]) else numeric(0)
Accel_it_mean <- if (length(Accel_it_values) > 0) mean(Accel_it_values) else NA
Accel_it_sd <- if (length(Accel_it_values) > 0) sd(Accel_it_values) else NA

OU_it <- list()
if (length(OU_adeq_list) > 0) {
  for (i in 1:length(OU_adeq_list)) {
    OU_it[[i]] <- OU_adeq_list[[i]]$tt[length(OU_adeq_list[[i]]$tt)]
  }
}
OU_it_values <- if (length(OU_it) > 0) sapply(OU_it, function(x) x[1]) else numeric(0)
OU_it_mean <- if (length(OU_it_values) > 0) mean(OU_it_values) else NA
OU_it_sd <- if (length(OU_it_values) > 0) sd(OU_it_values) else NA

OU_mov_opt_anc_it <- list()
if (length(OU_mov_opt_anc_adeq_list) > 0) {
  for (i in 1:length(OU_mov_opt_anc_adeq_list)) {
    OU_mov_opt_anc_it[[i]] <- OU_mov_opt_anc_adeq_list[[i]]$tt[length(OU_mov_opt_anc_adeq_list[[i]]$tt)]
  }
}
OU_mov_opt_anc_it_values <- if (length(OU_mov_opt_anc_it) > 0) sapply(OU_mov_opt_anc_it, function(x) x[1]) else numeric(0)
OU_mov_opt_anc_it_mean <- if (length(OU_mov_opt_anc_it_values) > 0) mean(OU_mov_opt_anc_it_values) else NA
OU_mov_opt_anc_it_sd <- if (length(OU_mov_opt_anc_it_values) > 0) sd(OU_mov_opt_anc_it_values) else NA

OU_mov_opt_it <- list()
if (length(OU_mov_opt_adeq_list) > 0) {
  for (i in 1:length(OU_mov_opt_adeq_list)) {
    OU_mov_opt_it[[i]] <- OU_mov_opt_adeq_list[[i]]$tt[length(OU_mov_opt_adeq_list[[i]]$tt)]
  }
}
OU_mov_opt_it_values <- if (length(OU_mov_opt_it) > 0) sapply(OU_mov_opt_it, function(x) x[1]) else numeric(0)
OU_mov_opt_it_mean <- if (length(OU_mov_opt_it_values) > 0) mean(OU_mov_opt_it_values) else NA
OU_mov_opt_it_sd <- if (length(OU_mov_opt_it_values) > 0) sd(OU_mov_opt_it_values) else NA

Stasis_Stasis_it <- list()
if (length(Stasis_Stasis_adeq_list) > 0) {
  for (i in 1:length(Stasis_Stasis_adeq_list)) {
    Stasis_Stasis_it[[i]] <- Stasis_Stasis_adeq_list[[i]]$tt[length(Stasis_Stasis_adeq_list[[i]]$tt)]
  }
}
Stasis_Stasis_it_values <- if (length(Stasis_Stasis_it) > 0) sapply(Stasis_Stasis_it, function(x) x[1]) else numeric(0)
Stasis_Stasis_it_mean <- if (length(Stasis_Stasis_it_values) > 0) mean(Stasis_Stasis_it_values) else NA
Stasis_Stasis_it_sd <- if (length(Stasis_Stasis_it_values) > 0) sd(Stasis_Stasis_it_values) else NA

Stasis_URW_adeq_list_2 <- c(Stasis_URW_adeq_list, URW_Stasis_adeq_list)
Stasis_URW_it <- list()
if (length(Stasis_URW_adeq_list_2) > 0) {
  for (i in 1:length(Stasis_URW_adeq_list_2)) {
    Stasis_URW_it[[i]] <- Stasis_URW_adeq_list_2[[i]]$tt[length(Stasis_URW_adeq_list_2[[i]]$tt)]
  }
}
Stasis_URW_it_values <- if (length(Stasis_URW_it) > 0) sapply(Stasis_URW_it, function(x) x[1]) else numeric(0)
Stasis_URW_it_mean <- if (length(Stasis_URW_it_values) > 0) mean(Stasis_URW_it_values) else NA
Stasis_URW_it_sd <- if (length(Stasis_URW_it_values) > 0) sd(Stasis_URW_it_values) else NA

Stasis_GRW_adeq_list_2 <- c(Stasis_GRW_adeq_list, GRW_Stasis_adeq_list)
Stasis_GRW_it <- list()
if (length(Stasis_GRW_adeq_list_2) > 0) {
  for (i in 1:length(Stasis_GRW_adeq_list_2)) {
    Stasis_GRW_it[[i]] <- Stasis_GRW_adeq_list_2[[i]]$tt[length(Stasis_GRW_adeq_list_2[[i]]$tt)]
  }
}
Stasis_GRW_it_values <- if (length(Stasis_GRW_it) > 0) sapply(Stasis_GRW_it, function(x) x[1]) else numeric(0)
Stasis_GRW_it_mean <- if (length(Stasis_GRW_it_values) > 0) mean(Stasis_GRW_it_values) else NA
Stasis_GRW_it_sd <- if (length(Stasis_GRW_it_values) > 0) sd(Stasis_GRW_it_values) else NA

Stasis_OU_adeq_list_2 <- c(Stasis_OU_adeq_list, OU_Stasis_adeq_list)
Stasis_OU_it <- list()
if (length(Stasis_OU_adeq_list_2) > 0) {
  for (i in 1:length(Stasis_OU_adeq_list_2)) {
    Stasis_OU_it[[i]] <- Stasis_OU_adeq_list_2[[i]]$tt[length(Stasis_OU_adeq_list_2[[i]]$tt)]
  }
}
Stasis_OU_it_values <- if (length(Stasis_OU_it) > 0) sapply(Stasis_OU_it, function(x) x[1]) else numeric(0)
Stasis_OU_it_mean <- if (length(Stasis_OU_it_values) > 0) mean(Stasis_OU_it_values) else NA
Stasis_OU_it_sd <- if (length(Stasis_OU_it_values) > 0) sd(Stasis_OU_it_values) else NA

URW_URW_it <- list()
if (length(URW_URW_adeq_list) > 0) {
  for (i in 1:length(URW_URW_adeq_list)) {
    URW_URW_it[[i]] <- URW_URW_adeq_list[[i]]$tt[length(URW_URW_adeq_list[[i]]$tt)]
  }
}
URW_URW_it_values <- if (length(URW_URW_it) > 0) sapply(URW_URW_it, function(x) x[1]) else numeric(0)
URW_URW_it_mean <- if (length(URW_URW_it_values) > 0) mean(URW_URW_it_values) else NA
URW_URW_it_sd <- if (length(URW_URW_it_values) > 0) sd(URW_URW_it_values) else NA

URW_GRW_adeq_list_2 <- c(URW_GRW_adeq_list, GRW_URW_adeq_list)
URW_GRW_it <- list()
if (length(URW_GRW_adeq_list_2) > 0) {
  for (i in 1:length(URW_GRW_adeq_list_2)) {
    URW_GRW_it[[i]] <- URW_GRW_adeq_list_2[[i]]$tt[length(URW_GRW_adeq_list_2[[i]]$tt)]
  }
}
URW_GRW_it_values <- if (length(URW_GRW_it) > 0) sapply(URW_GRW_it, function(x) x[1]) else numeric(0)
URW_GRW_it_mean <- if (length(URW_GRW_it_values) > 0) mean(URW_GRW_it_values) else NA
URW_GRW_it_sd <- if (length(URW_GRW_it_values) > 0) sd(URW_GRW_it_values) else NA

URW_OU_adeq_list_2 <- c(URW_OU_adeq_list, OU_URW_adeq_list)
URW_OU_it <- list()
if (length(URW_OU_adeq_list_2) > 0) {
  for (i in 1:length(URW_OU_adeq_list_2)) {
    URW_OU_it[[i]] <- URW_OU_adeq_list_2[[i]]$tt[length(URW_OU_adeq_list_2[[i]]$tt)]
  }
}
URW_OU_it_values <- if (length(URW_OU_it) > 0) sapply(URW_OU_it, function(x) x[1]) else numeric(0)
URW_OU_it_mean <- if (length(URW_OU_it_values) > 0) mean(URW_OU_it_values) else NA
URW_OU_it_sd <- if (length(URW_OU_it_values) > 0) sd(URW_OU_it_values) else NA

GRW_GRW_it <- list()
if (length(GRW_GRW_adeq_list) > 0) {
  for (i in 1:length(GRW_GRW_adeq_list)) {
    GRW_GRW_it[[i]] <- GRW_GRW_adeq_list[[i]]$tt[length(GRW_GRW_adeq_list[[i]]$tt)]
  }
}
GRW_GRW_it_values <- if (length(GRW_GRW_it) > 0) sapply(GRW_GRW_it, function(x) x[1]) else numeric(0)
GRW_GRW_it_mean <- if (length(GRW_GRW_it_values) > 0) mean(GRW_GRW_it_values) else NA
GRW_GRW_it_sd <- if (length(GRW_GRW_it_values) > 0) sd(GRW_GRW_it_values) else NA

GRW_OU_adeq_list_2 <- c(GRW_OU_adeq_list, OU_GRW_adeq_list)
GRW_OU_it <- list()
if (length(GRW_OU_adeq_list_2) > 0) {
  for (i in 1:length(GRW_OU_adeq_list_2)) {
    GRW_OU_it[[i]] <- GRW_OU_adeq_list_2[[i]]$tt[length(GRW_OU_adeq_list_2[[i]]$tt)]
  }
}
GRW_OU_it_values <- if (length(GRW_OU_it) > 0) sapply(GRW_OU_it, function(x) x[1]) else numeric(0)
GRW_OU_it_mean <- if (length(GRW_OU_it_values) > 0) mean(GRW_OU_it_values) else NA
GRW_OU_it_sd <- if (length(GRW_OU_it_values) > 0) sd(GRW_OU_it_values) else NA

OU_OU_it <- list()
if (length(OU_OU_adeq_list) > 0) {
  for (i in 1:length(OU_OU_adeq_list)) {
    OU_OU_it[[i]] <- OU_OU_adeq_list[[i]]$tt[length(OU_OU_adeq_list[[i]]$tt)]
  }
}
OU_OU_it_values <- if (length(OU_OU_it) > 0) sapply(OU_OU_it, function(x) x[1]) else numeric(0)
OU_OU_it_mean <- if (length(OU_OU_it_values) > 0) mean(OU_OU_it_values) else NA
OU_OU_it_sd <- if (length(OU_OU_it_values) > 0) sd(OU_OU_it_values) else NA

# get percentage of traits in the size and shape category for each model
metadata_trait <- read_delim("./timeseries/metadata_trait.txt", col_names = TRUE, delim = "\t")

metadata_trait_organized <- metadata_trait %>% # remove data with trait category = number
  filter(!trait_category %in% "number") 

metadata_trait_organized <- metadata_trait_organized %>% # classify trait_category in shape or size 
  mutate(size_or_shape = case_when(
    trait_category %in% c("size", "area") ~ "size",    # Assign "size" if ll is "sh"
    trait_category %in% c("shape", "angle") ~ "shape",))

metadata_trait_organized <- metadata_trait_organized %>% # remove lines without shape or size
  filter(!is.na(size_or_shape))

save(metadata_trait_organized, file = "./timeseries/metadata_trait_2.RData") # Save the new metadata_trait file with the column shape_or_size

GRW_size_nbr = 0
GRW_shape_nbr = 0
if (length(GRW_adeq_list) > 0) {
for (i in 1:length(GRW_adeq_list)) {
  for (j in 1:nrow(metadata_trait_organized)) {
    if (names(GRW_adeq_list[i]) == metadata_trait_organized$tsID[[j]]) {
      if (metadata_trait_organized$size_or_shape[[j]] == "size") {
        GRW_adeq_list[[i]]$trait_type <- "size"
        GRW_size_nbr = GRW_size_nbr + 1
      } else if (metadata_trait_organized$size_or_shape[[j]] == "shape") {
        GRW_adeq_list[[i]]$trait_type <- "shape"
        GRW_shape_nbr = GRW_shape_nbr + 1
      }
    }
  }
}
}

URW_size_nbr = 0
URW_shape_nbr = 0
if (length(URW_adeq_list) > 0) {
for (i in 1:length(URW_adeq_list)) {
  for (j in 1:nrow(metadata_trait_organized)) {
    if (names(URW_adeq_list[i]) == metadata_trait_organized$tsID[[j]]) {
      if (metadata_trait_organized$size_or_shape[[j]] == "size") {
        URW_adeq_list[[i]]$trait_type <- "size"
        URW_size_nbr = URW_size_nbr + 1
      } else if (metadata_trait_organized$size_or_shape[[j]] == "shape") {
        URW_adeq_list[[i]]$trait_type <- "shape"
        URW_shape_nbr = URW_shape_nbr + 1
      }
    }
  }
}
}

stasis_size_nbr = 0
stasis_shape_nbr = 0
if (length(stasis_adeq_list) > 0) {
for (i in 1:length(stasis_adeq_list)) {
  for (j in 1:nrow(metadata_trait_organized)) {
    if (names(stasis_adeq_list[i]) == metadata_trait_organized$tsID[[j]]) {
      if (metadata_trait_organized$size_or_shape[[j]] == "size") {
        stasis_adeq_list[[i]]$trait_type <- "size"
        stasis_size_nbr = stasis_size_nbr + 1
      } else if (metadata_trait_organized$size_or_shape[[j]] == "shape") {
        stasis_adeq_list[[i]]$trait_type <- "shape"
        stasis_shape_nbr = stasis_shape_nbr + 1
      }
    }
  }
}
}

strict_stasis_size_nbr = 0
strict_stasis_shape_nbr = 0
if (length(strict_stasis_adeq_list) > 0) {
for (i in 1:length(strict_stasis_adeq_list)) {
  for (j in 1:nrow(metadata_trait_organized)) {
    if (names(strict_stasis_adeq_list[i]) == metadata_trait_organized$tsID[[j]]) {
      if (metadata_trait_organized$size_or_shape[[j]] == "size") {
        strict_stasis_adeq_list[[i]]$trait_type <- "size"
        strict_stasis_size_nbr = strict_stasis_size_nbr + 1
      } else if (metadata_trait_organized$size_or_shape[[j]] == "shape") {
        strict_stasis_adeq_list[[i]]$trait_type <- "shape"
        strict_stasis_shape_nbr = strict_stasis_shape_nbr + 1
      }
    }
  }
}
}

decel_size_nbr = 0
decel_shape_nbr = 0
if (length(decel_adeq_list) > 0) {
for (i in 1:length(decel_adeq_list)) {
  for (j in 1:nrow(metadata_trait_organized)) {
    if (names(decel_adeq_list[i]) == metadata_trait_organized$tsID[[j]]) {
      if (metadata_trait_organized$size_or_shape[[j]] == "size") {
        decel_adeq_list[[i]]$trait_type <- "size"
        decel_size_nbr = decel_size_nbr + 1
      } else if (metadata_trait_organized$size_or_shape[[j]] == "shape") {
        decel_adeq_list[[i]]$trait_type <- "shape"
        decel_shape_nbr = decel_shape_nbr + 1
      }
    }
  }
}
}

accel_size_nbr = 0
accel_shape_nbr = 0
if (length(accel_adeq_list) > 0) {
for (i in 1:length(accel_adeq_list)) {
  for (j in 1:nrow(metadata_trait_organized)) {
    if (names(accel_adeq_list[i]) == metadata_trait_organized$tsID[[j]]) {
      if (metadata_trait_organized$size_or_shape[[j]] == "size") {
        accel_adeq_list[[i]]$trait_type <- "size"
        accel_size_nbr = accel_size_nbr + 1
      } else if (metadata_trait_organized$size_or_shape[[j]] == "shape") {
        accel_adeq_list[[i]]$trait_type <- "shape"
        accel_shape_nbr = accel_shape_nbr + 1
      }
    }
  }
}
}

OU_size_nbr = 0
OU_shape_nbr = 0
if (length(OU_adeq_list) > 0) {
for (i in 1:length(OU_adeq_list)) {
  for (j in 1:nrow(metadata_trait_organized)) {
    if (names(OU_adeq_list[i]) == metadata_trait_organized$tsID[[j]]) {
      if (metadata_trait_organized$size_or_shape[[j]] == "size") {
        OU_adeq_list[[i]]$trait_type <- "size"
        OU_size_nbr = OU_size_nbr + 1
      } else if (metadata_trait_organized$size_or_shape[[j]] == "shape") {
        OU_adeq_list[[i]]$trait_type <- "shape"
        OU_shape_nbr = OU_shape_nbr + 1
      }
    }
  }
}
}

OU_mov_opt_anc_size_nbr = 0
OU_mov_opt_anc_shape_nbr = 0
if (length(OU_mov_opt_anc_adeq_list) > 0) {
for (i in 1:length(OU_mov_opt_anc_adeq_list)) {
  for (j in 1:nrow(metadata_trait_organized)) {
    if (names(OU_mov_opt_anc_adeq_list[i]) == metadata_trait_organized$tsID[[j]]) {
      if (metadata_trait_organized$size_or_shape[[j]] == "size") {
        OU_mov_opt_anc_adeq_list[[i]]$trait_type <- "size"
        OU_mov_opt_anc_size_nbr = OU_mov_opt_anc_size_nbr + 1
      } else if (metadata_trait_organized$size_or_shape[[j]] == "shape") {
        OU_mov_opt_anc_adeq_list[[i]]$trait_type <- "shape"
        OU_mov_opt_anc_shape_nbr = OU_mov_opt_anc_shape_nbr + 1
      }
    }
  }
}
}

OU_mov_opt_size_nbr = 0
OU_mov_opt_shape_nbr = 0
if (length(OU_mov_opt_adeq_list) > 0) {
for (i in 1:length(OU_mov_opt_adeq_list)) {
  for (j in 1:nrow(metadata_trait_organized)) {
    if (names(OU_mov_opt_adeq_list[i]) == metadata_trait_organized$tsID[[j]]) {
      if (metadata_trait_organized$size_or_shape[[j]] == "size") {
        OU_mov_opt_adeq_list[[i]]$trait_type <- "size"
        OU_mov_opt_size_nbr = OU_mov_opt_size_nbr + 1
      } else if (metadata_trait_organized$size_or_shape[[j]] == "shape") {
        OU_mov_opt_adeq_list[[i]]$trait_type <- "shape"
        OU_mov_opt_shape_nbr = OU_mov_opt_shape_nbr + 1
      }
    }
  }
}
}

Stasis_Stasis_size_nbr = 0
Stasis_Stasis_shape_nbr = 0
if (length(Stasis_Stasis_adeq_list) > 0) {
  for (i in 1:length(Stasis_Stasis_adeq_list)) {
    for (j in 1:nrow(metadata_trait_organized)) {
      if (names(Stasis_Stasis_adeq_list[i]) == metadata_trait_organized$tsID[[j]]) {
        if (metadata_trait_organized$size_or_shape[[j]] == "size") {
          Stasis_Stasis_adeq_list[[i]]$trait_type <- "size"
          Stasis_Stasis_size_nbr = Stasis_Stasis_size_nbr + 1
        } else if (metadata_trait_organized$size_or_shape[[j]] == "shape") {
          Stasis_Stasis_adeq_list[[i]]$trait_type <- "shape"
          Stasis_Stasis_shape_nbr = Stasis_Stasis_shape_nbr + 1
        }
      }
    }
  }
}

Stasis_URW_size_nbr = 0
Stasis_URW_shape_nbr = 0
if (length(Stasis_URW_adeq_list) > 0) {
  for (i in 1:length(Stasis_URW_adeq_list)) {
    for (j in 1:nrow(metadata_trait_organized)) {
      if (names(Stasis_URW_adeq_list[i]) == metadata_trait_organized$tsID[[j]]) {
        if (metadata_trait_organized$size_or_shape[[j]] == "size") {
          Stasis_URW_adeq_list[[i]]$trait_type <- "size"
          Stasis_URW_size_nbr = Stasis_URW_size_nbr + 1
        } else if (metadata_trait_organized$size_or_shape[[j]] == "shape") {
          Stasis_URW_adeq_list[[i]]$trait_type <- "shape"
          Stasis_URW_shape_nbr = Stasis_URW_shape_nbr + 1
        }
      }
    }
  }
}

Stasis_GRW_size_nbr = 0
Stasis_GRW_shape_nbr = 0
if (length(Stasis_GRW_adeq_list) > 0) {
  for (i in 1:length(Stasis_GRW_adeq_list)) {
    for (j in 1:nrow(metadata_trait_organized)) {
      if (names(Stasis_GRW_adeq_list[i]) == metadata_trait_organized$tsID[[j]]) {
        if (metadata_trait_organized$size_or_shape[[j]] == "size") {
          Stasis_GRW_adeq_list[[i]]$trait_type <- "size"
          Stasis_GRW_size_nbr = Stasis_GRW_size_nbr + 1
        } else if (metadata_trait_organized$size_or_shape[[j]] == "shape") {
          Stasis_GRW_adeq_list[[i]]$trait_type <- "shape"
          Stasis_GRW_shape_nbr = Stasis_GRW_shape_nbr + 1
        }
      }
    }
  }
}

Stasis_OU_size_nbr = 0
Stasis_OU_shape_nbr = 0
if (length(Stasis_OU_adeq_list) > 0) {
  for (i in 1:length(Stasis_OU_adeq_list)) {
    for (j in 1:nrow(metadata_trait_organized)) {
      if (names(Stasis_OU_adeq_list[i]) == metadata_trait_organized$tsID[[j]]) {
        if (metadata_trait_organized$size_or_shape[[j]] == "size") {
          Stasis_OU_adeq_list[[i]]$trait_type <- "size"
          Stasis_OU_size_nbr = Stasis_OU_size_nbr + 1
        } else if (metadata_trait_organized$size_or_shape[[j]] == "shape") {
          Stasis_OU_adeq_list[[i]]$trait_type <- "shape"
          Stasis_OU_shape_nbr = Stasis_OU_shape_nbr + 1
        }
      }
    }
  }
}

URW_URW_size_nbr = 0
URW_URW_shape_nbr = 0
if (length(URW_URW_adeq_list) > 0) {
  for (i in 1:length(URW_URW_adeq_list)) {
    for (j in 1:nrow(metadata_trait_organized)) {
      if (names(URW_URW_adeq_list[i]) == metadata_trait_organized$tsID[[j]]) {
        if (metadata_trait_organized$size_or_shape[[j]] == "size") {
          URW_URW_adeq_list[[i]]$trait_type <- "size"
          URW_URW_size_nbr = URW_URW_size_nbr + 1
        } else if (metadata_trait_organized$size_or_shape[[j]] == "shape") {
          URW_URW_adeq_list[[i]]$trait_type <- "shape"
          URW_URW_shape_nbr = URW_URW_shape_nbr + 1
        }
      }
    }
  }
}

URW_GRW_size_nbr = 0
URW_GRW_shape_nbr = 0
if (length(URW_GRW_adeq_list) > 0) {
  for (i in 1:length(URW_GRW_adeq_list)) {
    for (j in 1:nrow(metadata_trait_organized)) {
      if (names(URW_GRW_adeq_list[i]) == metadata_trait_organized$tsID[[j]]) {
        if (metadata_trait_organized$size_or_shape[[j]] == "size") {
          URW_GRW_adeq_list[[i]]$trait_type <- "size"
          URW_GRW_size_nbr = URW_GRW_size_nbr + 1
        } else if (metadata_trait_organized$size_or_shape[[j]] == "shape") {
          URW_GRW_adeq_list[[i]]$trait_type <- "shape"
          URW_GRW_shape_nbr = URW_GRW_shape_nbr + 1
        }
      }
    }
  }
}

URW_OU_size_nbr = 0
URW_OU_shape_nbr = 0
if (length(URW_OU_adeq_list) > 0) {
  for (i in 1:length(URW_OU_adeq_list)) {
    for (j in 1:nrow(metadata_trait_organized)) {
      if (names(URW_OU_adeq_list[i]) == metadata_trait_organized$tsID[[j]]) {
        if (metadata_trait_organized$size_or_shape[[j]] == "size") {
          URW_OU_adeq_list[[i]]$trait_type <- "size"
          URW_OU_size_nbr = URW_OU_size_nbr + 1
        } else if (metadata_trait_organized$size_or_shape[[j]] == "shape") {
          URW_OU_adeq_list[[i]]$trait_type <- "shape"
          URW_OU_shape_nbr = URW_OU_shape_nbr + 1
        }
      }
    }
  }
}

GRW_GRW_size_nbr = 0
GRW_GRW_shape_nbr = 0
if (length(GRW_GRW_adeq_list) > 0) {
  for (i in 1:length(GRW_GRW_adeq_list)) {
    for (j in 1:nrow(metadata_trait_organized)) {
      if (names(GRW_GRW_adeq_list[i]) == metadata_trait_organized$tsID[[j]]) {
        if (metadata_trait_organized$size_or_shape[[j]] == "size") {
          GRW_GRW_adeq_list[[i]]$trait_type <- "size"
          GRW_GRW_size_nbr = GRW_GRW_size_nbr + 1
        } else if (metadata_trait_organized$size_or_shape[[j]] == "shape") {
          GRW_GRW_adeq_list[[i]]$trait_type <- "shape"
          GRW_GRW_shape_nbr = GRW_GRW_shape_nbr + 1
        }
      }
    }
  }
}

GRW_OU_size_nbr = 0
GRW_OU_shape_nbr = 0
if (length(GRW_OU_adeq_list) > 0) {
  for (i in 1:length(GRW_OU_adeq_list)) {
    for (j in 1:nrow(metadata_trait_organized)) {
      if (names(GRW_OU_adeq_list[i]) == metadata_trait_organized$tsID[[j]]) {
        if (metadata_trait_organized$size_or_shape[[j]] == "size") {
          GRW_OU_adeq_list[[i]]$trait_type <- "size"
          GRW_OU_size_nbr = GRW_OU_size_nbr + 1
        } else if (metadata_trait_organized$size_or_shape[[j]] == "shape") {
          GRW_OU_adeq_list[[i]]$trait_type <- "shape"
          GRW_OU_shape_nbr = GRW_OU_shape_nbr + 1
        }
      }
    }
  }
}

OU_OU_size_nbr = 0
OU_OU_shape_nbr = 0
if (length(OU_OU_adeq_list) > 0) {
  for (i in 1:length(OU_OU_adeq_list)) {
    for (j in 1:nrow(metadata_trait_organized)) {
      if (names(OU_OU_adeq_list[i]) == metadata_trait_organized$tsID[[j]]) {
        if (metadata_trait_organized$size_or_shape[[j]] == "size") {
          OU_OU_adeq_list[[i]]$trait_type <- "size"
          OU_OU_size_nbr = OU_OU_size_nbr + 1
        } else if (metadata_trait_organized$size_or_shape[[j]] == "shape") {
          OU_OU_adeq_list[[i]]$trait_type <- "shape"
          OU_OU_shape_nbr = OU_OU_shape_nbr + 1
        }
      }
    }
  }
}

# percentage of time series being size traits that are best fitted by the model over the total number of time series best fitted by the model
GRW_size_p = GRW_size_nbr/length(GRW_adeq_list)*100
URW_size_p = URW_size_nbr/length(URW_adeq_list)*100
stasis_size_p = stasis_size_nbr/length(stasis_adeq_list)*100
strict_stasis_size_p = strict_stasis_size_nbr/length(strict_stasis_adeq_list)*100
decel_size_p = decel_size_nbr/length(decel_adeq_list)*100
accel_size_p = accel_size_nbr/length(accel_adeq_list)*100
OU_size_p = OU_size_nbr/length(OU_adeq_list)*100
OU_mov_opt_anc_size_p = OU_mov_opt_anc_size_nbr/length(OU_mov_opt_anc_adeq_list)*100
OU_mov_opt_size_p = OU_mov_opt_size_nbr/length(OU_mov_opt_adeq_list)*100
Stasis_Stasis_size_p = Stasis_Stasis_size_nbr/length(Stasis_Stasis_adeq_list)*100
Stasis_URW_size_p = Stasis_URW_size_nbr/length(Stasis_URW_adeq_list)*100
Stasis_GRW_size_p = Stasis_GRW_size_nbr/length(Stasis_GRW_adeq_list)*100
Stasis_OU_size_p = Stasis_OU_size_nbr/length(Stasis_OU_adeq_list)*100
URW_URW_size_p = URW_URW_size_nbr/length(URW_URW_adeq_list)*100
URW_GRW_size_p = URW_GRW_size_nbr/length(URW_GRW_adeq_list)*100
URW_OU_size_p = URW_OU_size_nbr/length(URW_OU_adeq_list)*100
GRW_GRW_size_p = GRW_GRW_size_nbr/length(GRW_GRW_adeq_list)*100
GRW_OU_size_p = GRW_OU_size_nbr/length(GRW_OU_adeq_list)*100
OU_OU_size_p = OU_OU_size_nbr/length(OU_OU_adeq_list)*100

# percentage of time series being shape traits that are best fitted by the model over the total number of time series best fitted by the model
GRW_shape_p = GRW_shape_nbr/length(GRW_adeq_list)*100
URW_shape_p = URW_shape_nbr/length(URW_adeq_list)*100
stasis_shape_p = stasis_shape_nbr/length(stasis_adeq_list)*100
strict_stasis_shape_p = strict_stasis_shape_nbr/length(strict_stasis_adeq_list)*100
decel_shape_p = decel_shape_nbr/length(decel_adeq_list)*100
accel_shape_p = accel_shape_nbr/length(accel_adeq_list)*100
OU_shape_p = OU_shape_nbr/length(OU_adeq_list)*100
OU_mov_opt_anc_shape_p = OU_mov_opt_anc_shape_nbr/length(OU_mov_opt_anc_adeq_list)*100
OU_mov_opt_shape_p = OU_mov_opt_shape_nbr/length(OU_mov_opt_adeq_list)*100
Stasis_Stasis_shape_p = Stasis_Stasis_shape_nbr/length(Stasis_Stasis_adeq_list)*100
Stasis_URW_shape_p = Stasis_URW_shape_nbr/length(Stasis_URW_adeq_list)*100
Stasis_GRW_shape_p = Stasis_GRW_shape_nbr/length(Stasis_GRW_adeq_list)*100
Stasis_OU_shape_p = Stasis_OU_shape_nbr/length(Stasis_OU_adeq_list)*100
URW_URW_shape_p = URW_URW_shape_nbr/length(URW_URW_adeq_list)*100
URW_GRW_shape_p = URW_GRW_shape_nbr/length(URW_GRW_adeq_list)*100
URW_OU_shape_p = URW_OU_shape_nbr/length(URW_OU_adeq_list)*100
GRW_GRW_shape_p = GRW_GRW_shape_nbr/length(GRW_GRW_adeq_list)*100
GRW_OU_shape_p = GRW_OU_shape_nbr/length(GRW_OU_adeq_list)*100
OU_OU_shape_p = OU_OU_shape_nbr/length(OU_OU_adeq_list)*100

# make output table
resultshift_table <- data.frame(
  
  model = c("GRW", "URW", "Stasis", "Strict Stasis", "Decel", "Accel", "OU",
            "OU Mov. Optm. (Ancestral State)", "OU Mov. Optm.", "Stasis-Stasis", 
            "Stasis-URW", "Stasis-GRW", "Stasis-OU", "URW-URW", "URW-GRW", "URW-OU",
            "GRW-GRW", "GRW-OU", "OU-OU"),
  
  K = c(K_noshift, K_shift),
    
  relative_count = c(GRW_rc, URW_rc, stasis_rc, strict_stasis_rc, decel_rc, accel_rc, OU_rc, 
                       OU_mov_opt_anc_rc, OU_mov_opt_rc, Stasis_Stasis_rc, Stasis_URW_rc, Stasis_GRW_rc, 
                       Stasis_OU_rc, URW_URW_rc, URW_GRW_rc, URW_OU_rc, GRW_GRW_rc, GRW_OU_rc, 
                       OU_OU_rc),
  
  relative_percentage = c(GRW_rp, URW_rp, stasis_rp, strict_stasis_rp, decel_rp, accel_rp, OU_rp, 
                          OU_mov_opt_anc_rp, OU_mov_opt_rp, Stasis_Stasis_rp, Stasis_URW_rp, Stasis_GRW_rp, 
                          Stasis_OU_rp, URW_URW_rp, URW_GRW_rp, URW_OU_rp, GRW_GRW_rp, GRW_OU_rp, 
                          OU_OU_rp),
  
  percentage_passed = c(GRW_p, URW_p, stasis_p, strict_stasis_p, decel_p, accel_p, OU_p, 
                         OU_mov_opt_anc_p, OU_mov_opt_p, Stasis_Stasis_p, Stasis_URW_p2, Stasis_GRW_p2, 
                         Stasis_OU_p2, URW_URW_p, URW_GRW_p2, URW_OU_p2, GRW_GRW_p, GRW_OU_p2, 
                         OU_OU_p),
  
  absolute_count = c(GRW_c, URW_c, stasis_c, strict_stasis_c, decel_c, accel_c, OU_c, 
                    OU_mov_opt_anc_c, OU_mov_opt_c, Stasis_Stasis_c, Stasis_URW_c2, Stasis_GRW_c2, 
                    Stasis_OU_c2, URW_URW_c, URW_GRW_c2, URW_OU_c2, GRW_GRW_c, GRW_OU_c2, 
                    OU_OU_c),
  
  absolute_percentage = c(GRW_ap, URW_ap, stasis_ap, strict_stasis_ap, decel_ap, accel_ap, OU_ap, 
                          OU_mov_opt_anc_ap, OU_mov_opt_ap, Stasis_Stasis_ap, Stasis_URW_ap, Stasis_GRW_ap, 
                          Stasis_OU_ap, URW_URW_ap, URW_GRW_ap, URW_OU_ap, GRW_GRW_ap, GRW_OU_ap, 
                          OU_OU_ap),
  
  interval_time_passed = c(GRW_it_mean, URW_it_mean, Stasis_it_mean, Strict_stasis_it_mean, Decel_it_mean, Accel_it_mean, OU_it_mean, 
                    OU_mov_opt_anc_it_mean, OU_mov_opt_it_mean, Stasis_Stasis_it_mean, Stasis_URW_it_mean, Stasis_GRW_it_mean, 
                    Stasis_OU_it_mean, URW_URW_it_mean, URW_GRW_it_mean, URW_OU_it_mean, GRW_GRW_it_mean, GRW_OU_it_mean, 
                    OU_OU_it_mean),
  
  standard_deviation_passed = c(GRW_it_sd, URW_it_sd, Stasis_it_sd, Strict_stasis_it_sd, Decel_it_sd, Accel_it_sd, OU_it_sd, 
                         OU_mov_opt_anc_it_sd, OU_mov_opt_it_sd, Stasis_Stasis_it_sd, Stasis_URW_it_sd, Stasis_GRW_it_sd, 
                         Stasis_OU_it_sd, URW_URW_it_sd, URW_GRW_it_sd, URW_OU_it_sd, GRW_GRW_it_sd, GRW_OU_it_sd, 
                         OU_OU_it_sd),
  
  percentage_size_passed = c(GRW_size_p, URW_size_p, stasis_size_p, strict_stasis_size_p, decel_size_p, accel_size_p, OU_size_p, 
                      OU_mov_opt_anc_size_p, OU_mov_opt_size_p, Stasis_Stasis_size_p, Stasis_URW_size_p, Stasis_GRW_size_p, 
                      Stasis_OU_size_p, URW_URW_size_p, URW_GRW_size_p, URW_OU_size_p, GRW_GRW_size_p, GRW_OU_size_p, 
                      OU_OU_size_p),
  
  percentage_shape_passed = c(GRW_shape_p, URW_shape_p, stasis_shape_p, strict_stasis_shape_p, decel_shape_p, accel_shape_p, OU_shape_p, 
                             OU_mov_opt_anc_shape_p, OU_mov_opt_shape_p, Stasis_Stasis_shape_p, Stasis_URW_shape_p, Stasis_GRW_shape_p, 
                             Stasis_OU_shape_p, URW_URW_shape_p, URW_GRW_shape_p, URW_OU_shape_p, GRW_GRW_shape_p, GRW_OU_shape_p, 
                             OU_OU_shape_p)
  
)

Total_adeq_passed = sum(GRW_c, URW_c, stasis_c, strict_stasis_c, decel_c, accel_c, OU_c, 
                        OU_mov_opt_anc_c, OU_mov_opt_c, Stasis_Stasis_c, Stasis_URW_c, Stasis_GRW_c, 
                        Stasis_OU_c, URW_URW_c, URW_GRW_c, URW_OU_c, GRW_GRW_c, GRW_OU_c, 
                        OU_OU_c, OU_GRW_c, OU_URW_c, OU_Stasis_c, GRW_URW_c, GRW_Stasis_c, 
                        URW_Stasis_c)

# write to file
sink(file = "./results_paleoTS_v0.6.1/table_result_shift.txt")
resultshift_table
paste("Total number of time series investigated:", length(model_noshift_results))
paste("Total number of time series which passed adequacy tests:", Total_adeq_passed)
paste("Percentage of time series which passed adequacy tests:", (Total_adeq_passed*100)/length(model_noshift_results))
sink()

#----------------------------------------------------------------
# Plot of the number of time series best fitted for each models
#----------------------------------------------------------------

metadatalong_clear$best_model <- aicc_min

metadatalong_clear <- metadatalong_clear %>%
  mutate(best_model = case_when(
    aicc_min == 1 ~ 'GRW',
    aicc_min == 2 ~ 'URW',
    aicc_min == 3 ~ 'Stasis',
    aicc_min == 4 ~ 'Strict_stasis',
    aicc_min == 5 ~ 'Decel',
    aicc_min == 6 ~ 'Accel',
    aicc_min == 7 ~ 'OU',
    aicc_min == 8 ~ 'OU_mov_opt_anc',
    aicc_min == 9 ~ 'OU_mov_opt',
    aicc_min == 10 ~ 'Stasis_Stasis',
    aicc_min == 11 ~ 'Stasis_URW',
    aicc_min == 12 ~ 'Stasis_GRW',
    aicc_min == 13 ~ 'Stasis_OU',
    aicc_min == 14 ~ 'URW_URW',
    aicc_min == 15 ~ 'URW_GRW',
    aicc_min == 16 ~ 'URW_OU',
    aicc_min == 17 ~ 'GRW_GRW',
    aicc_min == 18 ~ 'GRW_OU',
    aicc_min == 19 ~ 'OU_OU',
    aicc_min == 20 ~ 'OU_GRW',
    aicc_min == 21 ~ 'OU_URW',
    aicc_min == 22 ~ 'OU_Stasis',
    aicc_min == 23 ~ 'GRW_URW',
    aicc_min == 24 ~ 'GRW_Stasis',
    aicc_min == 25 ~ 'URW_Stasis',
    TRUE ~ as.character(best_model)
  ))

# Adequate timeseries

bar_complete <- ggplot(metadatalong_adequate, aes(x = best_model)) +
  stat_count(geom = "bar", fill = "blue") +
  labs(title = "Barplot of the best model among adequate time series", x = "Best evolutionary model", y = "Number of time series")


# Inadequate timeseries

#remove inadequate time series before performing the analysis
metadatalong_adequate <- metadatalong_clear %>%
  filter(adequacy_status != "inadequate")



bar_adequate <- ggplot(metadatalong_adequate, aes(x = best_model)) +
  stat_count(geom = "bar", fill = "blue") +
  labs(title = "Barplot of the best model among adequate time series", x = "Best evolutionary model", y = "Number of time series")


# save the graphs
png("./results_paleoTS_v0.6.1/plot/results_count_bestfit.png", width = 2400, height = 1800)
grid.arrange(bar_complete, bar_adequate, ncol = 1)
dev.off()





