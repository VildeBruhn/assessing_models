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
library(devtools)
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
doParallel::registerDoParallel(cl = n_cores)



#-----------------
# IMPORT FILES
#-----------------

# import
timeseries <- read_delim("./timeseries/timeseries.txt", col_names = TRUE, delim = "\t")
metadata <- read_delim("./timeseries/metadata.txt", col_names = TRUE, delim = "\t")



##################
## RELATIVE FIT ##
##################

#---------------------------------------------------------------------------------------------
# EXCLUDING THE SHORTEST TIME SERIES (Keep time series only if containing 14 or more steps)
#---------------------------------------------------------------------------------------------

# join dataframes
df <- left_join(timeseries, metadata, by = c("tsID"))

# remove time series with less than 14 steps
df_shift <- subset(df, steps >= 14)

# remove modern time series
df_shift <- subset(df_shift, period_start != "Present")

# make list based on ID
df_shift <- lapply(split(df_shift, df_shift$tsID), function(x) as.list(x))

# process data
ln_data_meta_shift <- dt(df_shift, "tsID")

ln_data_shift <- lapply(ln_data_meta_shift, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$N, tt = x$tt, oldest = "first")
})

# Convert the time vector to unit length
ln_data_meta_shift <- lapply(ln_data_meta_shift, function(x) {
  x$tt <- x$tt/(max(x$tt))
  x
})

ln_data_shift <- lapply(ln_data_shift, function(x) {
  x$tt <- x$tt/(max(x$tt))
  x
})


#-------------------------------------------------
# FIT UNIVARIATE MODELS WITHOUT AND WITH A SHIFT
#-------------------------------------------------
# for paleoTS v0.6.2

###### MODEL FIT WAS RUN ON THE HPC IN TWO SUBSETS #######
# Create two subsets to make analysis faster
ln_data_shift_subset1 = ln_data_shift[0:(length(ln_data_shift)/2)] 
ln_data_shift_subset2 = ln_data_shift[(length(ln_data_shift)/2 + 1):length(ln_data_shift)]

# SUBSET 1
# test all possible univariate models with no shift from evoTS on every timeseries
model_noshift_results_subset1 <- foreach(i = seq_along(ln_data_shift_subset1),
                                         .packages = c("evoTS"),
                                         .combine = "c") %dopar% {
                                           nm <- names(ln_data_shift_subset1)[i]
                                           res <- tryCatch(
                                             fit.all.univariate(ln_data_shift_subset1[[i]]),
                                             error = function(e) {
                                               message(paste("Model", nm, "failed:", e$message))
                                               NULL
                                             }
                                           )
                                           setNames(list(res), nm)
                                         }

# function to iterate through each combination of model with a shift
fit_mode_shift <- function(ln_data_shift_subset1) {
  models_list <- c("Stasis", "URW", "GRW", "OU")
  store_results <- list()
  k <- 0
  for (i in 1:4) {
    model1 <- models_list[i]
    for (j in 1:4) {
      model2 <- models_list[j]
      # Safe fitting with error handling
      fit_result <- tryCatch(
        fit.mode.shift(ln_data_shift_subset1, model1, model2, minb = 7),
        error = function(e) {
          message(paste("Shift model", model1, "->", model2, "failed:", e$message))
          return(NULL)
        }
      )
      k <- k + 1
      store_results[[k]] <- fit_result
    }
  }
  return(store_results)
}

# test all possible univariate models with shift from evoTS on every timeseries
model_shift_results_subset1 <- foreach(i = seq_along(ln_data_shift_subset1),
                                       .packages = c("evoTS"),
                                       .combine = "c") %dopar% {
                                         nm <- names(ln_data_shift_subset1)[i]
                                         res <- tryCatch(
                                           fit_mode_shift(ln_data_shift_subset1[[i]]),
                                           error = function(e) {
                                             message(paste("Model", nm, "failed:", e$message))
                                             NULL
                                           }
                                         )
                                         setNames(list(res), nm)
                                       }

# Save the results and the data for subset 1
save(model_noshift_results_subset1, model_shift_results_subset1, file = "model_test_shift_subset1.RData")
save(ln_data_shift_subset1, file = "ln_data_shift_subset1.RData")


# SUBSET 2
# test all possible univariate models with no shift from evoTS on every timeseries
model_noshift_results_subset2 <- foreach(i = seq_along(ln_data_shift_subset2),
                                         .packages = c("evoTS"),
                                         .combine = "c") %dopar% {
                                           nm <- names(ln_data_shift_subset2)[i]
                                           if (is.null(nm) || nzchar(nm) == FALSE) nm <- as.character(i)
                                           
                                           res <- tryCatch(
                                             fit.all.univariate(ln_data_shift_subset2[[i]]),
                                             error = function(e) {
                                               message(paste("Model", nm, "failed:", e$message))
                                               NULL
                                             }
                                           )
                                           setNames(list(res), nm)
                                         }


fit_mode_shift <- function(ln_data_shift_subset2) {
  models_list <- c("Stasis", "URW", "GRW", "OU")
  store_results <- list()
  k <- 0
  for (i in 1:4) {
    model1 <- models_list[i]
    for (j in 1:4) {
      model2 <- models_list[j]
      # Safe fitting with error handling
      fit_result <- tryCatch(
        fit.mode.shift(ln_data_shift_subset2, model1, model2, minb = 7),
        error = function(e) {
          message(paste("Shift model", model1, "->", model2, "failed:", e$message))
          return(NULL)
        }
      )
      k <- k + 1
      store_results[[k]] <- fit_result
    }
  }
  return(store_results)
}

# test all possible univariate models with shift from evoTS on every timeseries
model_shift_results_subset2 <- foreach(i = seq_along(ln_data_shift_subset2),
                                       .packages = c("evoTS"),
                                       .combine = "c") %dopar% {
                                         nm <- names(ln_data_shift_subset2)[i]
                                         res <- tryCatch(
                                           fit_mode_shift(ln_data_shift_subset2[[i]]),
                                           error = function(e) {
                                             message(paste("Model", nm, "failed:", e$message))
                                             NULL
                                           }
                                         )
                                         setNames(list(res), nm)
                                       }


# Save the results and the data for subset 2
save(model_noshift_results_subset2, model_shift_results_subset2, file = "model_test_shift_subset2.RData")
save(ln_data_shift_subset2, file = "ln_data_shift_subset2.RData")

# Loading datasets and results from the HPC
load("ln_data_shift_subset1.RData") 
load("ln_data_shift_subset2.RData")
load("model_test_shift_subset1.RData") 
load("model_test_shift_subset2.RData") 

# Merging the result subsets
model_noshift_results = c(model_noshift_results_subset1, model_noshift_results_subset2)
model_shift_results = c(model_shift_results_subset1, model_shift_results_subset2)


#----------------------------
# ASSESS RELATIVE FIT (AICc)
#----------------------------

### Remove problematic timeseries ###
model_noshift_results <- Filter(function(iteration) !is.null(iteration), model_noshift_results)
model_shift_results <- Filter(\(iteration) length(iteration) == 16, model_shift_results)

model_noshift_results_TS = names(model_noshift_results)
model_shift_results <- model_shift_results[names(model_shift_results) %in% model_noshift_results_TS]

model_shift_results_TS = names(model_shift_results)
model_noshift_results <- model_noshift_results[names(model_noshift_results) %in% model_shift_results_TS]

ln_data_shift = ln_data_shift[names(ln_data_shift) %in% model_shift_results_TS]
ln_data_meta_shift = ln_data_meta_shift[names(ln_data_meta_shift) %in% model_shift_results_TS]

# Save the dataset used for the analysis with mode shift
save(ln_data_shift, file = "ln_data_shift.RData")
save(ln_data_meta_shift, file = "ln_data_meta_shift.RData")

# Extract AICcs for the models without shift
aicc_noshift <- lapply(model_noshift_results, function(x) x[(names(x) %in% c("AICc"))]) #USE model_noshift_results_clean if those time series are still problematic

# Extract AICc values of shift models on all results
aicc_shift_extraction <- lapply(model_shift_results, function(TS) { 
  sapply(TS, function(result) result$AICc)
})
modelshift_names <- sapply(model_shift_results[[1]], function(model) model$modelName)
aicc_shift_extraction <- lapply(aicc_shift_extraction, function(TS) {
  names(TS) <- modelshift_names
  return(TS)
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

names(aicc_results_complete) <- rownames(aicc[[1]])

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
paste("Total number of time series described by more than one model (threshold = 2):", length(aicc_filtered))
paste("Percentage of time series filtered (threshold = 2):", length(aicc_filtered)/length(aicc)*100)
sink()


  
##################
## ABSOLUTE FIT ##
##################

# Add a column with lowest AIC for each time series
for (i in 1:length(ln_data_shift)) {
  ln_data_shift[[i]]$Lowest_AICc <- aicc_min[[i]]
}


#----------------------------------------------------
# Filter the time series according to the best model
#----------------------------------------------------

categories <- c("GRW", "URW", "Stasis", "Strict_stasis", "Decel", "Accel", "OU",
                "OU_mov_opt_anc", "OU_mov_opt", "Stasis_Stasis", 
                "Stasis_URW", "Stasis_GRW", "Stasis_OU", "URW_Stasis", "URW_URW", "URW_GRW", "URW_OU",
                "GRW_Stasis","GRW_URW","GRW_GRW", "GRW_OU","OU_Stasis","OU_URW","OU_GRW", "OU_OU")

# Create a list to store the results
result_list <- list()

for (i in 1:length(categories)) {
  category <- categories[i]
  
  # Filter data for the current category
  filtered_data <- Filter(function(x) x[[10]] == i, ln_data_shift)
  filtered_data <- lapply(filtered_data, function(x) { x[[10]] <- NULL; x })
  filtered_data <- lapply(filtered_data, function(x) {
    as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt)
  })
  
  # Store the result in the result_list
  result_list[[category]] <- filtered_data
  assign(paste(category, sep = ""), filtered_data)
}

# Save the data
save(GRW, URW, Stasis,
     Strict_stasis, Decel, Accel, 
     OU, OU_mov_opt, OU_mov_opt_anc, 
     Stasis_Stasis, Stasis_URW, Stasis_GRW, Stasis_OU, 
     URW_URW, URW_GRW, URW_OU,
     GRW_GRW, GRW_OU, 
     OU_OU, OU_GRW, OU_URW, OU_Stasis, 
     GRW_URW, GRW_Stasis, URW_Stasis,
     file = "./aicc_shift_passed.RData")
                          
#-----------------------------------------------------------------------
# Splitting the time series best described by models with a shift models
#-----------------------------------------------------------------------
                          
### Stasis-Stasis ###
# Add a column with the best shift point to each timeseries
tsID_Stasis_Stasis = names(Stasis_Stasis)
model_results_Stasis_Stasis <- model_shift_results[tsID_Stasis_Stasis]

for (i in tsID_Stasis_Stasis) {
  Stasis_Stasis[[i]]$Shift_point <- model_shift_results[[i]][[1]]$parameters[["shift1"]]  #1 is for the model Stasis_Stasis, need to be changed for the other models
}

#Splitting the model
Stasis_Stasis_subset1 <- setNames(vector("list", length(Stasis_Stasis)), names(Stasis_Stasis))
Stasis_Stasis_subset2 <- setNames(vector("list", length(Stasis_Stasis)), names(Stasis_Stasis))

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
Stasis_URW_subset1 = setNames(vector("list", length(Stasis_URW)), names(Stasis_URW))
Stasis_URW_subset2 = setNames(vector("list", length(Stasis_URW)), names(Stasis_URW))

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
Stasis_GRW_subset1 = setNames(vector("list", length(Stasis_GRW)), names(Stasis_GRW))
Stasis_GRW_subset2 = setNames(vector("list", length(Stasis_GRW)), names(Stasis_GRW))

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
Stasis_OU_subset1 = setNames(vector("list", length(Stasis_OU)), names(Stasis_OU))
Stasis_OU_subset2 = setNames(vector("list", length(Stasis_OU)), names(Stasis_OU))

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
URW_URW_subset1 = setNames(vector("list", length(URW_URW)), names(URW_URW))
URW_URW_subset2 = setNames(vector("list", length(URW_URW)), names(URW_URW))

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
URW_GRW_subset1 = setNames(vector("list", length(URW_GRW)), names(URW_GRW))
URW_GRW_subset2 = setNames(vector("list", length(URW_GRW)), names(URW_GRW))

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
URW_OU_subset1 = setNames(vector("list", length(URW_OU)), names(URW_OU))
URW_OU_subset2 = setNames(vector("list", length(URW_OU)), names(URW_OU))

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
GRW_GRW_subset1 = setNames(vector("list", length(GRW_GRW)), names(GRW_GRW))
GRW_GRW_subset2 = setNames(vector("list", length(GRW_GRW)), names(GRW_GRW))

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
GRW_OU_subset1 = setNames(vector("list", length(GRW_OU)), names(GRW_OU))
GRW_OU_subset2 = setNames(vector("list", length(GRW_OU)), names(GRW_OU))

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
OU_OU_subset1 = setNames(vector("list", length(OU_OU)), names(OU_OU))
OU_OU_subset2 = setNames(vector("list", length(OU_OU)), names(OU_OU))

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
OU_GRW_subset1 = setNames(vector("list", length(OU_GRW)), names(OU_GRW))
OU_GRW_subset2 = setNames(vector("list", length(OU_GRW)), names(OU_GRW))

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
OU_URW_subset1 = setNames(vector("list", length(OU_URW)), names(OU_URW))
OU_URW_subset2 = setNames(vector("list", length(OU_URW)), names(OU_URW))

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
OU_Stasis_subset1 = setNames(vector("list", length(OU_Stasis)), names(OU_Stasis))
OU_Stasis_subset2 = setNames(vector("list", length(OU_Stasis)), names(OU_Stasis))

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
GRW_URW_subset1 = setNames(vector("list", length(GRW_URW)), names(GRW_URW))
GRW_URW_subset2 = setNames(vector("list", length(GRW_URW)), names(GRW_URW))

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
GRW_Stasis_subset1 = setNames(vector("list", length(GRW_Stasis)), names(GRW_Stasis))
GRW_Stasis_subset2 = setNames(vector("list", length(GRW_Stasis)), names(GRW_Stasis))

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
URW_Stasis_subset1 = setNames(vector("list", length(URW_Stasis)), names(URW_Stasis))
URW_Stasis_subset2 = setNames(vector("list", length(URW_Stasis)), names(URW_Stasis))

if (length(URW_Stasis) > 0) {
  for (i in 1:length(URW_Stasis)) {
    gg <- rep(1:2, c(URW_Stasis[[i]]$Shift_point, length(URW_Stasis[[i]]$mm) - URW_Stasis[[i]]$Shift_point))
    URW_Stasis_split = paleoTS:::split4punc(URW_Stasis[[i]],gg, overlap=TRUE)
    URW_Stasis_subset1[[i]] = URW_Stasis_split[[1]]
    URW_Stasis_subset2[[i]] = URW_Stasis_split[[2]]
  }
}

# Saving OU subsets to work on them in another script (the adequacy tests of OU models is not working in parallel so need to be implemented in a loop)
# see assessing_models_shift_OU_adequacy.R for the code
save(OU, OU_mov_opt_anc, OU_mov_opt, Stasis_OU_subset2, URW_OU_subset2, GRW_OU_subset2,
     OU_OU_subset1, OU_OU_subset2, OU_GRW_subset1, OU_URW_subset1, OU_Stasis_subset1,
     file = "./OU_shift.RData")

#------------------------------------
# Testing the adequacy of the models
#------------------------------------

# test adequacy
GRW_adeq <- mclapply(GRW, fit3adequacy.trend, plot = FALSE)
URW_adeq <- mclapply(URW, fit3adequacy.RW, plot = FALSE)
stasis_adeq <- mclapply(Stasis, fit4adequacy.stasis, plot = FALSE) 

strict_stasis_adeq <- mclapply(
  Strict_stasis,
  function(x) {
    tryCatch(
      fit4adequacy.stasis(x, plot = FALSE),
      error = function(e) NA
    )
  }
)
adeq_issues_stasis <- which(sapply(strict_stasis_adeq, function(x) is.na(x)[1]))
strict_stasis_adeq <- Filter(function(x) !is.na(x)[1], strict_stasis_adeq)

decel_adeq <- mclapply(Decel, fit3adequacy.decel, plot = FALSE)

# reverse accelerated to become decelerated
Accel_Decel <- Accel
Accel_Decel <- lapply(Accel_Decel, function(x) {
  x$mm <- rev(x$mm)
  x$vv <- rev(x$vv)
  x$nn <- rev(x$nn)
  x$tt <- rev(x$tt)
  for (i in 1:length(x$tt)){
    x$tt[i] <- 1 - x$tt[i]
  }
  return(x)
})
accel_adeq <- mclapply(Accel_Decel, fit3adequacy.decel, plot = FALSE)

#OU_adeq <- mclapply(OU, fit3adequacy.OU, plot = FALSE)
#OU_mov_opt_anc_adeq <- mclapply(OU_mov_opt_anc, fit3adequacy.OU, plot = FALSE)
#OU_mov_opt_adeq <- mclapply(OU_mov_opt, fit3adequacy.OU, plot = FALSE)

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
#names(Stasis_OU_subset2_adeq) = tsID_Stasis_OU

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
#names(URW_OU_subset2_adeq) = tsID_URW_OU

GRW_GRW_subset1_adeq <- mclapply(GRW_GRW_subset1, fit3adequacy.trend, plot = FALSE)
GRW_GRW_subset2_adeq <- mclapply(GRW_GRW_subset2, fit3adequacy.trend, plot = FALSE)
names(GRW_GRW_subset1_adeq) = tsID_GRW_GRW
names(GRW_GRW_subset2_adeq) = tsID_GRW_GRW

GRW_OU_subset1_adeq <- mclapply(GRW_OU_subset1, fit3adequacy.trend, plot = FALSE)
#GRW_OU_subset2_adeq <- mclapply(GRW_OU_subset2, fit3adequacy.OU, plot = FALSE)
names(GRW_OU_subset1_adeq) = tsID_GRW_OU
#names(GRW_OU_subset2_adeq) = tsID_GRW_OU

#OU_OU_subset1_adeq <- mclapply(OU_OU_subset1, fit3adequacy.OU, plot = FALSE)
#OU_OU_subset2_adeq <- mclapply(OU_OU_subset2, fit3adequacy.OU, plot = FALSE)
#names(OU_OU_subset1_adeq) = tsID_OU_OU 
#names(OU_OU_subset2_adeq) = tsID_OU_OU

#OU_GRW_subset1_adeq <- mclapply(OU_GRW_subset1, fit3adequacy.OU, plot = FALSE)
OU_GRW_subset2_adeq <- mclapply(OU_GRW_subset2, fit3adequacy.trend, plot = FALSE)
#names(OU_GRW_subset1_adeq) = tsID_OU_GRW
names(OU_GRW_subset2_adeq) = tsID_OU_GRW

#OU_URW_subset1_adeq <- mclapply(OU_URW_subset1, fit3adequacy.OU, plot = FALSE)
OU_URW_subset2_adeq <- mclapply(OU_URW_subset2, fit3adequacy.RW, plot = FALSE)
#names(OU_URW_subset1_adeq) = tsID_OU_URW
names(OU_URW_subset2_adeq) = tsID_OU_URW

#OU_Stasis_subset1_adeq <- mclapply(OU_Stasis_subset1, fit3adequacy.OU, plot = FALSE)
#names(OU_Stasis_subset1_adeq) = tsID_OU_Stasis 
OU_Stasis_subset2_adeq <- mclapply(
  OU_Stasis_subset2,
  function(x) {
    tryCatch(
      fit4adequacy.stasis(x, plot = FALSE),
      error = function(e) NA
    )
  }
)
names(OU_Stasis_subset2_adeq) = tsID_OU_Stasis
adeq_issues_OU_Stasis_subset2 <- which(sapply(OU_Stasis_subset2_adeq, function(x) is.na(x)[1]))
OU_Stasis_subset2_adeq <- Filter(function(x) !is.na(x)[1], OU_Stasis_subset2_adeq)

GRW_URW_subset1_adeq <- mclapply(GRW_URW_subset1, fit3adequacy.trend, plot = FALSE)
names(GRW_URW_subset1_adeq) = tsID_GRW_URW
GRW_URW_subset2_adeq <- mclapply(
  GRW_URW_subset2,
  function(x) {
    tryCatch(
      fit3adequacy.RW(x, plot = FALSE),
      error = function(e) NA
    )
  }
)
names(GRW_URW_subset2_adeq) = tsID_GRW_URW
adeq_issues_GRW_URW_subset2 <- which(sapply(GRW_URW_subset2_adeq, function(x) is.na(x)[1]))
GRW_URW_subset2_adeq <- Filter(function(x) !is.na(x)[1], GRW_URW_subset2_adeq)

GRW_Stasis_subset1_adeq <- mclapply(GRW_Stasis_subset1, fit3adequacy.trend, plot = FALSE)
GRW_Stasis_subset2_adeq <- mclapply(GRW_Stasis_subset2, fit4adequacy.stasis, plot = FALSE)
names(GRW_Stasis_subset1_adeq) = tsID_GRW_Stasis
names(GRW_Stasis_subset2_adeq) = tsID_GRW_Stasis

URW_Stasis_subset1_adeq <- mclapply(URW_Stasis_subset1, fit3adequacy.RW, plot = FALSE)
URW_Stasis_subset2_adeq <- mclapply(URW_Stasis_subset2, fit4adequacy.stasis, plot = FALSE)
names(URW_Stasis_subset1_adeq) = tsID_URW_Stasis
names(URW_Stasis_subset2_adeq) = tsID_URW_Stasis

# Loading results of the OU adequacy (the OU models are not working in parallel)
load("./OU_shift_adeq.RData")
adeq_issues = c(adeq_issues_stasis, adeq_issues_OU_Stasis_subset2, adeq_issues_GRW_URW_subset2, adeq_issues_OU)

# get adequacy results for only adequate time series
GRW_adeq_passed <- adequate3tests(GRW_adeq)
URW_adeq_passed <- adequate3tests(URW_adeq)
Stasis_adeq_passed <- adequate4tests(stasis_adeq)
Strict_stasis_adeq_passed <- adequate4tests(strict_stasis_adeq)
Decel_adeq_passed <- adequate3tests(decel_adeq)
Accel_adeq_passed <- adequate3tests(accel_adeq)
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
stasis_list_adequate <- names(Stasis_adeq_passed)
stasis_adeq_list <- Stasis[names(Stasis) %in% stasis_list_adequate]
strict_stasis_list_adequate <- names(Strict_stasis_adeq_passed)
strict_stasis_adeq_list <- Strict_stasis[names(Strict_stasis) %in% strict_stasis_list_adequate]
decel_list_adequate <- names(Decel_adeq_passed)
decel_adeq_list <- Decel[names(Decel) %in% decel_list_adequate]
accel_list_adequate <- names(Accel_adeq_passed)
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

# Save the time series which passed adequacy tests
save(GRW_adeq_passed, URW_adeq_passed, Stasis_adeq_passed,
     Strict_stasis_adeq_passed, Decel_adeq_passed, Accel_adeq_passed, 
     OU_adeq_passed, OU_mov_opt_adeq_passed, OU_mov_opt_anc_adeq_passed, 
     Stasis_Stasis_adeq_passed, Stasis_URW_adeq_passed, Stasis_GRW_adeq_passed, Stasis_OU_adeq_passed, 
     URW_URW_adeq_passed, URW_GRW_adeq_passed, URW_OU_adeq_passed,
     GRW_GRW_adeq_passed, GRW_OU_adeq_passed, 
     OU_OU_adeq_passed, OU_GRW_adeq_passed, OU_URW_adeq_passed, OU_Stasis_adeq_passed, 
     GRW_URW_adeq_passed, GRW_Stasis_adeq_passed, URW_Stasis_adeq_passed,
file = "adeq_shift_passed.Rdata")


#---------------------------------------------------------------------------------------------------------------------------------
# Table S2. Detail results of the models including mode-shift explaining time series best according to relative and absolute fit.
#---------------------------------------------------------------------------------------------------------------------------------

# get the number of parameters for each model
K_noshift <- model_noshift_results[[1]]$K
names(K_noshift) <- rownames(model_noshift_results[[1]])

K_shift <- sapply(model_shift_results[[1]], `[[`, "K")
names(K_shift) <- sapply(model_shift_results[[1]], `[[`, "modelName")

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
Stasis_URW_rc = length(Stasis_URW) 
Stasis_GRW_rc = length(Stasis_GRW) 
Stasis_OU_rc = length(Stasis_OU) 
URW_URW_rc = length(URW_URW)
URW_GRW_rc = length(URW_GRW) 
URW_OU_rc = length(URW_OU) 
GRW_GRW_rc = length(GRW_GRW)
GRW_OU_rc = length(GRW_OU) 
OU_OU_rc = length(OU_OU)
URW_Stasis_rc = length(URW_Stasis)
GRW_Stasis_rc = length(GRW_Stasis)
OU_Stasis_rc = length(OU_Stasis)
GRW_URW_rc = length(GRW_URW)
OU_URW_rc = length(OU_URW)
OU_GRW_rc = length(OU_GRW)

total_rc = sum(GRW_rc, URW_rc, stasis_rc, strict_stasis_rc, decel_rc, accel_rc, OU_rc, 
              OU_mov_opt_anc_rc, OU_mov_opt_rc, Stasis_Stasis_rc, Stasis_URW_rc, Stasis_GRW_rc, 
              Stasis_OU_rc, URW_URW_rc, URW_GRW_rc, URW_OU_rc, GRW_GRW_rc, GRW_OU_rc, 
              OU_OU_rc, OU_GRW_rc, OU_URW_rc, OU_Stasis_rc, GRW_URW_rc, GRW_Stasis_rc, 
              URW_Stasis_rc)

total_notrad_rc = total_rc - (GRW_rc + URW_rc + stasis_rc + strict_stasis_rc)
total_notrad_rp = total_notrad_rc/total_rc*100

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
URW_Stasis_rp = URW_Stasis_rc/length(model_noshift_results)*100
GRW_Stasis_rp = GRW_Stasis_rc/length(model_noshift_results)*100
OU_Stasis_rp = OU_Stasis_rc/length(model_noshift_results)*100
GRW_URW_rp = GRW_URW_rc/length(model_noshift_results)*100
OU_URW_rp = OU_URW_rc/length(model_noshift_results)*100
OU_GRW_rp = OU_GRW_rc/length(model_noshift_results)*100

# Absolute count - get counts of time series which passed adequacy tests
GRW_c <- length(GRW_adeq_passed)
URW_c <- length(URW_adeq_passed)
stasis_c <- length(Stasis_adeq_passed)
strict_stasis_c <- length(Strict_stasis_adeq_passed)
decel_c <- length(Decel_adeq_passed)
accel_c <- length(Accel_adeq_passed)
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

total_c = sum(GRW_c, URW_c, stasis_c, strict_stasis_c, decel_c, accel_c, OU_c, 
                        OU_mov_opt_anc_c, OU_mov_opt_c, Stasis_Stasis_c, Stasis_URW_c, Stasis_GRW_c, 
                        Stasis_OU_c, URW_URW_c, URW_GRW_c, URW_OU_c, GRW_GRW_c, GRW_OU_c, 
                        OU_OU_c, OU_GRW_c, OU_URW_c, OU_Stasis_c, GRW_URW_c, GRW_Stasis_c, 
                        URW_Stasis_c)

total_notrad_c = total_c - (GRW_c + URW_c + stasis_c + strict_stasis_c)
total_notrad_p = total_notrad_c/total_c*100

# Adequacy success rate - get counts of time series which passed adequacy tests on the total of time series for the model investigated
GRW_p <- (length(GRW_adeq_passed)/length(GRW_adeq))*100
URW_p <- (length(URW_adeq_passed)/length(URW_adeq))*100
stasis_p <- (length(Stasis_adeq_passed)/length(stasis_adeq))*100
strict_stasis_p <- (length(Strict_stasis_adeq_passed)/length(strict_stasis_adeq))*100
decel_p <- (length(Decel_adeq_passed)/length(decel_adeq))*100
accel_p <- (length(Accel_adeq_passed)/length(accel_adeq))*100
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

# Time series which could not be evaluated in the adequacy tests
TS_noneval_adeq = length(adeq_issues)

# make output table
adeq_table <- data.frame(
  
  parameters = c(K_noshift, K_shift),
  
  relative_count = c(GRW_rc, URW_rc, stasis_rc, strict_stasis_rc, decel_rc, accel_rc, OU_rc, 
                      OU_mov_opt_anc_rc, OU_mov_opt_rc, Stasis_Stasis_rc, Stasis_URW_rc, Stasis_GRW_rc, 
                      Stasis_OU_rc, URW_Stasis_rc, URW_URW_rc, URW_GRW_rc, URW_OU_rc, GRW_Stasis_rc, 
                      GRW_URW_rc, GRW_GRW_rc, GRW_OU_rc, OU_Stasis_rc, OU_URW_rc, OU_GRW_rc, OU_OU_rc),
  
  relative_percentage = c(GRW_rp, URW_rp, stasis_rp, strict_stasis_rp, decel_rp, accel_rp, OU_rp, 
                          OU_mov_opt_anc_rp, OU_mov_opt_rp, Stasis_Stasis_rp, Stasis_URW_rp, Stasis_GRW_rp, 
                          Stasis_OU_rp, URW_Stasis_rp, URW_URW_rp, URW_GRW_rp, URW_OU_rp, GRW_Stasis_rp, 
                          GRW_URW_rp, GRW_GRW_rp, GRW_OU_rp, OU_Stasis_rp, OU_URW_rp, OU_GRW_rp, OU_OU_rp),
  
  count_passed = c(GRW_c, URW_c, stasis_c, strict_stasis_c, decel_c, accel_c, OU_c, 
                   OU_mov_opt_anc_c, OU_mov_opt_c, Stasis_Stasis_c, Stasis_URW_c, Stasis_GRW_c, 
                   Stasis_OU_c, URW_Stasis_c, URW_URW_c, URW_GRW_c, URW_OU_c, GRW_Stasis_c, 
                   GRW_URW_c, GRW_GRW_c, GRW_OU_c, OU_Stasis_c, OU_URW_c, OU_GRW_c, OU_OU_c),
  
  percentage_passed = c(GRW_p, URW_p, stasis_p, strict_stasis_p, decel_p, accel_p, OU_p, 
                        OU_mov_opt_anc_p, OU_mov_opt_p, Stasis_Stasis_p, Stasis_URW_p, Stasis_GRW_p, 
                        Stasis_OU_p, URW_Stasis_p, URW_URW_p, URW_GRW_p, URW_OU_p, GRW_Stasis_p, 
                        GRW_URW_p, GRW_GRW_p, GRW_OU_p, OU_Stasis_p, OU_URW_p, OU_GRW_p, OU_OU_p))


# write to file
sink(file = "./results_paleoTS_v0.6.1/Results_fits_with_shift.txt")
adeq_table
paste("Total count    ", total_rc, "    ", total_c)
paste("Time series not explained by stasis, unbiased random walk, or general random walk    ", round(total_notrad_rp,2), "    ", round(total_notrad_p,2))
paste("Time series from which absolute fit could not be assessed    ", TS_noneval_adeq)
sink()



#------------------------------------------------------------------------------------------------------------------------------------
# Table 3. Results of the types of models including mode-shift explaining time series best according to relative and absolute fit.
#------------------------------------------------------------------------------------------------------------------------------------

singlemode_trad_rc = sum(GRW_rc, URW_rc, stasis_rc, strict_stasis_rc)
singlemode_else_rc = sum(decel_rc, accel_rc, OU_rc, OU_mov_opt_anc_rc, OU_mov_opt_rc)
modeshift_rc = sum(Stasis_Stasis_rc, Stasis_URW_rc, Stasis_GRW_rc, Stasis_OU_rc, URW_URW_rc, 
                   URW_GRW_rc, URW_OU_rc, GRW_GRW_rc, GRW_OU_rc, OU_OU_rc, OU_GRW_rc, 
                   OU_URW_rc, OU_Stasis_rc, GRW_URW_rc, GRW_Stasis_rc, URW_Stasis_rc)
  
singlemode_trad_rp = sum(GRW_rp, URW_rp, stasis_rp, strict_stasis_rp)
singlemode_else_rp = sum(decel_rp, accel_rp, OU_rp, OU_mov_opt_anc_rp, OU_mov_opt_rp)
modeshift_rp = sum(Stasis_Stasis_rp, Stasis_URW_rp, Stasis_GRW_rp, Stasis_OU_rp, URW_URW_rp, 
                   URW_GRW_rp, URW_OU_rp, GRW_GRW_rp, GRW_OU_rp, OU_OU_rp, OU_GRW_rp, 
                   OU_URW_rp, OU_Stasis_rp, GRW_URW_rp, GRW_Stasis_rp, URW_Stasis_rp)
  
singlemode_trad_c = sum(GRW_c, URW_c, stasis_c, strict_stasis_c)
singlemode_else_c = sum(decel_c, accel_c, OU_c, OU_mov_opt_anc_c, OU_mov_opt_c)
modeshift_c = sum(Stasis_Stasis_c, Stasis_URW_c, Stasis_GRW_c, Stasis_OU_c, URW_URW_c, 
                   URW_GRW_c, URW_OU_c, GRW_GRW_c, GRW_OU_c, OU_OU_c, OU_GRW_c, 
                   OU_URW_c, OU_Stasis_c, GRW_URW_c, GRW_Stasis_c, URW_Stasis_c)

singlemode_trad_p = singlemode_trad_c/singlemode_trad_rc*100
singlemode_else_p = singlemode_else_c/singlemode_else_rc*100
modeshift_p = modeshift_c/modeshift_rc*100

# make output table
adeq_table_summary <- data.frame(
  
  parameters = c("2-3", "3-5", "4-8"),
  
  relative_count = c(singlemode_trad_rc, singlemode_else_rc, modeshift_rc),
  
  relative_percentage = c(singlemode_trad_rp, singlemode_else_rp, modeshift_rp),
  
  count_passed = c(singlemode_trad_c, singlemode_else_c, modeshift_c),
  
  percentage_passed = c(singlemode_trad_p, singlemode_else_p, modeshift_p))

# write to file
sink(file = "./results_paleoTS_v0.6.1/Results_fits_with_shift_summary.txt")
adeq_table_summary
paste("Total count    ", total_rc, "    ", total_c)
paste("Time series from which absolute fit could not be assessed    ", TS_noneval_adeq)
sink()