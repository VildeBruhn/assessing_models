## packages are loaded in the main script "Assessing models of evolution - 
## UNIVARIATE MODELS"

#------------------------
# Data manipulation (dt)
#------------------------

dt <- function(df, tsID){
  
  # split into two lists; one with log transformed and one without
  df_log <- list()
  df_not_log <- list()
  
  df_log <- lapply(df, function(x){
    if(isTRUE(x$log_transformed[1] == "yes") == TRUE){
      df_log$x = x
    }
  })
  
  df_log = df_log[-which(sapply(df_log, is.null))]
  
  df_not_log <- lapply(df, function(x){
    if(isTRUE(x$log_transformed[1] == "no") == TRUE){
      df_not_log$x = x
    }
  })
  
  df_not_log = df_not_log[-which(sapply(df_not_log, is.null))]
  
  
  ## CREATE paleoTS OBJECTS AND TRANSFORM DATA ##
  
  # apply the as.paleoTS function to all data frames in the list
  data_log <- lapply(df_log, function(x) {
    as.paleoTS(mm = x$trait_mean, vv = x$trait_var, nn = x$N, tt = x$age_MY, oldest = "first")
  })
  
  data_not_log <- lapply(df_not_log, function(x) {
    as.paleoTS(mm = x$trait_mean, vv = x$trait_var, nn = x$N, tt = x$age_MY, oldest = "first")
  })
  
  # approximate log-transformation of the data sets that do not have log transformed values
  data_not_log_log <- lapply(data_not_log, ln.paleoTS)
  #nan <- as.data.frame(which(rapply(data_not_log_log, is.nan)))
  
  # join log and not log data sets again now that all are on log scale
  data_log <- mapply(c, data_log, df_log, SIMPLIFY = FALSE)
  data_not_log_log <- mapply(c, data_not_log_log, df_not_log, SIMPLIFY = FALSE)
  ln_data <- c(data_log, data_not_log_log)
  
  return(ln_data)
}

#---------------------------------------------------
# Get timeseries that adequate
#---------------------------------------------------

adequate3tests <- function(data){
  
  # get only passed or failed adequacy test into list
  data_adeq <- lapply(data, function(x) {
    x$result = toString(x$summary[5])
    return(x)
  })
  
  # filter out those that did not pass the adequacy test
  data_adeq_passed <- Filter(function(x) x$result ==  "c(\"PASSED\", \"PASSED\", \"PASSED\")", data_adeq)
  
  return(data_adeq_passed)
  
}

adequate4tests <- function(data){
  
  # get only passed or failed adequacy test into list
  data_adeq <- lapply(data, function(x) {
    x$results = toString(x$summary[5])
    return(x)
  })
  
  # filter out those that did not pass the adequacy test
  data_adeq_passed <- Filter(function(x) x$results ==  "c(\"PASSED\", \"PASSED\", \"PASSED\", \"PASSED\")", data_adeq)
  
  return(data_adeq_passed)
  
}

# this one does not include the last adequacy test in adePEM::fit3adequacy
adequate2tests <- function(data){
  
  # get only passed or failed adequacy test into list
  data_adeq <- lapply(data, function(x) {
    x$results = toString(x$summary[5])
    return(x)
  })
  
  # filter out those that did not pass the adequacy test
  data_adeq_passed1 <- Filter(function(x) x$results ==  "c(\"PASSED\", \"PASSED\", \"PASSED\")", data_adeq)
  data_adeq_passed2 <- Filter(function(x) x$results ==  "c(\"PASSED\", \"PASSED\", \"FAILED\")", data_adeq)
  data_adeq_passed3 <- Filter(function(x) x$results ==  "c(\"PASSED\", \"PASSED\", \"NA\")", data_adeq)
  data_adeq_passed <- c(data_adeq_passed1, data_adeq_passed2, data_adeq_passed3)
  
  return(data_adeq_passed)
  
}

### Add model info to metadata ###

model_aicc <- function(data, aicc_results, model){
  data <- lapply(data, function(x){
    if (any(names(aicc_results) == x$tsID[1]) == TRUE){
      x$model_aicc <- model
    } else {
      data <- data
    }
    return(x)
  })
}

model_type <- function(data, aicc_results, model){
  data <- lapply(data, function(x){
    if (any(names(aicc_results) == x$tsID[1]) == TRUE){
      x$model_type <- model
    } else {
      data <- data
    }
    return(x)
  })
}

model_adeq <- function(data, adeq_results, model){
  data <- lapply(data, function(x){
    if (any(names(adeq_results) == x$tsID[1]) == TRUE){
      x$model_adequate <- model
    } else {
      data <- data 
    }
    return(x)
  })
}

### binding list to dataframe ###
bind <- function(data, unit_list){
  binded <- lapply(data, function(x) x[(names(x) %in% unit_list)])
  binded <- bind_rows(binded, .id="data_frame")
  # remove duplicated data
  binded <- binded[!duplicated(binded),]
}
