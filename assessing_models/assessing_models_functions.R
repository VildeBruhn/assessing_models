library(evoTS)
library(adePEM)
library(tidyverse)
library(data.table)

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
  
  # check if oldest is first, change those with youngest first to oldest first
  df_log <- lapply(df_log, function(x) {
    if (isTRUE(x$oldest_first[1] == "no") == TRUE){
      arrange(x, desc(age_MY))
    } else if (isTRUE(x$oldest_first[1] == "yes") == TRUE){
      x <- x
    }
  })
  
  df_not_log <- lapply(df_not_log, function(x) {
    if (isTRUE(x$oldest_first[1] == "no") == TRUE){
      arrange(x, desc(age_MY))
    } else if (isTRUE(x$oldest_first[1] == "yes") == TRUE){
      x <- x
    }
  })
  
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

#-------------------------------------------------------------------
# fit4adequacy.stasis (doesn't work on strict stasis through adePEM)
#-------------------------------------------------------------------

auto_corr_test_stasis <- function (y, nrep = 1000, conf = 0.95, plot = TRUE, save.replicates = TRUE, 
                                   omega = NULL) 
{
  x <- y$mm
  v <- y$vv
  n <- y$nn
  tt <- y$tt
  lower <- (1 - conf)/2
  upper <- (1 + conf)/2
  theta <- opt.joint.Stasis(y)$parameters[1]
  if (is.null(omega)) 
    omega <- opt.joint.Stasis(y)$parameters[2]
  obs.auto.corr <- auto.corr(x, model = "stasis")
  bootstrap.matrix <- matrix(data = NA, nrow = nrep, ncol = 1)
  for (i in 1:nrep) {
    x_sim <- sim.Stasis(ns = length(x), theta = theta, omega = omega, 
                        vp = v, nn = n, tt = tt)
    bootstrap.matrix[i, 1] <- auto.corr(x_sim$mm, model = "stasis")
  }
  bootstrap.auto.corr <- length(bootstrap.matrix[, 1][bootstrap.matrix[, 
                                                                       1] > obs.auto.corr])/nrep
  if (bootstrap.auto.corr > round(upper, 3) | bootstrap.auto.corr < 
      round(lower, 3)) 
    pass.auto.corr.test <- "FAILED"
  else pass.auto.corr.test <- "PASSED"
  if (bootstrap.auto.corr > 0.5) 
    bootstrap.auto.corr <- 1 - bootstrap.auto.corr
  if (plot == TRUE) {
    layout(1:1)
    plotting.distributions(bootstrap.matrix[, 1], obs.auto.corr, 
                           test = "auto.corr", xlab = "Simulated data", main = "Autocorrelation")
  }
  output <- as.data.frame(cbind(round(obs.auto.corr, 5), round(min(bootstrap.matrix), 
                                                               5), round(max(bootstrap.matrix), 5), bootstrap.auto.corr/0.5, 
                                pass.auto.corr.test), nrow = 5, byrow = TRUE)
  rownames(output) <- "auto.corr"
  colnames(output) <- c("estimate", "min.sim", "max.sim", "'p-value'", 
                        "result")
  summary.out <- as.data.frame(c(nrep, conf))
  rownames(summary.out) <- c("replications", "confidence level")
  colnames(summary.out) <- ("Value")
  if (save.replicates == FALSE) {
    out <- list(info = summary.out, summary = output)
    return(out)
  }
  else {
    out <- list(replicates = bootstrap.matrix, info = summary.out, 
                summary = output)
    return(out)
  }
}

sim.Stasis <- function (ns = 20, theta = 0, omega = 0, vp = 1, nn = rep(20, 
                                                                        ns), tt = 0:(ns - 1)) 
{
  xmu <- rnorm(ns, mean = theta, sd = sqrt(omega))
  xobs <- xmu + rnorm(ns, 0, sqrt(vp/nn))
  gp <- c(theta, omega)
  names(gp) <- c("theta", "omega")
  x <- as.paleoTS(mm = xobs, vv = rep(vp, 1), nn = nn, tt = tt, 
                  MM = xmu, genpars = gp, label = "Created by sim.Stasis", 
                  reset.time = FALSE)
  return(x)
}

adeq_stasis <- function (y, nrep = 1000, conf = 0.95, plot = FALSE, omega = NULL) 
{
  x <- y$mm
  v <- y$vv
  n <- y$nn
  tt <- y$tt
  theta <- opt.joint.Stasis(y)$parameters[1]
  if (is.null(omega)) 
    omega <- opt.joint.Stasis(y)$parameters[2]
  lower <- (1 - conf)/2
  upper <- (1 + conf)/2
  obs.auto.corr <- auto.corr(x, model = "stasis")
  obs.runs.test <- runs.test(x, model = "stasis", theta = theta)
  obs.slope.test <- slope.test(x, tt, model = "stasis", theta = theta)
  obs.net.change.test <- net.change.test(x, model = "stasis")
  out.auto <- auto_corr_test_stasis(y, nrep, conf, plot = FALSE, 
                                    theta, omega)
  out.runs <- runs.test.stasis(y, nrep, conf, plot = FALSE, 
                               theta, omega)
  out.slope <- slope.test.stasis(y, nrep, conf, plot = FALSE, 
                                 theta, omega)
  out.net <- net.change.test.stasis(y, nrep, conf, plot = FALSE, 
                                    theta, omega)
  output <- c(as.vector(matrix(unlist(out.auto[[3]]), ncol = 5, 
                               byrow = FALSE)), as.vector(matrix(unlist(out.runs[[3]]), 
                                                                 ncol = 5, byrow = FALSE)), as.vector(matrix(unlist(out.slope[[3]]), 
                                                                                                             ncol = 5, byrow = FALSE)), as.vector(matrix(unlist(out.net[[3]]), 
                                                                                                                                                         ncol = 5, byrow = FALSE)))
  output <- as.data.frame(cbind(c(output[c(1, 6, 11, 16)]), 
                                c(output[c(2, 7, 12, 17)]), c(output[c(3, 8, 13, 18)]), 
                                c(output[c(4, 9, 14, 19)]), c(output[c(5, 10, 15, 20)])), 
                          ncol = 5)
  rownames(output) <- c("auto.corr", "runs.test", "slope.test", 
                        "net.change.test")
  colnames(output) <- c("estimate", "min.sim", "max.sim", "p-value", 
                        "result")
  if (plot == TRUE) {
    par(mfrow = c(2, 2))
    model.names <- c("auto.corr", "runs.test", "slope.test", 
                     "net.change.test")
    plotting.distributions(out.auto$replicates, obs.auto.corr, 
                           model.names[1], xlab = "Simulated data", main = "Autocorrelation")
    plotting.distributions(out.runs$replicates, obs.runs.test, 
                           model.names[2], xlab = "Simulated data", main = "Runs")
    plotting.distributions(out.slope$replicates, obs.slope.test, 
                           model.names[3], xlab = "Simulated data", main = "Fixed variance")
    plotting.distributions(out.net$replicates, obs.net.change.test, 
                           model.names[4], xlab = "Simulated data", main = "Net evolution")
  }
  summary.out <- as.data.frame(c(nrep, conf))
  rownames(summary.out) <- c("replications", "confidence level")
  colnames(summary.out) <- ("Value")
  out <- list(info = summary.out, summary = output)
  return(out)
}

#------------------------------------------------
# fit3adequacy.OU (doesn't work through adePEM)
#------------------------------------------------

adeq_OU <- function (y, nrep = 1000, conf = 0.95, cutoff = 0.8, plot = FALSE) 
{
  x <- y$mm
  v <- y$vv
  n <- y$nn
  tt <- y$tt
  anc <- opt.joint.OU(y)$parameters[1]
  vstep <- opt.joint.OU(y)$parameters[2]
  theta <- opt.joint.OU(y)$parameters[3]
  alpha <- opt.joint.OU(y)$parameters[4]
  tmp_OU <- opt.joint.OU(y)
  pred_OU <- est.OU(y, tmp_OU, tt = tt)
  detrended_OU <- x - pred_OU$ee
  lower <- (1 - conf)/2
  upper <- (1 + conf)/2
  obs.auto.corr <- auto.corr(detrended_OU, model = "OU", tt)
  obs.runs.test <- runs.test(detrended_OU, model = "OU", tt)
  obs_sum_of_residuals <- 0
  out.auto <- auto.corr.test.OU(y, nrep, conf, plot = FALSE, 
                                save.replicates = TRUE)
  out.runs <- runs.test.OU(y, nrep, conf, plot = FALSE, save.replicates = TRUE)
  out.var <- variance.test.OU(y, nrep, cutoff, plot = FALSE, 
                              save.replicates = TRUE)
  output <- c(as.vector(matrix(unlist(out.auto[[3]]), ncol = 5, 
                               byrow = FALSE)), as.vector(matrix(unlist(out.runs[[3]]), 
                                                                 ncol = 5, byrow = FALSE)), as.vector(matrix(unlist(out.var[[3]]), 
                                                                                                             ncol = 5, byrow = FALSE)))
  output <- as.data.frame(cbind(c(output[c(1, 6, 11)]), c(output[c(2, 
                                                                   7, 12)]), c(output[c(3, 8, 13)]), c(output[c(4, 9, 14)]), 
                                c(output[c(5, 10, 15)])), ncol = 5)
  rownames(output) <- c("auto.corr", "runs.test", "slope.test")
  colnames(output) <- c("estimate", "min.sim", "max.sim", "p-value", 
                        "result")
  if (plot == TRUE) {
    par(mfrow = c(1, 3))
    model.names <- c("auto.corr", "runs.test", "slope.test")
    plotting.distributions(out.auto$replicates, obs.auto.corr, 
                           model.names[1], xlab = "Simulated data", main = "Autocorrelation")
    plotting.distributions(out.runs$replicates, obs.runs.test, 
                           model.names[2], xlab = "Simulated data", main = "Runs")
    plotting.distributions(out.var$replicates, obs_sum_of_residuals, 
                           model.names[3], xlab = "Simulated data", main = "Initial rapid change")
  }
  summary.out <- as.data.frame(c(nrep, conf, cutoff))
  rownames(summary.out) <- c("replications", "confidence level", 
                             "cut-off faster evolution")
  colnames(summary.out) <- ("Value")
  out <- list(info = summary.out, summary = output)
  return(out)
}
