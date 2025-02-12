library(tidyverse)
library(BLR)
source("src/process_data.R")


## Function to fit BLRs ----

fit_model_BLR <- function(data) {

  data_list <- create_data_list(data)
  
  invisible(capture.output({ # Suppress output from BLR calls
    # 1. Hahn & Linero
    naive <- BLR(y = data_list$Y, XR = data_list$X, XF = as.matrix(data_list$D))
    fitted_ps <- BLR(y = data_list$D, XR = data_list$X)
    hahn <- BLR(y = data_list$Y, XF = cbind(data_list$D - fitted_ps$yHat), XR = data_list$X)
    linero <- BLR(y = data_list$Y, XF = cbind(data_list$D, fitted_ps$yHat), XR = data_list$X)
  })) # end of silenced BLR calls

  # Ensure the return object is not printed
  invisible(list(naive = naive, hahn = hahn, linero = linero))
}


fit_model_FDML <- function(data) {

  data_list <- create_data_list(data)

  # 2. FDML
  fitted_ps <- BLR(y = data_list$D, XR = data_list$X)
  fitted_a_step1 <- fitted_ps
  fitted_y_step1 <- BLR(y = data_list$Y, XR = data_list$X)
     
  a_res_step2 <- data_list$D - fitted_a_step1$mu - data_list$X %*% fitted_a_step1$bR
  y_res_step2 <- data_list$Y - fitted_y_step1$mu - data_list$X %*% fitted_y_step1$bR

  fitted_dml_full <- lm(y_res_step2 ~ a_res_step2)

  n <- dim(data_list$X)[1] # number of observations
  ix_step1 <- 1:(floor(n/2))
  fitted_a_step1 <- BLR(y = data_list$D[ix_step1], XR = data_list$X[ix_step1,])
  fitted_y_step1 <- BLR(y = data_list$Y[ix_step1], XR = data_list$X[ix_step1,])

  ix_step2 <- (floor(n/2) + 1):n
  a_res_step2 <- data_list$D[ix_step2] - fitted_a_step1$mu - data_list$X[ix_step2,] %*% fitted_a_step1$bR
  y_res_step2 <- data_list$Y[ix_step2] - fitted_y_step1$mu - data_list$X[ix_step2,] %*% fitted_y_step1$bR

  fitted_dml_split <- lm(y_res_step2 ~ a_res_step2)

  return(list(FDML_full = fitted_dml_full, FDML_split = fitted_dml_split))

}

extract_results_BLR <- function(fit, name) {
  gamma_hat <- fit$bF[1]
  interval <- gamma_hat + c(-1.96,1.96) * fit$SD.bF[1] # cf. get_interval in Linero/src/get_interval.R
  # Todo: might be interesting to plot the posteriors of gamma_hat
  interval_width <- interval[2] - interval[1]
  LCL <- interval[1]
  UCL <- interval[2]

  table <- data.frame(
    method = name,
    gamma_hat = gamma_hat, 
    LCL = LCL,
    UCL = UCL,
    interval_width = interval_width
  )
}

extract_results_FDML <- function(fit, name) {
  gamma_hat <- summary(fit)$coefficients[2,1]
  interval <- unname(confint(fit)[2,])
  interval_width <- interval[2] - interval[1]
  LCL <- interval[1]
  UCL <- interval[2]

  table <- data.frame(
    method = name,
    gamma_hat = gamma_hat, 
    LCL = LCL,
    UCL = UCL,
    interval_width = interval_width
  )
}