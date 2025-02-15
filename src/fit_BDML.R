## Load necessary libraries ----
library(cmdstanr)
library(tidyverse)
source("src/process_data.R")

## Compile the Stan models ----
model_dml_b <- cmdstan_model("dml_b.stan")
model_dml_b2 <- cmdstan_model("dml_b2.stan")
model_dml_r2d2 <- cmdstan_model("dml_r2d2.stan")

## Function to fit model B ----

fit_model_dml_b <- function(raw_data, variational=TRUE) {
    
    data_list <- create_data_list(raw_data)
    
    if (variational) {
        fitted_dml_b <- model_dml_b$variational(
            data = list(K = 2, J = data_list$J, N = data_list$N, x = data_list$X, y = cbind(data_list$Y, data_list$D)),
            iter = 5e5,
            output_samples = 5000
        )
    } else {
        fitted_dml_b <- model_dml_b$sample(
            data = list(K = 2, J = data_list$J, N = data_list$N, x = data_list$X, y = cbind(data_list$Y, data_list$D)),
            chains = 4,
            parallel_chains = 4,
            refresh = 0,
            show_messages = TRUE,
            show_exceptions = TRUE,
        )
    }
}

## Function to fit model B2 (hierarchical) ----

fit_model_dml_b2 <- function(raw_data, variational=TRUE) {

    data_list <- create_data_list(raw_data)

    if (variational) {
        fitted_dml_b2 <- model_dml_b2$variational(
            data = list(K = 2, J = data_list$J, N = data_list$N, x = data_list$X, y = cbind(data_list$Y, data_list$D)),
            iter = 5e5,
            output_samples = 5000
        )
    } else {
        fitted_dml_b2 <- model_dml_b2$sample(
            data = list(K = 2, J = data_list$J, N = data_list$N, x = data_list$X, y = cbind(data_list$Y, data_list$D)),
            chains = 4,
            parallel_chains = 4,
            refresh = 0,
            show_messages = TRUE,
            show_exceptions = TRUE,
        )
    }
}

## Function to fit model R2D2 ----

fit_model_dml_r2d2 <- function(raw_data) {

    data_list <- create_data_list(raw_data)
  
    fitted_dml_r2d2 <- model_dml_r2d2$sample(
        data = list(J = data_list$J, N = data_list$N, x = data_list$X, y = cbind(data_list$Y, data_list$D), b = 0.5),
        chains = 4,
        parallel_chains = 4,
        refresh = 0,
        show_messages = TRUE,
        show_exceptions = TRUE,
    )
}


extract_results_BDML <- function(fit, name) {
    draws <- fit$draws("alpha")
    fit_stats <- data.frame(
        mean = mean(draws),
        q2.5 = quantile(draws, 0.025),
        q97.5 = quantile(draws, 0.975)
    )
    gamma_hat <- fit_stats$mean
    interval <- c(fit_stats$q2.5, fit_stats$q97.5)
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