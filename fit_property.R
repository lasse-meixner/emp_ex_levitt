library(tidyverse)
source("src/fit_BDML.R")
source("src/fit_BLR.R")
source("src/plot_results.R")

# load data
property_data <- read.csv("data/property_data-2.csv")

### fit models on property data

# STAN BDML MODELS
property_fit_b <- fit_model_dml_b(property_data, variational = FALSE) |> 
  extract_results_BDML("BDML-Basic")
property_fit_b2 <- fit_model_dml_b2(property_data, variational = FALSE) |> 
  extract_results_BDML("BDML-Hier")

# BLR MODELS
property_fit_BRL <- fit_model_BLR(property_data)
property_fit_BRL_naive <- property_fit_BRL$naive |> extract_results_BLR("Naive")
property_fit_BRL_hahn <- property_fit_BRL$hahn  |> extract_results_BLR("Hahn")
property_fit_BRL_linero <- property_fit_BRL$linero |> extract_results_BLR("Linero")

# Frequentist DML models
property_fit_FDML <- fit_model_FDML(property_data)
property_fit_FDML_split <- property_fit_FDML$FDML_split |> extract_results_FDML("FDML-Split")
property_fit_FDML_full  <- property_fit_FDML$FDML_full  |> extract_results_FDML("FDML-Full")

# Combine rows into a results table
property_result_rows <- list(
  property_fit_b, 
  property_fit_b2, 
  property_fit_BRL_naive, 
  property_fit_BRL_hahn, 
  property_fit_BRL_linero, 
  property_fit_FDML_split, 
  property_fit_FDML_full
)
property_results_table <- do.call(rbind, property_result_rows)

# Save to CSV
write.csv(property_results_table, "results/property_results_table.csv")

# Plot and save results
plot_gamma_estimates(property_results_table, "property_gamma_estimates")
