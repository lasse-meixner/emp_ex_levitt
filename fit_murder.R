library(tidyverse)
source("src/fit_BDML.R")
source("src/fit_BLR.R")
source("src/plot_results.R")

# load data
murder_data <- read.csv("data/murder_data-2.csv")

### fit models on murder data

# STAN BDML MODELS
murder_fit_b <- fit_model_dml_b(murder_data, variational = FALSE) |> 
  extract_results_BDML("BDML-Basic") 
murder_fit_b2 <- fit_model_dml_b2(murder_data, variational = FALSE) |> 
  extract_results_BDML("BDML-Hier")

# BLR MODELS
murder_fit_BRL <- fit_model_BLR(murder_data)
murder_fit_BRL_naive <- murder_fit_BRL$naive |> extract_results_BLR("Naive")
murder_fit_BRL_hahn <- murder_fit_BRL$hahn  |> extract_results_BLR("Hahn")
murder_fit_BRL_linero <- murder_fit_BRL$linero |> extract_results_BLR("Linero")

# Frequentist DML models
murder_fit_FDML <- fit_model_FDML(murder_data)
murder_fit_FDML_split <- murder_fit_FDML$FDML_split |> extract_results_FDML("FDML-Split")
murder_fit_FDML_full  <- murder_fit_FDML$FDML_full  |> extract_results_FDML("FDML-Full")

# Combine rows into a results table
murder_result_rows <- list(
  murder_fit_b, 
  murder_fit_b2, 
  murder_fit_BRL_naive, 
  murder_fit_BRL_hahn, 
  murder_fit_BRL_linero, 
  murder_fit_FDML_split, 
  murder_fit_FDML_full
)
murder_results_table <- do.call(rbind, murder_result_rows)

# Save to CSV
write.csv(murder_results_table, "results/murder_results_table.csv")

# Plot and save results
plot_gamma_estimates(murder_results_table, "murder_gamma_estimates")
