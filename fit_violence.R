library(tidyverse)
source("src/fit_BDML.R")
source("src/fit_BLR.R")
source("src/plot_results.R")

# load data
violence_data <- read.csv("data/violence_data-2.csv")

### fit models on violence data

# STAN BDML MODELS
violence_fit_b <- fit_model_dml_b(violence_data, variational = FALSE) |> 
  extract_results_BDML("BDML-Basic")
violence_fit_b2 <- fit_model_dml_b2(violence_data, variational = FALSE) |> 
  extract_results_BDML("BDML-Hier")

# BLR MODELS
violence_fit_BRL <- fit_model_BLR(violence_data)
violence_fit_BRL_naive <- violence_fit_BRL$naive |> extract_results_BLR("Naive")
violence_fit_BRL_hahn <- violence_fit_BRL$hahn  |> extract_results_BLR("Hahn")
violence_fit_BRL_linero <- violence_fit_BRL$linero |> extract_results_BLR("Linero")

# Frequentist DML models
violence_fit_FDML <- fit_model_FDML(violence_data)
violence_fit_FDML_split <- violence_fit_FDML$FDML_split |> extract_results_FDML("FDML-Split")
violence_fit_FDML_full  <- violence_fit_FDML$FDML_full  |> extract_results_FDML("FDML-Full")

# Combine rows into a results table
violence_result_rows <- list(
  violence_fit_b, 
  violence_fit_b2, 
  violence_fit_BRL_naive, 
  violence_fit_BRL_hahn, 
  violence_fit_BRL_linero, 
  violence_fit_FDML_split, 
  violence_fit_FDML_full
)
violence_results_table <- do.call(rbind, violence_result_rows)

# Save to CSV
write.csv(violence_results_table, "results/violence_results_table.csv")

# Plot and save results
plot_gamma_estimates(violence_results_table, "violence_gamma_estimates")
