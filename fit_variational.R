library(tidyverse)
source("src/fit_BDML.R")
source("src/fit_BLR.R")
source("src/plot_results.R")

# load data

murder_data <- read.csv("data/murder_data-2.csv")
violence_data <- read.csv("data/violence_data-2.csv")
property_data <- read.csv("data/property_data-2.csv")

### fit models on murder data

# STAN BDML MODELS
murder_fit_b <- fit_model_dml_b(murder_data, variational = TRUE) |> extract_results_BDML("BDML-Basic") # takes ~18h on my machine if not using variational approximation
murder_fit_b2 <- fit_model_dml_b2(murder_data, variational = TRUE) |> extract_results_BDML("BDML-Hier")
# BLR MODELS
murder_fit_BRL <- fit_model_BLR(murder_data)
murder_fit_BRL_naive <- murder_fit_BRL$naive |> extract_results_BLR("Naive")
murder_fit_BRL_hahn <- murder_fit_BRL$hahn |> extract_results_BLR("Hahn")
murder_fit_BRL_linero <- murder_fit_BRL$linero |> extract_results_BLR("Linero")
# Frequentist DML models
murder_fit_FDML <- fit_model_FDML(murder_data)
murder_fit_FDML_split <- murder_fit_FDML$FDML_split |> extract_results_FDML("FDML-Split")
murder_fit_FDML_full <- murder_fit_FDML$FDML_full |> extract_results_FDML("FDML-Full")

# append rows to table
murder_result_rows <- list(murder_fit_b, murder_fit_b2, murder_fit_BRL_naive, murder_fit_BRL_hahn, murder_fit_BRL_linero, murder_fit_FDML_split, murder_fit_FDML_full)
murder_results_table <- do.call(rbind, murder_result_rows)
# save to csvff
write.csv(murder_results_table, "results/murder_results_table.csv")

### fit models on violence data

# STAN BDML MODELS
violence_fit_b <- fit_model_dml_b(violence_data, variational = TRUE) |> extract_results_BDML("BDML-Basic")
violence_fit_b2 <- fit_model_dml_b2(violence_data, variational = TRUE) |> extract_results_BDML("BDML-Hier")
# BLR MODELS
violence_fit_BRL <- fit_model_BLR(violence_data)
violence_fit_BRL_naive <- violence_fit_BRL$naive |> extract_results_BLR("Naive")
violence_fit_BRL_hahn <- violence_fit_BRL$hahn |> extract_results_BLR("Hahn")
violence_fit_BRL_linero <- violence_fit_BRL$linero |> extract_results_BLR("Linero")
# Frequentist DML models
violence_fit_FDML <- fit_model_FDML(violence_data)
violence_fit_FDML_split <- violence_fit_FDML$FDML_split |> extract_results_FDML("FDML-Split")
violence_fit_FDML_full <- violence_fit_FDML$FDML_full |> extract_results_FDML("FDML-Full")

# append rows to table
violence_result_rows <- list(violence_fit_b, violence_fit_b2, violence_fit_BRL_naive, violence_fit_BRL_hahn, violence_fit_BRL_linero, violence_fit_FDML_split, violence_fit_FDML_full)
violence_results_table <- do.call(rbind, violence_result_rows)
# save to csv
write.csv(violence_results_table, "results/violence_results_table.csv")


### fit models on property data

# STAN BDML MODELS
property_fit_b <- fit_model_dml_b(property_data, variational = TRUE) |> extract_results_BDML("BDML-Basic")
property_fit_b2 <- fit_model_dml_b2(property_data, variational = TRUE) |> extract_results_BDML("BDML-Hier")
# BLR MODELS
property_fit_BRL <- fit_model_BLR(property_data)
property_fit_BRL_naive <- property_fit_BRL$naive |> extract_results_BLR("Naive")
property_fit_BRL_hahn <- property_fit_BRL$hahn |> extract_results_BLR("Hahn")
property_fit_BRL_linero <- property_fit_BRL$linero |> extract_results_BLR("Linero")
# Frequentist DML models
property_fit_FDML <- fit_model_FDML(property_data)
property_fit_FDML_split <- property_fit_FDML$FDML_split |> extract_results_FDML("FDML-Split")
property_fit_FDML_full <- property_fit_FDML$FDML_full |> extract_results_FDML("FDML-Full")

# append rows to table
property_result_rows <- list(property_fit_b, property_fit_b2, property_fit_BRL_naive, property_fit_BRL_hahn, property_fit_BRL_linero, property_fit_FDML_split, property_fit_FDML_full)
property_results_table <- do.call(rbind, property_result_rows)
# save to csv
write.csv(property_results_table, "results/property_results_table.csv")

# plot and save results
plot_gamma_estimates(property_results_table, "property_gamma_estimates_variational")
plot_gamma_estimates(murder_results_table, "property_gamma_estimates_variational")
plot_gamma_estimates(violence_results_table, "property_gamma_estimates_variational")


