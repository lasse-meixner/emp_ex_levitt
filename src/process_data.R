library(tidyverse)

## Function to create data list ----
create_data_list <- function(raw_data) {
  # first extract Y, and D as outcomes
    Y <- raw_data[, 1] # outcome
    D <- raw_data[, 2] # treatment
    X <- as.matrix(raw_data[, -c(1, 2)]) # large set of controls
    list(K = 2, J = ncol(X), N = nrow(X), X = X, Y = Y, D = D)
}