## Application of the Bayesian Double Machine Learning Method to the empirical example from Belloni, Chernozhukov, and Hansen (2014)

### Structure

The original code from BCH(2014) can be found in `BCH-2014-replication/`. 

The pre-processed data from their code is saved in `data/`. A snippet of the results table from their double-variable-selection demonstration on this data from their paper is saved as `BCH2014_table2.png`.

The main entry point for fitting our models to this data is `fit.R`. This depends on the data files in `data/`, the STAN model files saved in the root folder, and the auxiliary source scripts in `src/`. Results of our model fits are saved in `results/`.

### Old

The files in `R_replication/` are deprecate files not used in the final analysis.

