# Code for "Finite sample performance of optimal treatment rule estimators with right-censored outcomes."

## Author: Michael Jetsupphasuk

The materials presented here can be used to reproduce the results in the manuscript. See below for a description of each file.

- `estimators.R`: compares value estimators
- `eval.R`: computes and aggregates results after estimation of optimal treatment rules
- `get_treat_cv_sc*.R`: estimates optimal treatment rules with sample splitting estimators for specified scenario
- `get_treat_whole_sc.R*`: estimates optimal treatment rules without sample splitting estimators for specified scenario
- `sims_create.R`: creates simulation datasets
- `sims_summaries.R`: computes summary statistics for simulation datasets
- `utils_sim.R`: functions used in other files

Download the R script `RISTfunctions.r` from https://sites.google.com/site/teazrq/software. 
