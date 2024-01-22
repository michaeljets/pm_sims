# Computes survival estimates from treatment rules

source('code/utils_sim.R')

library(survival)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# set ggplot theme
theme_set(theme_bw())

# read in results for a particular sim
scenario = 0

# sample splitting or no sample splitting
cv_or_whole = 'cv'

# read in data
data_string = paste0('data/simlist_sc', scenario, '.rds')
datalist = readRDS(data_string)

# place to save results
save_string = paste0('results/sim_res_', cv_or_whole, '_sc', scenario, '.rds')


# Get results -------------------------------------------------------------


# initialize matrix to store results
num_sims = 500
est_0 = est_1 = est_d =
  true_0 = true_1 = true_d =
  matrix(NA, nrow=num_sims, ncol=60)
  
est_dhat = true_dhat = misclass = 
  est_utility = true_utility =
  matrix(NA, nrow=num_sims, ncol=12)

est_prac_utility = true_prac_utility =
  numeric(num_sims)

# loop over sims
for (sim in 1:num_sims) {
  
  # filename prefix
  prefix = ifelse(cv_or_whole=='cv', 'results/cv/sc', 'results/whole/whole_sc')
  
  # read in dtr_objs
  if (scenario!=3) {
    filename = paste0(prefix, scenario, '_sim', sim, '.rds')
  } else {
    filename = paste0(prefix, 5, '_sim', sim, '.rds')
    
  }
  dtr_objs = readRDS(filename)
  if (cv_or_whole=='cv') dtr_objs = dtr_objs[-grep('whole', names(dtr_objs))]
  names(dtr_objs) = gsub('owl_lin', 'owl_rist', names(dtr_objs))
  names(dtr_objs) = gsub('owl_rad', 'owl_rist_rad', names(dtr_objs))
  
  # accidentally saved bart twice for some
  if (sum('bart'==names(dtr_objs)) == 2) (dtr_objs[[12]] = NULL)
  
  # read in data
  data = datalist[[sim]]
  
  # get estimated survival
  est_surv_all = get_est_curvs(data, dtr_objs, get_true_dtr(scenario), cv_or_whole)
  
  est_0[sim,] = est_surv_all$est_0
  est_1[sim,] = est_surv_all$est_1
  est_d[sim,] = est_surv_all$est_d
  est_dhat[sim,] = sapply(est_surv_all$est_dhat, function(x) x[60])
  
  # get true survival
  true_surv_all = get_true_curvs(data, scenario, est_surv_all$dhat)
  
  true_0[sim,] = true_surv_all$true_0
  true_1[sim,] = true_surv_all$true_1
  true_d[sim,] = true_surv_all$true_d
  true_dhat[sim,] = sapply(true_surv_all$true_dhat, function(x) x[60])
  misclass[sim,] = true_surv_all$misclass
  
  # calculate utilities
  est_utility[sim,] = est_dhat[sim,] - true_d[sim, 60]
  true_utility[sim,] = true_dhat[sim,] - true_d[sim, 60]
  est_prac_utility[sim] = max(est_dhat[sim,]) - true_d[sim, 60]
  true_prac_utility[sim] = max(true_dhat[sim,]) - true_d[sim, 60]
  
  # track progress
  if (sim %% 10 == 0) {
    print(sim)
    print(Sys.time())
  }
  
}


# get method names
method_name = names(dtr_objs)

# save results
saveRDS(mget(c('est_0', 'est_1', 'est_d',
               'true_0', 'true_1', 'true_d',
               'est_dhat', 'true_dhat', 'misclass',
               'est_utility', 'true_utility',
               'est_prac_utility', 'true_prac_utility',
               'method_name')),
        save_string)


# Summaries ---------------------------------------------------------------


# basic checks
boxplot(est_0-true_0, ylim=c(-0.02,0.02)); abline(h=0,col='red')
boxplot(est_1-true_1, ylim=c(-0.02,0.02)); abline(h=0,col='red')
boxplot(est_d-true_d, ylim=c(-0.02,0.02)); abline(h=0,col='red')


# descriptive method names
desc_names = c('ridge' = 'Ridge Cox regression',
               'lasso' = 'Lasso Cox regression',
               'elastic' = 'Elastic net Cox regression',
               'csf' = 'Causal survival forests',
               'bart' = 'Bayesian additive regression trees',
               'owl_rist' = 'OWL, linear kernel',
               'owl_rist_rad' = 'OWL, radial kernel',
               'rwl' = 'Residual weighted learning',
               'genetic' = 'Genetic',
               'genetic_a' = 'Genetic, augment',
               'genetic_s' = 'Genetic, smooth',
               'genetic_as' = 'Genetic, augment + smooth',
               'Maximum' = 'Maximum')


# creating one plot with facets
est_utility_lst = true_utility_lst = misclass_lst = list()
for (i in 0:3) {
  
  # read in results
  read_string = paste0('results/sim_res_', cv_or_whole, '_sc', i, '.rds')
  sim_results = readRDS(read_string)
  est_utility = sim_results$est_utility
  true_utility = sim_results$true_utility
  est_prac_utility = sim_results$est_prac_utility
  true_prac_utility = sim_results$true_prac_utility
  misclass = sim_results$misclass
  method_name = sim_results$method_name
  
  # store
  est_utility_lst[[i+1]] = format_results(est_utility, method_name, desc_names, est_prac_utility)
  true_utility_lst[[i+1]] = format_results(true_utility, method_name, desc_names, true_prac_utility)
  misclass_lst[[i+1]] = format_results(misclass, method_name, desc_names)
  
}


# plot results
gg_est = plot_results(combine_results(est_utility_lst), zero_line=TRUE)
gg_true = plot_results(combine_results(true_utility_lst))
gg_mis = plot_results(combine_results(misclass_lst))

# save plots
ggsave(paste0('results/est_utility_', cv_or_whole, '.tiff'), 
       gg_est,
       width = 7, height = 6,
       units = 'in')

ggsave(paste0('results/true_utility_', cv_or_whole, '.tiff'), 
       gg_true,
       width = 7, height = 6,
       units = 'in')

ggsave(paste0('results/misclass_', cv_or_whole, '.tiff'), 
       gg_mis,
       width = 7, height = 6,
       units = 'in')



