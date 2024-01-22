# creates / documents simulated datasets

source('code/utils_sim.R')

library(survival)
library(haven)
library(readr)
library(tidyverse)


# Universal parts ---------------------------------------------------------

# censoring parameters
# censor.func = function(data) -2.5 + 0.1*data$X1 + 0.05*data$X18
censor.func = function(data) -4 + 0.25*data$X16 + 0.2*data$X18
censor.lambda = 0.01  # scale
censor.gamma = 0.03  # shape

# number of sims
num_sims = 1000




# Scenario 0 --------------------------------------------------------------

# sample size
n = 2500

# generate survival times
lambda = 0.002  # scale
gamma = 0.02  # shape
coarsen = TRUE

# linear part under different treatment rules
event.func.0 = function(data) -1 - 2*data$X3 + 0.1*data$X18 + 0.5*data$X19
event.func.1 = function(data) -1.1 - 2*data$X3 + 0.1*data$X18 + 0.5*data$X19
event.func = function(data) ifelse(data$treat==0, event.func.0(data), event.func.1(data))

# create sims
set.seed(6651709)
simlist = lapply(1:num_sims, function(b) GenSims(n, censor.func, censor.lambda, censor.gamma,
                                                 event.func, lambda, gamma,
                                                 end=61, coarsen=coarsen, corr=TRUE))

# save data
saveRDS(simlist, 'data/simlist_sc0.rds')


# large sample version
n = 10000
set.seed(21347891)
simlist = lapply(1:num_sims, function(b) GenSims(n, censor.func, censor.lambda, censor.gamma,
                                                 event.func, lambda, gamma,
                                                 end=61, coarsen=coarsen, corr=TRUE))
saveRDS(simlist, 'data/simlist_sc4.rds')


# Scenario 1 --------------------------------------------------------------

# sample size
n = 2500

# generate survival times
lambda = 0.002  # scale
gamma = 0.02  # shape
coarsen = TRUE

# linear part under different treatment rules
# event.func.0 = function(data) 2 - 2*data$X3
# event.func.1 = function(data) 2 - 2*data$X3 - 7.5 + 0.005*data$X18 + 1.25*data$X19
event.func.0 = function(data) -1 - 2*data$X3 + 0.1*data$X18 + 0.5*data$X19
event.func.1 = function(data) -1 - 2*data$X3 - 1.9 - 0.1*data$X18 + 1.2*data$X19
event.func = function(data) ifelse(data$treat==0, event.func.0(data), event.func.1(data))

# create sims
set.seed(11028755)
simlist = lapply(1:num_sims, function(b) GenSims(n, censor.func, censor.lambda, censor.gamma,
                                                 event.func, lambda, gamma,
                                                 end=61, coarsen=coarsen, corr=TRUE))

# save data
saveRDS(simlist, 'data/simlist_sc1.rds')


# descriptions
data = simlist[[1]]
true_treat_rule = as.numeric(-1.9 - 0.2*data$X18 + 0.7*data$X19 < 0)
table(true_treat_rule, as.numeric(event.func.1(data) < event.func.0(data)))  # check if treat rule above correct


# large sample version
n = 10000
set.seed(4720922)
simlist = lapply(1:num_sims, function(b) GenSims(n, censor.func, censor.lambda, censor.gamma,
                                                 event.func, lambda, gamma,
                                                 end=61, coarsen=coarsen, corr=TRUE))
saveRDS(simlist, 'data/simlist_sc5.rds')


# Scenario 2 --------------------------------------------------------------

# sample size
n = 2500

# generate survival times
lambda = 0.002  # scale
gamma = 0.02  # shape
coarsen = TRUE

# linear part under different treatment rules
event.func.0 = function(data) -1.5 - 0.2*data$X16 + 0.2*data$X17 - 0.2*data$X18 + 0.5*data$X19 + 1.1*rowSums(data %>% select(X1:X15))
event.func.1 = function(data) -1.6 - 0.1*data$X16 + 0.1*data$X17 - 0.1*data$X18 + 0.4*data$X19 + 1.4*rowSums(data %>% select(X1:X6, X13:X15)) - 0.1*rowSums(data %>% select(X7:X12))
event.func = function(data) ifelse(data$treat==0, event.func.0(data), event.func.1(data))

# create sims
set.seed(776812)
simlist = lapply(1:num_sims, function(b) GenSims(n, censor.func, censor.lambda, censor.gamma,
                                                 event.func, lambda, gamma,
                                                 end=61, coarsen=coarsen, corr=TRUE))

# save data
saveRDS(simlist, 'data/simlist_sc2.rds')


# descriptions
data = simlist[[1]]
true_treat_rule = as.numeric(-0.1 + 0.1*data$X16 - 0.1*data$X17 + 0.1*data$X18 - 0.1*data$X19 + 0.3*rowSums(data %>% select(X1:X6, X13:X15)) - 1.2*rowSums(data %>% select(X7:X12)) < 0)
table(true_treat_rule, as.numeric(event.func.1(data) < event.func.0(data)))  # check if treat rule above correct


# large sample version
n = 10000
set.seed(381939)
simlist = lapply(1:num_sims, function(b) GenSims(n, censor.func, censor.lambda, censor.gamma,
                                                 event.func, lambda, gamma,
                                                 end=61, coarsen=coarsen, corr=TRUE))
saveRDS(simlist, 'data/simlist_sc6.rds')


# Scenario 3 --------------------------------------------------------------


# sample size
n = 2500

# generate survival times
lambda = 0.002  # scale
gamma = 0.02  # shape
coarsen = TRUE

# linear part under different treatment rules
event.func.0 = function(data) -0.4 - 0.3*data$X18 + 0.8*data$X19
event.func.1 = function(data) -0.4 - 0.3*data$X18 + 0.8*data$X19 - 0.95 + (data$X18-6)^2 + (data$X19-4)^2
event.func = function(data) ifelse(data$treat==0, event.func.0(data), event.func.1(data))

# create sims
set.seed(723791)
simlist = lapply(1:num_sims, function(b) GenSims(n, censor.func, censor.lambda, censor.gamma,
                                                 event.func, lambda, gamma,
                                                 end=61, coarsen=coarsen, corr=TRUE))

# save data
saveRDS(simlist, 'data/simlist_sc3.rds')


# descriptions
data = simlist[[1]]
true_treat_rule = as.numeric(- 0.95 + (data$X18-6)^2 + (data$X19-4)^2 < 0)
table(true_treat_rule, as.numeric(event.func.1(data) < event.func.0(data)))  # check if treat rule above correct


# large sample version
n = 10000
set.seed(8889320)
simlist = lapply(1:num_sims, function(b) GenSims(n, censor.func, censor.lambda, censor.gamma,
                                                 event.func, lambda, gamma,
                                                 end=61, coarsen=coarsen, corr=TRUE))
saveRDS(simlist, 'data/simlist_sc7.rds')

