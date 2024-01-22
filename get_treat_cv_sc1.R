source('code/utils_sim.R')

library(survival)
library(readr)
library(dplyr)
library(grf)
require(WeightSVM)
require(modelObj)
require(DynTxRegime)
require(BART)


# collect arguments passed in from job array
args = commandArgs(trailingOnly=TRUE)
sim = as.numeric(args[1])
scenario = 1
data_string = paste0('data/simlist_sc', scenario, '.rds')
save_string = paste0('results/cv/cv_sc', scenario, '_sim', sim, '.rds')

# set seed
set.seed(54623 + sim)

# read in data
data = readRDS(data_string)[[sim]]

# covariate matrix
covs = paste0('X', 1:19)
X = data[, covs] %>% as.matrix()

# for RMST, last time point is considered not censored
data = data %>%
  mutate(r_delta = ifelse(obs_times>=60, 1, delta))


# RIST --------------------------------------------------------------------

# already computed in get_treat_whole_*.R files

times_pred = readRDS(paste0('results/rist/sc', scenario, '_sim', sim, '_rist.rds'))
data$rist = times_pred


# Set up ------------------------------------------------------------------

# set folds and repetitions
M = 1  # number of repetitions
K = 5  # number of folds


# Estimate propensity score and censoring model ---------------------------

# estimate propensity scores
prop_frmla = as.formula(paste0('treat ~', paste0(covs, collapse='+')))
prop_model = glm(prop_frmla, data = data, family = binomial)
prop_scores = predict(prop_model, type = 'response')

# censoring weights
frmla = paste0('Surv(obs_times, ltfu) ~ treat +', paste0(covs, collapse='+'))
censor_frmla = as.formula(frmla)
censor_model = coxph(censor_frmla, data = data)
tmp = survfit(censor_model, newdata = data)
cens_weights = summary(tmp, times=1:60, extend=TRUE)$surv %>% t()


# Ridge -------------------------------------------------------------------

ridge = cv_estimator(data, covs=covs, method='cox', nfolds=K, reps=M, extra_args='ridge')


# Lasso -------------------------------------------------------------------

lasso = cv_estimator(data, covs=covs, method='cox', nfolds=K, reps=M, extra_args='lasso')


# Elastic net -------------------------------------------------------------

elastic = cv_estimator(data, covs=covs, method='cox', nfolds=K, reps=M, extra_args='elastic')


# CSF ---------------------------------------------------------------------

csf = cv_estimator(data, covs=covs, method='csf', nfolds=K, reps=M, prop_scores=prop_scores)


# Genetic algorithm -------------------------------------------------------

# format data
dat = cbind(treat=data$treat, X, delta=data$delta, obs_times=data$obs_times) %>% as.data.frame()

# read whole data results to get initial values
whole_genetic = readRDS(paste0('results/whole/whole_sc', scenario, '_sim', sim, '.rds'))$genetic

# set initial values
init_vals = whole_genetic$eta %>% round(0)

# estimate DTR
genetic = cv_estimator(data, covs=covs, method='genetic', nfolds=K, reps=M,
                       prop_scores=prop_scores, censor=cens_weights,
                       extra_args=list(init_vals=init_vals, max.generations=50))



# Genetic algorithm, augmented --------------------------------------------

# format data
dat = cbind(treat=data$treat, X, delta=data$delta, obs_times=data$obs_times) %>% as.data.frame()

# read whole data results to get initial values
whole_genetic_a = readRDS(paste0('results/whole/whole_sc', scenario, '_sim', sim, '.rds'))$genetic_a

# set initial values
init_vals_a = whole_genetic_a$eta %>% round(0)

# estimate DTR
genetic_a = cv_estimator(data, covs=covs, method='genetic', nfolds=K, reps=M,
                         prop_scores=prop_scores, censor=cens_weights,
                         extra_args=list(init_vals=init_vals_a, augmented=TRUE, max.generations=50))



# Genetic algorithm, smooth -----------------------------------------------

# format data
dat = cbind(treat=data$treat, X, delta=data$delta, obs_times=data$obs_times) %>% as.data.frame()

# read whole data results to get initial values
whole_genetic_s = readRDS(paste0('results/whole/whole_sc', scenario, '_sim', sim, '.rds'))$genetic_s

# set initial values
init_vals_s = whole_genetic_s$eta %>% round(0)

# estimate DTR
genetic_s = cv_estimator(data, covs=covs, method='genetic', nfolds=K, reps=M,
                         prop_scores=prop_scores, censor=cens_weights,
                         extra_args=list(init_vals=init_vals_s, smooth=TRUE, max.generations=50))


# Genetic algorithm, augmented + smooth -----------------------------------

# format data
dat = cbind(treat=data$treat, X, delta=data$delta, obs_times=data$obs_times) %>% as.data.frame()

# read whole data results to get initial values
whole_genetic_as = readRDS(paste0('results/whole/whole_sc', scenario, '_sim', sim, '.rds'))$genetic_as

# set initial values
init_vals_as = whole_genetic_as$eta %>% round(0)

# estimate DTR
genetic_as = cv_estimator(data, covs=covs, method='genetic', nfolds=K, reps=M,
                          prop_scores=prop_scores, censor=cens_weights,
                          extra_args=list(init_vals=init_vals_as, smooth=TRUE, augmented=TRUE, max.generations=50))


# OWL ---------------------------------------------------------------------

# estimate DTR
owl_rist = cv_estimator(data, method='owl', nfolds=K, reps=M, covs=covs,
                        prop_scores=prop_scores, extra_args=list('rist', kernel='linear'))
owl_rist_rad = cv_estimator(data, method='owl', nfolds=K, reps=M, covs=covs,
                            prop_scores=prop_scores, extra_args=list('rist', kernel='radial'))


# RWL ---------------------------------------------------------------------

# estimate DTR
rwl = cv_estimator(data, method='rwl', covs=covs, nfolds=K, reps=M)


# BART --------------------------------------------------------------------

bart = cv_estimator(data, method='bart', covs=covs, nfolds=K, reps=M)

# Save results ------------------------------------------------------------

saveRDS(mget(c('ridge', 'lasso', 'elastic',
               'genetic', 'genetic_a', 'genetic_s', 'genetic_as',
               'owl_rist', 'owl_rist_rad', 'rwl', 'bart')),
        file = save_string)
