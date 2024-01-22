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
scenario = 2
data_string = paste0('data/simlist_sc', scenario, '.rds')
save_string = paste0('results/whole/whole_sc', scenario, '_sim', sim, '.rds')

# set seed
set.seed(24221 + sim)

# read in data
data = readRDS(data_string)[[sim]]

# covariate matrix
covs = paste0('X', 1:19)
X = data[, covs] %>% as.matrix()

# for RMST, last time point is considered not censored
data = data %>%
  mutate(r_delta = ifelse(obs_times>=60, 1, delta))


# RIST --------------------------------------------------------------------

# impute censored values using RIST
source('code/RISTfunctions.r')

# set up data
frmla = as.formula(paste0('~-1+treat*(', paste0(covs, collapse='+'), ')') )
rist_data = model.matrix(frmla, data = data)
rist_data = as.data.frame(cbind(rist_data, delta=data$r_delta, obs_times=data$obs_times))

# parameters
P = ncol(rist_data)-2   # number of dimension for X
K = floor(P/3)  # number of covariate considered per spilt
nmin = 150   # minimum number of observed data in each node
M = 100     # number of trees in each fold
L = 2    # number of folds
tao = 61    # length of study

# split dataset into censored
rist_data_c = rist_data[data$r_delta==0, ]

# RIST
R_Muti_ERT_build = Muti_ERT_fit(as.matrix(rist_data), M, K, L, nmin, SupLogRank=1, tao=tao, impute="random")
R_Muti_ERT_predict= Muti_ERT_Predict(as.matrix(rist_data_c), R_Muti_ERT_build$Forest_seq[[L]], R_Muti_ERT_build$SurvMat_seq[[L]], R_Muti_ERT_build$time_intrest)

# predictions
times_pred_c = R_Muti_ERT_predict$TMT_predict
times_pred = data$obs_times
times_pred[data$r_delta==0] = times_pred_c


# add imputed survival times (i.e. now complete case) to dataset
data$rist = times_pred

# save rist
saveRDS(times_pred, paste0('results/rist/sc', scenario, '_sim', sim, '_rist.rds'))

# times_pred = readRDS(paste0('results/rist/sc', scenario, '_sim', sim, '_rist.rds'))
# data$rist = times_pred



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



# Cox ---------------------------------------------------------------------


# formatting / prep work
cox_frmla = as.formula(paste0('~treat*(', paste0(covs, collapse='+'), ')') )

# get treatment rule
whole_ridge = get_cox_dtr(data, frmla=cox_frmla, lambda='min', penalty='ridge', standardize=TRUE, tpt=60, test_data=NULL)
whole_lasso = get_cox_dtr(data, frmla=cox_frmla, lambda='min', penalty='lasso', standardize=TRUE, tpt=60, test_data=NULL)
whole_elastic = get_cox_dtr(data, frmla=cox_frmla, lambda='min', penalty='elastic', standardize=TRUE, tpt=60, test_data=NULL)



# Genetic -----------------------------------------------------------------


# format data
dat = cbind(treat=data$treat, X, delta=data$delta, obs_times=data$obs_times) %>% as.data.frame()

# run on whole data
whole_genetic = get_genetic_dtr(data=dat, smooth=FALSE, augmented=FALSE, covs=covs,
                                prop_scores=prop_scores, censor=cens_weights,
                                max.generations=50)
whole_genetic_a = get_genetic_dtr(data=dat, smooth=FALSE, augmented=TRUE, covs=covs,
                                  prop_scores=prop_scores, censor=cens_weights,
                                  max.generations=50)
whole_genetic_s = get_genetic_dtr(data=dat, smooth=TRUE, augmented=FALSE, covs=covs,
                                  prop_scores=prop_scores, censor=cens_weights,
                                  max.generations=50)
whole_genetic_as = get_genetic_dtr(data=dat, smooth=TRUE, augmented=TRUE, covs=covs,
                                   prop_scores=prop_scores, censor=cens_weights,
                                   max.generations=50)


# CSF ---------------------------------------------------------------------


whole_csf = get_csf_dtr(data, X, prop_scores, tpt=60, test_data=NULL)



# OWL ---------------------------------------------------------------------


whole_owl_lin = get_owl_dtr(data, covs, prop_scores, kernel='linear')
whole_owl_rad = get_owl_dtr(data, covs, prop_scores, kernel='radial')

whole_rwl = get_rwl_dtr(data, covs)

# don't include svm model
whole_owl_lin[['svm_model']] = NULL
whole_owl_rad[['svm_model']] = NULL
whole_rwl[['svm_model']] = NULL


# BART --------------------------------------------------------------------

x.train = cbind(treat=data$treat, X)
whole_bart = get_bart_dtr(data, x.train)


# Save --------------------------------------------------------------------

# put results in a list
results = mget(grep('^whole', ls(), value=TRUE))
names(results) = gsub('whole_', '', names(results))

# save
saveRDS(results, file=save_string)
