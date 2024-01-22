
# Simulation --------------------------------------------------------------


# generates simulation dataset
GenSims <- function(n, censor.func, censor.lambda, censor.gamma,
                    event.func, event.lambda, event.gamma, 
                    end=61, coarsen=TRUE, corr=TRUE, seed=NULL) {
  
  require(dplyr)
  
  # set seed
  if (!is.null(seed)) (set.seed(seed))
  
  # two ways of generating covariate data:
  # all features independent, or latent variable construction
  data = GenCovs(n, corr=corr, seed=seed)
  
  
  # gen censoring times
  censor.linear.part = censor.func(data)
  data$censor_times = GenSimCens(censor.linear.part, censor.lambda, censor.gamma, coarsen=coarsen)
  
  
  # gen event data
  event.linear.part = event.func(data)
  event_data = GenSimTimes(event.linear.part, event.lambda, event.gamma,
                           data$censor_times, end=end, coarsen=coarsen)
  data$true_times = event_data$true_times
  data$delta = event_data$delta
  data$obs_times = event_data$obs_times
  data$ltfu = as.numeric((ceiling(data$censor_times) <= ceiling(data$true_times)) & (data$censor_times <= end))
  
  # add id column
  data$id = 1:nrow(data)
  
  # sort by observed time
  data = data %>% arrange(obs_times)
  
  # output data frame
  return(data)
  
}


# generates covariates and treatment
GenCovs <- function(n, corr, seed=NULL) {
  
  require(stringr)
  require(dplyr)
  
  # set seed
  if (!is.null(seed)) (set.seed(seed))
  
  if (corr==FALSE) {
    
    # create covariate data
    # match margins on na-accord data
    data = tibble(X1 = rbinom(n, size=1, prob=0.1),
                  X2 = rbinom(n, size=1, prob=0.484),
                  X3 = rbinom(n, size=1, prob=0.211),
                  X4 = rbinom(n, size=1, prob=0.076),
                  X5 = rbinom(n, size=1, prob=0.04),
                  X6 = rbinom(n, size=1, prob=0.106),
                  X7 = rbinom(n, size=1, prob=0.119),
                  X8 = rbinom(n, size=1, prob=0.089),
                  X9 = rbinom(n, size=1, prob=0.052),
                  X10 = rbinom(n, size=1, prob=0.163),
                  X11 = rbinom(n, size=1, prob=0.072),
                  X12 = rbinom(n, size=1, prob=0.046),
                  X13 = rbinom(n, size=1, prob=0.125),
                  X14 = rbinom(n, size=1, prob=0.435),
                  X15 = rbinom(n, size=1, prob=0.13),
                  X16 = rnorm(n, mean=log(40), sd=0.5),
                  X17 = rnorm(n, mean=log(25.1), sd=0.2),
                  X18 = rnorm(n, mean=log(332), sd=1),
                  X19 = rnorm(n, mean=log10(38440), sd=1))
    
  }
  
  if (corr==TRUE) {
    
    # latent variable construction of correlated discrete covs
    uh_prob = 0.2
    UH = rbinom(n, size=1, prob=uh_prob)
    UH[UH==0] = -1
    UH_factor = ifelse(UH==1, UH / (uh_prob/(1-uh_prob)), UH / ((1-uh_prob)/uh_prob))
    UH_factor = UH*plogis(UH_factor)
    data = tibble(X1 = rbinom(n, size=1, prob=plogis(qlogis(0.1) + UH_factor*1.25)),
                  X2 = rbinom(n, size=1, prob=plogis(qlogis(0.484) + UH_factor*1.25)),
                  X3 = rbinom(n, size=1, prob=plogis(qlogis(0.211) + UH_factor*1.25)),
                  X4 = rbinom(n, size=1, prob=plogis(qlogis(0.076) + UH_factor*1.25)),
                  X5 = rbinom(n, size=1, prob=plogis(qlogis(0.04) + UH_factor*1.25)),
                  X6 = rbinom(n, size=1, prob=plogis(qlogis(0.106) + UH_factor*1.25)),
                  X7 = rbinom(n, size=1, prob=plogis(qlogis(0.119) + UH_factor*1.25)),
                  X8 = rbinom(n, size=1, prob=plogis(qlogis(0.089) + UH_factor*1.25)),
                  X9 = rbinom(n, size=1, prob=plogis(qlogis(0.052) + UH_factor*1.25)),
                  X10 = rbinom(n, size=1, prob=plogis(qlogis(0.163) + UH_factor*1.25)),
                  X11 = rbinom(n, size=1, prob=plogis(qlogis(0.072) + UH_factor*1.25)),
                  X12 = rbinom(n, size=1, prob=plogis(qlogis(0.046) + UH_factor*1.25)))
    
    # add unhealthy indicator to data and sort by it
    data$UH = UH
    data = data %>% arrange(UH)
    
    # latent variable construction of correlated continuous covs
    corrMat = matrix(c(1, 0.2, -0.2,
                       0.2, 1, -0.4,
                       -0.2, -0.4, 1),
                     nrow = 3)
    sdevs = diag(c(0.2, 0.5, 1))
    Sigma = sdevs %*% corrMat %*% sdevs
    Xc_uh = MASS::mvrnorm(sum(UH==1), mu=c(log(25.1 - 5), log(332 - 200), log10(38440 + 10000)), Sigma=Sigma)
    Xc_h = MASS::mvrnorm(sum(UH==-1), mu=c(log(25.1 + 1.25), log(332 + 50), log10(38440 - 2000)), Sigma=Sigma)
    Xc = rbind(Xc_h, Xc_uh)
    colnames(Xc) = c('X17', 'X18', 'X19')
    
    # put together
    data = bind_cols(data, as.data.frame(Xc))
    
    # add in other covariates
    data = data %>%
      mutate(X16 = rnorm(n, mean=log(40), sd=0.4),
             X13 = rbinom(n, size=1, prob=0.125),
             X14 = rbinom(n, size=1, prob=0.435),
             X15 = rbinom(n, size=1, prob=0.13))
    
    # put columns in order
    data = data %>% select(str_sort(colnames(data), numeric=TRUE))
    
  }
  
  # treatment assignment mechanism
  data$treat = rbinom(n, size=1, prob=plogis(data$X13 - data$X14 - data$X15 - data$X1))
  
  # output data frame
  return(data)
  
}


GenSimTimes <- function(linear.part, lambda, gamma, censor, end, coarsen=TRUE, seed=NULL) {
  
  # simulates survival times using a linear function of data
  # generates data using a gompertz distribution with user-given parameters
  
  # set seed
  if (!is.null(seed)) (set.seed(seed))
  
  # generate true survival times
  true_times = (1/gamma) * log(((-gamma*log(runif(length(linear.part),0,1)))/(lambda*exp(linear.part))) + 1)
  
  # indicators
  obs_times = pmin(true_times, censor, end)
  
  # coarsen
  if (coarsen) {
    delta = as.numeric((ceiling(true_times) <= ceiling(censor)) & (true_times <= end))
    obs_times = ceiling(obs_times)
  }
  else {
    delta = as.numeric((true_times <= censor) & (true_times <= end))
  }
  
  # return
  return(list(true_times = true_times, delta = delta, obs_times = obs_times))
}


GenSimCens <- function(linear.part, lambda, gamma, coarsen=TRUE, seed=NULL) {
  
  # simulates censoring times using a linear function of data
  # generates data using a gompertz distribution with user-given parameters
  
  # set seed
  if (!is.null(seed)) (set.seed(seed))
  
  # generate true survival times
  censor_times = (1/gamma) * log(((-gamma*log(runif(length(linear.part),0,1)))/(lambda*exp(linear.part))) + 1)
  
  return(censor_times)
}


GenSurvCurv <- function(linear.part, lambda, gamma, tpt=1:60) {
  
  # generate true survival probabilities
  sapply(tpt, function(x) exp( (-lambda*(exp(gamma*x)-1) / gamma) * exp(linear.part) ) )
}


create_plot <- function(n, censor.func, censor.lambda, censor.gamma,
                        event.func.0, event.func.1, lambda, gamma,
                        end=60, coarsen=TRUE, corr=TRUE,
                        plot_all_cens=FALSE) {
  
  # create covariate data
  data = GenCovs(n, corr)
  
  
  # generate censor times from Gompertz distribution
  censor.linear.part = censor.func(data)
  
  # plot true censor curves
  true_cens = GenSurvCurv(censor.linear.part, censor.lambda, censor.gamma)
  plot(1:60, apply(true_cens, 2, mean), type='l', xlab='time', ylab='survival', main='Censoring survival curve', ylim=c(0,1))
  lines(1:60, true_cens[first(which(data$treat==0)),], type='l', col='blue')
  lines(1:60, true_cens[first(which(data$treat==1)),], type='l', col='blue')
  
  if (plot_all_cens) {
    for (i in 1:nrow(true_cens)) {
      lines(1:60, true_cens[i,], type='l', col='blue')
    }
  }
  
  # plot true risk curves
  linear.part.0 = event.func.0(data)
  linear.part.1 = event.func.1(data)
  linear.part.d = pmin(linear.part.0, linear.part.1)
  
  true_surv_0 = GenSurvCurv(linear.part.0, lambda, gamma) %>% apply(2, mean)
  true_surv_1 = GenSurvCurv(linear.part.1, lambda, gamma) %>% apply(2, mean)
  true_surv_d = GenSurvCurv(linear.part.d, lambda, gamma) %>% apply(2, mean)
  
  plot.df = tibble(true_surv_0, true_surv_1, true_surv_d, time=1:60) %>%
    pivot_longer(cols = !time,
                 names_to = 'dtr',
                 values_to = 'surv_prob',
                 names_prefix = 'true_surv_')
  gg = ggplot(plot.df, aes(x=time, y=surv_prob, color=dtr)) +
    geom_line() +
    lims(y = c(0,1)) +
    labs(x = 'Time',
         y = 'Survival probability') +
    scale_color_grey('', labels = c('0', '1', 'Optimal')) +
    theme(legend.position = 'bottom')
  
  print(gg)
  
  # plot(1:60, true_surv_0, col='red', type='l', xlab='time', ylab='risk', main='True risk curves', ylim=c(0,1))
  # lines(1:60, true_surv_1, col='blue', type='l')
  # lines(1:60, true_surv_d, col='purple', type='l')
  # legend('bottomright', legend = c('treat 0 only', 'treat 1 only', 'dtr'), col = c('red', 'blue', 'purple'), lty = 1)
  
  print('treat dtr minus treat 0')
  print(round(true_surv_d[60] - true_surv_0[60], 2))
  print('treat dtr minus treat 1')
  print(round(true_surv_d[60] - true_surv_1[60], 2))
  
  return(plot.df)
  
}



# Dynamic treatment rules -------------------------------------------------

# functions in this section take in data and output treatment rules
# the output is a list where `treat_rule` is the treatment recommendations
# and `eta` are the coefficients to the treatment rule, if applicable


get_cox_dtr <- function(data, frmla, lambda='1se', penalty='lasso', standardize=TRUE, tpt=60, test_data=NULL) {
  
  require(survival)
  require(glmnet)
  require(dplyr)
  
  # set alpha by penalty
  if (penalty=='lasso') (alpha=1)
  else if (penalty=='elastic') (alpha=0.5)
  else if (penalty=='ridge') (alpha=0)
  else (stop('invalid penalty'))
  
  # covariate matrix
  X = model.matrix(frmla, data = data)
  X = X[,-1]  # remove intercept
  
  # Cox model
  y = Surv(data$obs_times, data$delta)
  cox.fit.reg = glmnet(X, y, family = 'cox', alpha = alpha, standardize = standardize)
  
  # cv to tune lambda
  cox.cvfit.reg = cv.glmnet(X, y, family = 'cox', type.measure = 'C', alpha = alpha, standardize = standardize)
  if (lambda=='1se') (lambda = cox.cvfit.reg$lambda.1se)
  else if (lambda=='min') (lambda = cox.cvfit.reg$lambda.min)
  else (stop('invalid lambda'))
  
  # get coefficients
  coefs = as.numeric(coef(cox.fit.reg, s=lambda))
  names(coefs) = row.names(coef(cox.fit.reg, s=lambda))
  
  # if input test data then get treatment recommendations on test data
  if (!is.null(test_data)) (data = test_data)
  
  # covariate matrices for ZOM
  X0 = model.matrix(frmla, data = data %>% mutate(treat=0))
  X0 = X0[,-1]  # remove intercept
  
  X1 = model.matrix(frmla, data = data %>% mutate(treat=1))
  X1 = X1[,-1]  # remove intercept
  
  # survival curves under ZOMs
  tmp0 = survfit(cox.fit.reg, s = lambda, x = X, y = y, newx = X0)
  tmp1 = survfit(cox.fit.reg, s = lambda, x = X, y = y, newx = X1)
  cox.surv0 = summary(tmp0, times = 1:tpt, extend = TRUE)$surv
  cox.surv1 = summary(tmp1, times = 1:tpt, extend = TRUE)$surv
  
  # get treatment rule
  if (is.null(nrow(cox.surv0))) {
    treat_rule = as.numeric(cox.surv1[tpt] >= cox.surv0[tpt])
  }
  else {
    treat_rule = as.numeric(cox.surv1[tpt,] >= cox.surv0[tpt,])
  }
  
  # return treatment rule and coefficients
  return(list(treat_rule=treat_rule, eta=coefs))
  
}


get_csf_dtr <- function(data, X, prop_scores, tpt=60, test_data=NULL) {
  
  require(grf)
  
  # propensity score
  prop_scores = data$treat*prop_scores + (1-data$treat)*(1-prop_scores)
  
  # train CSF
  cs.forest = causal_survival_forest(X, data$obs_times, data$treat, data$r_delta,
                                     W.hat = prop_scores,
                                     num.trees = 5000, 
                                     target = 'survival.probability', 
                                     horizon = tpt,
                                     honesty=TRUE)
  
  # if input test data then get treatment recommendations on test data
  if (!is.null(test_data)) (X = test_data)
  
  # predict
  cs.pred = predict(cs.forest, newdata = X, estimate.variance = FALSE)
  
  # get treatment assignment
  treat_rule = as.numeric(cs.pred$predictions >= 0)
  
  return(list(treat_rule=treat_rule, eta=NA))
  
}


get_bart_dtr <- function(data, x.train, ndpost=100, thin=15, burn=500, tpt=60, test_data=NULL) {
  
  # treatment indicator required to be first column of x.train
  
  require(BART)
  
  n = nrow(data)
  num_times = length(unique(data$obs_times))
  num_obstimes = n * num_times
  
  # run BART
  post = surv.bart(x.train=x.train, times=data$obs_times, delta=data$delta,
                   x.test=x.train, ndpost=ndpost, keepevery=thin, nskip=burn)
  
  # if input test data then get treatment recommendations on test data
  if (!is.null(test_data)) {
    
    # get prediction covariate matrices
    # covariates for test
    x.test = model.matrix(as.formula(paste0('~ ', paste0(covs, collapse='+'))),
                          data = test_data)
    x.test = x.test[,-1]
    x.test = cbind(newid=1:nrow(x.test), treat=test_data$treat, x.test)
    # pre = surv.pre.bart(times=test_data$obs_times, delta=test_data$delta,
    #                     x.train=x.test, x.test=x.test)
    
    tx.test.0 = do.call(rbind, replicate(num_times, x.test, simplify=FALSE))
    tx.test.0 = cbind(t=sort(rep(1:num_times, nrow(test_data))), tx.test.0)
    tx.test.0 = tx.test.0[order(tx.test.0[,'newid']),]
    tx.test.0 = tx.test.0[,colnames(tx.test.0)!='newid']
    
    tx.test.1 = tx.test.0
    
    
    # reset num_obstimtes
    num_obstimes = nrow(test_data) * num_times
    
  }
  
  else {
    
    x.train = cbind(newid=1:nrow(x.train), x.train)
    # # get prediction covariate matrices
    # pre = surv.pre.bart(times=data$obs_times, delta=data$delta,
    #                     x.train=x.train, x.test=x.train)
    
    tx.test.0 = do.call(rbind, replicate(num_times, x.train, simplify=FALSE))
    tx.test.0 = cbind(t=sort(rep(1:num_times, n)), tx.test.0)
    tx.test.0 = tx.test.0[order(tx.test.0[,'newid']),]
    tx.test.0 = tx.test.0[,colnames(tx.test.0)!='newid']
    
    tx.test.1 = tx.test.0
    
  }
  
  # tx.test.0 = tx.test.1 = pre$tx.test
  tx.test.0[, 2] = 0
  tx.test.1[, 2] = 1
  
  # make predictions
  pred.0 = predict(post, newdata=tx.test.0)
  pred.1 = predict(post, newdata=tx.test.1)
  
  # estimate DTR
  h = seq(tpt, num_obstimes, by=num_times)
  treat_rule = as.numeric(pred.1$surv.test.mean[h] - pred.0$surv.test.mean[h] > 0)
  
  if (any(is.na(treat_rule))) {
    saveRDS(list(data, post, treat_rule, test_data), 'temp.rds')
    stop
  }
  
  return(list(treat_rule=treat_rule, eta=NA))
  
}


get_owl_dtr <- function(data, covs, prop_scores, kernel='linear', 
                        costs=NULL, gamma=NULL, tpt=60, test_data=NULL) {
  
  require(survival)
  require(WeightSVM)
  require(dplyr)
  
  # svm formula
  svm_frmla = as.formula(paste0('~ ', paste0(covs, collapse='+')))
  
  # costs to consider
  if (is.null(costs)) (costs = 2^(seq(-8, 2, by=4)))
  else (costs = costs)
  
  # gamma to consider
  if (is.null(gamma) & kernel=='radial') (gamma = c(1/nrow(data), 2^(-5), 2^(-10)))
  else (gamma = gamma)
  
  # propensity score
  prop_scores = data$treat*prop_scores + (1-data$treat)*(1-prop_scores)
  
  # weights
  times_pred = data$rist
  w = as.numeric(times_pred > tpt) / prop_scores
  
  # prep work
  P = ncol(data)-2
  Y = 2*data$treat - 1
  X = model.matrix(svm_frmla, data = data)
  X = scale(X[,-1])
  
  # svm
  if (length(costs)>1 | length(gamma)>1) {
    
    if (!is.null(gamma)) (search.grid = list(cost = costs, gamma = gamma))
    else (search.grid = list(cost = costs))
    
    svm_model = best.tune_wsvm(train.x=X, train.y=as.factor(data$treat), weight=w, ranges=search.grid,
                               tunecontrol=tune.control(sampling='cross', cross=3), 
                               kernel=kernel, scale=FALSE)
  }
  
  else {
    svm_model = wsvm(x=X, y=as.factor(data$treat), weight=w,
                     kernel=kernel, scale=FALSE, cost=costs, gamma=gamma)
  }
  
  # treatment rule
  treat_rule = predict(svm_model, newdata=X) %>% as.character() %>% as.numeric()
  
  # get coefficients of decision rule
  beta = drop(t(svm_model$coefs) %*% svm_model$SV)
  beta0 = -svm_model$rho
  eta = c(beta0, beta)
  
  # apply to test set
  if (!is.null(test_data)) {
    X_test = model.matrix(svm_frmla, data = test_data)
    X_test = scale(X_test[,-1])
    treat_rule = predict(svm_model, newdata=X_test) %>% as.character() %>% as.numeric()
    return(list(treat_rule=treat_rule, eta=eta))
  }
  
  return(list(treat_rule=treat_rule, eta=eta, svm_model=svm_model))
  
}


get_rwl_dtr <- function(data, covs, costs=NULL, tpt=60, test_data=NULL) {
  
  require(survival)
  require(WeightSVM)
  require(modelObj)
  require(DynTxRegime)
  require(dplyr)
  
  # formula
  frmla = as.formula(paste0('~ ', paste0(covs, collapse='+')))
  
  # costs to consider
  if (is.null(costs)) (costs = 2^(seq(-8, 2, by=2)))
  else (costs = costs)
  
  # reward
  reward = as.numeric(data$rist > tpt)
  
  # define regression models
  moPropen <- modelObj::buildModelObj(model = frmla, 
                                      solver.method = 'glm', 
                                      solver.args = list(family='binomial'), 
                                      predict.method = 'predict.glm', 
                                      predict.args = list(type='response'))
  
  moMain <- modelObj::buildModelObj(model = frmla,
                                    solver.method = 'lm')
  
  # specify class of treatment regimes (linear)
  kernel = 'linear'
  
  # scale data
  data = data %>%
    mutate_at(covs, function(x) as.numeric(scale(x)))
  
  data = as.data.frame(data)
  
  if (length(costs)>1) {
    svm_model <- rwl(moPropen = moPropen,
                     moMain = moMain,
                     data = data,
                     response = reward,
                     txName = 'treat',
                     regime = frmla,
                     lambdas = costs,
                     cvFolds = 5L,
                     kernel = kernel,
                     responseType = 'continuous',
                     verbose = 0L)
  }
  
  else {
    svm_model <- rwl(moPropen = moPropen,
                     moMain = moMain,
                     data = data,
                     response = reward,
                     txName = 'treat',
                     regime = frmla,
                     lambdas = costs,
                     cvFolds = 0L,
                     kernel = kernel,
                     responseType = 'continuous',
                     verbose = 0L)
  }
  
  # apply to test set
  if (!is.null(test_data)) {
    eta = regimeCoef(object = svm_model)
    X_test = model.matrix(frmla, data = test_data)
    X_test = scale(X_test[,-1])
    treat_rule = as.numeric(cbind(1, X_test) %*% eta >= 0)
    return(list(treat_rule=treat_rule, eta=eta))
  }
  
  # return treatment rule
  treat_rule = optTx(x=svm_model, newdata=data)$optimalTx
  
  # get coefficients
  eta = regimeCoef(object = svm_model)
  
  return(list(treat_rule=treat_rule, eta=eta, svm_model=svm_model))
}


get_genetic_dtr <- function(data, prop_scores, censor, covs, smooth=FALSE, augmented=FALSE, 
                            standardize=TRUE, max.generations=100, pop.size=1000, 
                            tpt=60, init_vals=NULL, test_data=NULL) {
  
  # data must be ordered X (treat first column), delta, times
  
  library(compiler)
  library(rgenoud)
  library(survival)
  library(glmnet)
  require(dplyr)
  
  
  # make sure data is sorted
  stopifnot(!is.unsorted(data$obs_times))
  
  # X2 is treatment regime covariates
  X2 = as.matrix(data[,2:(ncol(data)-2)])
  
  
  # augmented
  if (augmented==TRUE) {
    
    # failure model
    cox_frmla = as.formula(paste0('~ treat * (', paste0(covs, collapse='+'), ')'))
    
    X = model.matrix(cox_frmla, data = data)
    X = X[,-1]  # remove intercept
    
    # penalized Cox model
    y = Surv(data$obs_times, data$delta)
    cox.fit.reg = glmnet(X, y, family = 'cox', alpha = 0.5)
    
    # cv to tune lambda
    cox.cvfit.reg = cv.glmnet(X, y, family = 'cox', type.measure = 'C', alpha = 0.5)
    lambda = cox.cvfit.reg$lambda.min
    
    # covariate matrices
    X0 = model.matrix(cox_frmla, data = data %>% mutate(treat=0))
    X0 = X0[,-1]  # remove intercept
    
    X1 = model.matrix(cox_frmla, data = data %>% mutate(treat=1))
    X1 = X1[,-1]  # remove intercept
    
    # surv fit
    failure.fit.0 = survfit(cox.fit.reg, s = lambda, x = X, y = y, newx = X0)
    failure.surv.0 = summary(failure.fit.0, times = 1:tpt, extend = TRUE)$surv %>% t()
    
    failure.fit.1 = survfit(cox.fit.reg, s = lambda, x = X, y = y, newx = X1)
    failure.surv.1 = summary(failure.fit.1, times = 1:tpt, extend = TRUE)$surv %>% t()
    
    
    # failure model hazard
    failure.cumhaz = summary(failure.fit.0, times = 1:tpt, extend = TRUE)$cumhaz %>% t()
    failure.haz0 = failure.cumhaz[,1]
    failure.haz = apply(failure.cumhaz, 1, diff) %>% t()
    failure.haz.0 = cbind(failure.haz0, failure.haz)
    
    failure.cumhaz = summary(failure.fit.1, times = 1:tpt, extend = TRUE)$cumhaz %>% t()
    failure.haz0 = failure.cumhaz[,1]
    failure.haz = apply(failure.cumhaz, 1, diff) %>% t()
    failure.haz.1 = cbind(failure.haz0, failure.haz)
    
    aug_terms = list(failure.surv.0=failure.surv.0,
                     failure.surv.1=failure.surv.1,
                     lambda.surv.0=failure.haz.0,
                     lambda.surv.1=failure.haz.1)
    
    # use algorithm to get eta
    eta = Genetic.IPWE(fn=aipwe_km_opt, X=X2, u=data$obs_times, z=data$treat, delta=data$delta,
                       prop_scores=prop_scores, censor=censor, smooth=smooth, 
                       nvars=ncol(X2)+1, max.generations=max.generations, pop.size=pop.size,
                       aug_terms=aug_terms, init_vals=init_vals, standardize=standardize)
    val_hat = eta[length(eta)]
    eta = eta[1:(length(eta)-1)]
    names(eta) = c('intercept', colnames(X2))
    
  }
  
  else {
    
    # use algorithm to get eta
    eta = Genetic.IPWE(fn=ipwe_km_opt, X=X2, u=data$obs_times, z=data$treat, delta=data$delta,
                       prop_scores=prop_scores, censor=censor, smooth=smooth, 
                       nvars=ncol(X2)+1, max.generations=max.generations, pop.size=pop.size,
                       init_vals=init_vals, standardize=standardize)
    val_hat = eta[length(eta)]
    eta = eta[1:(length(eta)-1)]
    names(eta) = c('intercept', colnames(X2))
    
  }
  
  # treatment rule
  if (standardize) (X2 = scale(X2))
  treat_rule = as.numeric(cbind(1,X2) %*% eta >= 0)
  
  if (!is.null(test_data)) {
    test_data = test_data[, !grepl('drugstart', colnames(test_data))]
    if (standardize) (test_data = scale(test_data))
    treat_rule = as.numeric(cbind(1,test_data) %*% eta >= 0)
    return(list(treat_rule=treat_rule, eta=eta))
  }
  
  return(list(treat_rule=treat_rule, eta=eta, val_hat=val_hat))
}


# CV estimator -----------------------------------------------------


cv_estimator <- function(data, method, covs, prop_scores=NULL, censor=NULL, nfolds=NULL, reps=NULL, extra_args=NULL) {
  
  # computes CV treatment rules
  # see `get_treat_rules.R` for examples of usage 
  
  # Inputs:
  #   - data: data frame with survival times (`obs_times`), event indicator (`delta`), 
  #           treatment indicator (`treat`), and covariates (specified in `covs` argument)
  #   - method: character vector of length 1, which precision medicine method? options:
  #         - "cox": Cox regression
  #         - "genetic": genetic algorithm
  #         - "owl": outcome weighted learning, requires `rist` in `data`, imputed survival times
  #         - "rwl": residual weighted learning, requires `rist` in `data`, imputed survival times
  #         - "csf": causal survival forests
  #   - covs: vector of covariates, assumed to be same for all of propensity score model,
  #           probability of being uncensored model, and treatment rule inputs
  #   - prop_scores: some methods require user to pass in propensity scores
  #   - censor: matrix (length of data x number of unique times), some methods require user to pass 
  #             in matrix of censoring survival probabilities
  #   - nfolds: number of folds in cross-validation
  #   - reps: number of repetitions for cross-validation
  #   - extra_args: named list, some methods require extra arguments
  
  N = nrow(data)
  
  # make sure sorted by obs times
  stopifnot(!is.unsorted(data$obs_times))
  
  # default to LOOCV
  if (is.null(nfolds)) (nfolds = N)
  if (is.null(reps)) (reps = 1)
  
  # for each repetition
  treat_rules_stsh = list()
  eta = list()
  
  for (i in 1:reps) {
    
    # CV indices
    if (nfolds==N) (folds = sample(1:N, size=N, replace=FALSE))
    else (folds = sample(1:nfolds, size = N, replace = TRUE))
    
    # treatment rules within each rep
    treat_rules_stsh[[i]] = list()
    eta[[i]] = list()
    
    for (j in 1:nfolds) {
      
      # set up data
      data.train = data[folds != j, ]
      data.test = data[folds == j, ]
      
      # data.test.order = order(data.test$obs_times)
      # data.test = data.test[data.test.order, ]
      
      # compute treatment rule using training data
      if (method=='cox') {
        
        # type of penalization
        penalties = c('lasso', 'elastic', 'ridge')
        if (any(penalties %in% extra_args)) (penalty=penalties[penalties %in% extra_args])
        else (penalty = 'lasso')
        
        # formatting / prep work
        cox_frmla = as.formula(paste0('~ treat * (', paste0(covs, collapse='+'), ')'))
        
        # get treatment rule
        cox_res = get_cox_dtr(data.train, frmla=cox_frmla, lambda='min', penalty=penalty, tpt=60, test_data=data.test)
        treat_rule = cox_res$treat_rule
        eta[[i]][[j]] = cox_res$eta
        treat_rules_stsh[[i]][[j]] = data.frame(id=data.test$id, treat_rule=treat_rule)
        
      }
      
      if (method=='genetic') {
        
        # requires treatment and censoring weights
        stopifnot(!is.null(prop_scores))
        stopifnot(!is.null(censor))
        
        # smooth or not smooth
        if ('smooth' %in% extra_args) (smooth=TRUE)
        else (smooth=FALSE)
        
        # augmented or not augmented
        if ('augmented' %in% extra_args) (augmented=TRUE)
        else (augmented=FALSE)
        
        # initial values to genetic algorithm
        if ('init_vals' %in% names(extra_args)) (init_vals = extra_args$init_vals)
        else (init_vals=NULL)
        
        # max generations to genetic algorithm
        if ('max.generations' %in% names(extra_args)) (max.generations = extra_args$max.generations)
        else (max.generations=100)
        
        # population size to genetic algorithm
        if ('pop.size' %in% names(extra_args)) (pop.size = extra_args$pop.size)
        else (pop.size=1000)
        
        # formatting / prep work
        X = model.matrix(as.formula(paste0('~ ', paste0(covs, collapse='+'))),
                         data = data)
        X = X[,-1]
        
        X.train = X[folds !=j, ]
        X.test = X[folds == j, ]
        
        dat = cbind(treat=data.train$treat, X.train, delta=data.train$delta, obs_times=data.train$obs_times) %>% as.data.frame()
        
        # get treatment rule
        genetic_res = get_genetic_dtr(data=dat, smooth=smooth, augmented=augmented, prop_scores=prop_scores[folds!=j], censor=censor[folds!=j,],
                                      covs=covs, max.generations=max.generations, pop.size=pop.size, test_data=X.test, init_vals=init_vals)
        treat_rule = genetic_res$treat_rule
        eta[[i]][[j]] = genetic_res$eta
        treat_rules_stsh[[i]][[j]] = data.frame(id=data.test$id, treat_rule=treat_rule)
        
      }
      
      
      if (method=='owl') {
        
        # requires treatment weights
        stopifnot(!is.null(prop_scores))
        
        # kernel
        if ('kernel' %in% names(extra_args)) (kernel = extra_args$kernel)
        else (kernel='linear')
        
        # cost
        if ('costs' %in% names(extra_args)) (costs = extra_args$costs)
        else (costs=NULL)
        
        # run owl w/ rist
        owl = get_owl_dtr(data=data.train, kernel=kernel, covs=covs, prop_scores=prop_scores[folds!=j],
                          costs=costs, tpt=60, test_data=data.test)
        
        # return treatment rule
        treat_rule = owl$treat_rule
        treat_rules_stsh[[i]][[j]] = data.frame(id=data.test$id, treat_rule=treat_rule)
        
        # get coefficients
        eta[[i]][[j]] = as.vector(owl$eta)
        
      }
      
      
      if (method=='rwl') {
        
        # cost
        if ('costs' %in% names(extra_args)) (costs = extra_args$costs)
        else (costs=NULL)
        
        # run rwl w/ rist
        rwl = get_rwl_dtr(data=data.train, covs=covs, costs=costs, tpt=60, test_data=data.test)
        
        # return treatment rule
        treat_rule = rwl$treat_rule
        treat_rules_stsh[[i]][[j]] = data.frame(id=data.test$id, treat_rule=treat_rule)
        
        # get coefficients
        eta[[i]][[j]] = as.vector(rwl$eta)
        
      }
      
      
      if (method=='csf') {
        
        # requires treatment weights
        stopifnot(!is.null(prop_scores))
        
        # covariates for training
        X = model.matrix(as.formula(paste0('~ ', paste0(covs, collapse='+'))),
                         data = data.train)
        X = X[,-1]
        
        # covariates for test
        X.test = model.matrix(as.formula(paste0('~ ', paste0(covs, collapse='+'))),
                              data = data.test)
        X.test = X.test[,-1]
        
        # run csf
        csf = get_csf_dtr(data.train, X, prop_scores=prop_scores[folds!=j], test_data=X.test)
        
        # get treatment rule
        treat_rule = csf$treat_rule
        treat_rules_stsh[[i]][[j]] = data.frame(id=data.test$id, treat_rule=treat_rule)
        
        # leave coefficients blank
        eta[[i]][[j]] = NA
        
      }
      
      if (method=='bart') {
        
        # covariates for training
        X = model.matrix(as.formula(paste0('~ ', paste0(covs, collapse='+'))),
                         data = data.train)
        X = X[,-1]
        
        # run bart
        bart = get_bart_dtr(data.train, cbind(treat=data.train$treat, X), test_data=data.test)
        
        # get treatment rule
        treat_rule = bart$treat_rule
        treat_rules_stsh[[i]][[j]] = data.frame(id=data.test$id, treat_rule=treat_rule)
        
        # leave coefficients blank
        eta[[i]][[j]] = NA
        
      }
      
    }
    
  }
  
  return(list(treat_rule = treat_rules_stsh,
              folds = folds,
              eta = eta))
  
}



# Value estimators --------------------------------------------------------


ipwe_km <- function(u, z, delta, prop_scores, censor, treat_rule, tpt) {
  
  # K-M value estimator from jiang et al.
  
  # Inputs:
  #   - u: event times
  #   - z: treatment indicator
  #   - delta: event indicator
  #   - prop_scores: propensity scores
  #   - censor: matrix, probability of being uncensored
  #   - treat_rule: treatment recommendations
  
  # Outputs:
  #   - vector of survival probabilities with length
  #     equal to number of unique event times
  
  
  # assume all vectors sorted by u
  stopifnot(!is.unsorted(u))
  
  # inputed propensity scores should be P(Z=1|X)
  # change to P(Z=z|X)
  prop_scores2 = prop_scores*(z) + (1-prop_scores)*(1-z)
  follow_treat = z*treat_rule + (1-z)*(1-treat_rule)
  
  # compute quantity inside \prod
  w = follow_treat / (prop_scores2 * censor)  # weights
  num = c()
  den = c()
  for (i in 1:length(tpt)) {
    index = (u==tpt[i])
    num[i] = sum(w[index & delta, i])
    den[i] = sum(w[u >= tpt[i], i])
    if (num[i]==0 & den[i]==0) (den[i]=1)
    # index = (u>=tpt[i]-1 & u<tpt[i])
    # num[i] = sum(w[index & delta, i])
    # den[i] = sum(w[u >= tpt[i]-1, i])
    # if (num[i]==0 & den[i]==0) (den[i]=1)
  }
  
  # compute K-M curve
  curv = cumprod(1-num/den)
  curv[is.na(curv)] = 0

  # if (se.fit) {
  #   Lambda = cumsum(num/den)
  # }

  return(curv)
}


aipwe_km <- function(u, z, delta, prop_scores, censor, failure.surv, lambda.surv, treat_rule) {
  
  # same as ipwe_km but with augmentation
  
  # assume all vectors sorted by u
  stopifnot(!is.unsorted(u))
  
  # inputed propensity scores should be P(Z=1|X)
  # change to P(Z=z|X)
  prop_scores2 = prop_scores*(z) + (1-prop_scores)*(1-z)
  follow_treat = z*treat_rule + (1-z)*(1-treat_rule)
  
  # compute quantity inside \prod
  Ntimes = ncol(censor)
  unique_times = 1:Ntimes
  w = follow_treat / (prop_scores2 * censor)  # weights
  num = c()
  den = c()
  for (i in 1:Ntimes) {
    index = (u==unique_times[i])
    num[i] = sum(w[index, i] * delta[index]) + sum((1-w[index, i]) * failure.surv[index,i] * censor[index,i] * lambda.surv[index,i])
    den[i] = sum(w[u>unique_times[i] | (index & delta==0), i]) + num[i] + sum((1-w[index, i]) * failure.surv[index,i] * censor[index,i])
    if (num[i]==0 & den[i]==0) (den[i]=1)
  }
  
  # compute K-M curve
  curv = cumprod(1-num/den)
  curv[is.na(curv)] = 0
  
  return(curv)
}



ipwe_km_opt <- function(eta, X, obs_times, z, delta, prop_scores, censor, smooth=FALSE) {
  
  # K-M value estimator from jiang et al. - to be passed into Genetic.IPWE
  
  # same as ipwe_km except takes as a function eta instead of treatment rule
  # eta are coefficients to a linear treatment rule with covariates X
  # outputs estimated survival probability at last time point only
  
  # assume all vectors sorted by u
  stopifnot(!is.unsorted(obs_times))
  
  len_eta = length(eta)
  
  # smooth treatment rule?
  if (smooth) {
    sd.etaX <- sd(c(cbind(1,X) %*% eta))
    if (!is.finite(sd.etaX)) return(-1000)
    if (sd.etaX > 0) eta <- eta/sd.etaX else eta <- c(ifelse(eta[1] >= 0, 1, -1), rep(0, len_eta-1))
    treat_rule = pnorm(c(cbind(1,X) %*% eta)/((length(z)/4)^(-1/3)))
  }
  
  else {
    treat_rule = as.numeric((cbind(1,X) %*% eta) >= 0)
  }
  
  # inputed propensity scores should be P(Z=1|X)
  # change to P(Z=z|X)
  prop_scores2 = prop_scores*(z) + (1-prop_scores)*(1-z)
  follow_treat = z*treat_rule + (1-z)*(1-treat_rule)
  
  # compute quantity inside \prod
  Ntimes = ncol(censor)
  unique_times = 1:Ntimes
  w = follow_treat / (prop_scores2 * censor)  # weights
  num = c()
  den = c()
  for (i in 1:Ntimes) {
    index = (obs_times==unique_times[i])
    num[i] = sum(w[index, i] * delta[index])
    den[i] = sum(w[obs_times>unique_times[i] | (index & delta==0), i]) + num[i]
  }
  
  # compute K-M curve
  curv = cumprod(1-num/den)
  curv[is.na(curv)] = 0
  
  # output only last time point
  return(curv[length(curv)])
}


aipwe_km_opt <- function(eta, X, obs_times, z, delta, prop_scores, censor, failure.surv.0, failure.surv.1, lambda.surv.0, lambda.surv.1, smooth=FALSE) {
  
  # K-M value estimator from jiang et al. - to be passed into Genetic.IPWE
  
  # same as aipwe_km except takes as a function eta instead of treatment rule
  # eta are coefficients to a linear treatment rule with covariates X
  # outputs estimated survival probability at last time point only
  
  # assume all vectors sorted by u
  stopifnot(!is.unsorted(obs_times))
  
  len_eta = length(eta)
  
  # smooth treatment rule?
  if (smooth) {
    sd.etaX <- sd(c(cbind(1,X) %*% eta))
    if (!is.finite(sd.etaX)) return(-1000)
    if (sd.etaX > 0) eta <- eta/sd.etaX else eta <- c(ifelse(eta[1] >= 0, 1, -1), rep(0, len_eta-1))
    treat_rule = pnorm(c(cbind(1,X) %*% eta)/((length(z)/4)^(-1/3)))
  }
  
  else {
    treat_rule = as.numeric((cbind(1,X) %*% eta) >= 0)
  }
  
  # inputed propensity scores should be P(Z=1|X)
  # change to P(Z=z|X)
  prop_scores2 = prop_scores*(z) + (1-prop_scores)*(1-z)
  follow_treat = z*treat_rule + (1-z)*(1-treat_rule)
  
  failure.surv = failure.surv.1*treat_rule + failure.surv.0*(1-treat_rule)
  lambda.surv = lambda.surv.1*treat_rule + lambda.surv.0*(1-treat_rule)
  
  # compute quantity inside \prod
  Ntimes = ncol(censor)
  unique_times = 1:Ntimes
  w = follow_treat / (prop_scores2 * censor)  # weights
  num = c()
  den = c()
  for (i in 1:Ntimes) {
    index = (obs_times==unique_times[i])
    num[i] = sum(w[index, i] * delta[index]) + sum((1-w[index, i]) * failure.surv[index,i] * censor[index,i] * lambda.surv[index,i])
    den[i] = sum(w[obs_times>unique_times[i] | (index & delta==0), i]) + num[i] + sum((1-w[index, i]) * failure.surv[index,i] * censor[index,i])
    if (num[i]==0 & den[i]==0) (den[i]=1)
  }
  
  # compute K-M curve
  curv = cumprod(1-num/den)
  curv[is.na(curv)] = 0
  
  # output only last time point
  return(curv[length(curv)])
}


ipw_estimator <- function(u, z, delta, prop_scores, censor, treat_rule, tpt, hajek=FALSE) {
  
  prop_scores2 = prop_scores^(treat_rule) * (1-prop_scores)^(1-treat_rule)
  follow_treat = z^(treat_rule) * (1-z)^(1-treat_rule)
  
  num = follow_treat * delta * as.numeric(u > tpt)
  den = prop_scores2 * censor
  
  num2 = follow_treat
  den2 = prop_scores2 
  
  if (hajek==FALSE) (return(mean(num/den)))
  else (return(sum(num/den) / sum(num2/den2)))
}


# Weighted SVM with L1 norm -----------------------------------------------


wsvmL1 <- function(x, y, w, gamma) {
  
  # code from https://doi.org/10.1007/s10985-016-9376-x
  
  require("Rglpk")
  # pp is the dimension of covariate
  pp <- ncol(x)
  # nn is the sample size
  nn <- nrow(x)
  ## Solve: (see original for more comments)
  ## Create objective vector
  obj <- c(c(w), rep(gamma, pp), rep(0, 1+pp))
  
  G <- matrix(0, nrow=nn+pp+pp, ncol=nn+pp+1+pp)
  G[1:nn, ] <- cbind(diag(nn), matrix(0, nn, pp), y, x*y)
  G[(nn+1):(nn+pp), ] <- cbind(matrix(0, nrow=pp, ncol=nn),
                               cbind(diag(pp), matrix(0, nrow=pp, ncol=1), -diag(pp)))
  G[(nn+pp+1):(nn+2*pp), ] <- cbind(matrix(0, nrow=pp, ncol=nn),
                                    cbind(diag(pp), matrix(0, nrow=pp, ncol=1), diag(pp)))
  ## Set rhs of Gx <> h
  h <- c(rep(1, nn), rep(0, 2*pp))
  ## Set directions of the linear ineq/eq
  dir <- c(rep(">=", nn+2*pp))
  ## Set bounds on u, s, b,
  ## note lower bound defaults to zero and
  ## upper defaults to +Inf for each variable
  ## so we only need to specify the lower bounds on b_0, ..., b_p is -Inf
  lwr = list (ind=c((nn+pp+1):(nn+pp+1+pp)),
              val=rep(-Inf,1+pp))
  bds = list (lower = lwr)
  soln = Rglpk_solve_LP (obj, G, dir, h, bounds=bds, max=F,
                         types=rep("C", nn+1));
  b0 <- soln$solution[(nn+pp+1)]
  b <- soln$solution[(nn+pp+1+1):(nn+pp+1+pp)]
  M <- matrix(c(b0,b), pp+1, 1)
  rownames(M)<- c('intercept', colnames(x))
  sv <- list(sv = M)
  return(sv)
  
}

cv.wsvmL1 <- function(lambdas,dat,nfold) {
  
  # code from https://doi.org/10.1007/s10985-016-9376-x
  
  n <- length(dat$y);
  folds <- split(sample(seq(n)), rep(1:nfold, length = n));
  emat <- matrix(NA, length(lambdas), nfold);
  for (j in 1:length(lambdas)){
    for (i in seq(nfold)) {
      omit <- folds[[i]];
      model <- wsvmL1(x = dat$Xs[-omit,], y = dat$y[-omit],
                      w = dat$weight[-omit], gamma=lambdas[j])
      emat[j, i] <- sum(dat$weight[omit]*
                          (sign(cbind(1,dat$Xs[omit,])%*%model$sv)-dat$y[omit])^2)
    }
  }
  cv <- apply(emat, 1, mean);
  lam <- lambdas[which.min(cv)];
  list(lam=lam, emat=emat);
  
}



# Evaluation --------------------------------------------------------------


extract_surv <- function(treat_list, data, prop_scores, censor, tpt) {
  
  # extract survival probability estimates (do this for each imputation)
  
  # Inputs:
  # treat_list: treat_rule output from cv_estimator function
  # data: data, sorted by obs time
  # prop_scores: propensity scores
  # censor: censoring weights
  # Outputs:
  # single survival curve
  
  
  # index to arrange data in original order (i.e. by id) -> sorted obs_times order
  u.order = rank(data$id, ties='first')
  
  # combine treatment recommendations for each M and arrange
  # i.e. get list of df where each df is treat rule for all data
  treat_list_flat = lapply(treat_list, function(x) bind_rows(lapply(x, function(y) as.data.frame(y))))
  treat_list_flat = lapply(treat_list_flat, function(x) x %>% arrange(id) %>% slice(u.order))
  
  # for each M, compute risk estimate
  surv_M = lapply(treat_list_flat, function(x) ipwe_km(data$obs_times, data$treat, data$delta, prop_scores, censor, x$treat_rule, tpt))
  
  # format
  surv_M_df = surv_M %>% as.data.frame() %>% t()
  
  # return average across the M
  return(colMeans(surv_M_df))
}


extract_treat <- function(treat_list, data) {
  
  # extract risk estimates (do this for each imputation)
  
  # Inputs:
  # treat_list: treat_rule output from cv_estimator function
  # data: data, sorted by obs time
  # Outputs:
  # vector of treatment rules, in same order as `data`
  
  
  # only works for M=1
  stopifnot('M must be 1' = length(treat_list)==1)
  
  # index to arrange data in original order (i.e. by id) -> sorted obs_times order
  u.order = rank(data$id, ties='first')
  
  # combine treatment recommendations for each M and arrange
  # i.e. get list of df where each df is treat rule for all data
  treat_list_flat = lapply(treat_list, function(x) bind_rows(lapply(x, function(y) as.data.frame(y))))
  treat_list_flat = lapply(treat_list_flat, function(x) x %>% arrange(id) %>% slice(u.order))
  
  
  # treatment rule as named vector
  treat_rule = treat_list_flat[[1]]$treat_rule
  names(treat_rule) = treat_list_flat[[1]]$id
  return(treat_rule)
}


format_coef <- function(eta) {
  
  # formats coefficients from CV treatment rules into data frame
  
  # Inputs:
  # eta: the list of lists coefficients output from jackknife_estimator
  # Outputs:
  # data frame with colnames as coefficient names and K and M variables
  # to represent folds and repeats, respectively
  
  eta_list_df = lapply(eta, function(x) as.data.frame(x) %>% t() %>% as.data.frame() %>% mutate(K = 1:n()))
  n = nrow(eta_list_df[[1]])
  eta_df = bind_rows(eta_list_df) %>% mutate(M = rep(1:length(eta_list_df), n) %>% sort())
  row.names(eta_df) = NULL
  
  return(eta_df)
  
}


get_true_dtr <- function(scenario) {
  
  if (scenario == 0) {
    true_treat_rule = rep(1, nrow(data))
  }
  
  if (scenario == 1) {
    true_treat_rule = as.numeric(-1.9 - 0.2*data$X18 + 0.7*data$X19 < 0)
  }
  
  if (scenario == 2) {
    true_treat_rule = as.numeric(-0.1 + 0.1*data$X16 - 0.1*data$X17 + 0.1*data$X18 - 0.1*data$X19 + 0.3*rowSums(data %>% select(X1:X6, X13:X15)) - 1.2*rowSums(data %>% select(X7:X12)) < 0)
  }
  
  if (scenario == 3) {
    true_treat_rule = as.numeric(- 0.95 + (data$X18-6)^2 + (data$X19-4)^2 < 0)
  }
  
  return(true_treat_rule)
}

get_est_curvs <- function(data, dtr_objs, true_treat_rule, cv_or_whole) {
  
  # prop scores
  prop_frmla = as.formula(paste0('treat ~', paste0('X', 1:19, collapse='+')))
  prop_model = glm(prop_frmla, data = data, family = binomial)
  prop_scores = predict(prop_model, type = 'response')
  
  # censoring weights
  frmla = paste0('Surv(obs_times, ltfu) ~ treat +', paste0('X', 1:19, collapse='+'))
  censor_frmla = as.formula(frmla)
  cens_model = coxph(censor_frmla, data = data)
  cens_weights = survfit(cens_model, newdata = data)$surv %>% t()
  
  # estimated curves
  
  # V_hat(d*), d* = 0,1,d
  est_0 = ipwe_km(data$obs_times, data$treat, data$delta, prop_scores, cens_weights, 0, 1:60)
  est_1 = ipwe_km(data$obs_times, data$treat, data$delta, prop_scores, cens_weights, 1, 1:60)
  est_d = ipwe_km(data$obs_times, data$treat, data$delta, prop_scores, cens_weights, true_treat_rule, 1:60)
  
  if (cv_or_whole=='cv') {
    
    # V_hat(d_hat)
    est_dhat = lapply(dtr_objs, function(x) extract_surv(x$treat_rule, data, prop_scores, cens_weights, 1:60))
    
    # d_hat
    dhat = lapply(dtr_objs, function(x) extract_treat(x$treat_rule, data))
    
  }
  
  if (cv_or_whole=='whole') {
    
    # V_hat(d_hat)
    est_dhat = lapply(dtr_objs, function(x) ipwe_km(data$obs_times, data$treat, data$delta, prop_scores, cens_weights, x$treat_rule, 1:60))
    
    # d_hat
    dhat = lapply(dtr_objs, function(x) x$treat_rule)
    
  }
  
  
  # return results
  return(mget(c('est_0', 'est_1', 'est_d',
                'est_dhat', 'dhat')))
  
}


get_true_curvs <- function(data, scenario, dhat) {
  
  # DGP based on scenario
  if (scenario == 0) {
    
    # oracle
    true_treat_rule = rep(1, nrow(data))
    
    linear.part.0 = -1 - 2*data$X3 + 0.1*data$X18 + 0.5*data$X19
    linear.part.1 = -1.1 - 2*data$X3 + 0.1*data$X18 + 0.5*data$X19
    linear.part.d = pmin(linear.part.0, linear.part.1)
    
  }
  
  if (scenario == 1) {
    
    # oracle
    true_treat_rule = as.numeric(-1.9 - 0.2*data$X18 + 0.7*data$X19 < 0)
    
    linear.part.0 = -1 - 2*data$X3 + 0.1*data$X18 + 0.5*data$X19
    linear.part.1 = -1 - 2*data$X3 - 1.9 - 0.1*data$X18 + 1.2*data$X19
    linear.part.d = pmin(linear.part.0, linear.part.1)
    
  }
  
  if (scenario == 2) {
    
    # oracle
    true_treat_rule = as.numeric(-0.1 + 0.1*data$X16 - 0.1*data$X17 + 0.1*data$X18 - 0.1*data$X19 + 0.3*rowSums(data %>% select(X1:X6, X13:X15)) - 1.2*rowSums(data %>% select(X7:X12)) < 0)
    
    linear.part.0 = -1.5 - 0.2*data$X16 + 0.2*data$X17 - 0.2*data$X18 + 0.5*data$X19 + 1.1*rowSums(data %>% select(X1:X15))
    linear.part.1 = -1.6 - 0.1*data$X16 + 0.1*data$X17 - 0.1*data$X18 + 0.4*data$X19 + 1.4*rowSums(data %>% select(X1:X6, X13:X15)) - 0.1*rowSums(data %>% select(X7:X12))
    linear.part.d = pmin(linear.part.0, linear.part.1)
    
  }
  
  if (scenario == 3) {
    
    # oracle
    true_treat_rule = as.numeric(- 0.95 + (data$X18-6)^2 + (data$X19-4)^2 < 0)
    
    linear.part.0 = -0.4 - 0.3*data$X18 + 0.8*data$X19
    linear.part.1 = -0.4 - 0.3*data$X18 + 0.8*data$X19 - 0.95 + (data$X18-6)^2 + (data$X19-4)^2
    linear.part.d = pmin(linear.part.0, linear.part.1)
    
  }
  
  # common parameters
  lambda = 0.002  # scale
  gamma = 0.02  # shape
  
  # get true survival
  true_0 = GenSurvCurv(linear.part.0, lambda, gamma) %>% apply(2, mean)
  true_1 = GenSurvCurv(linear.part.1, lambda, gamma) %>% apply(2, mean)
  true_d = GenSurvCurv(linear.part.d, lambda, gamma) %>% apply(2, mean)
  
  # V(d_hat)
  linear.part.dhat = lapply(dhat, function(x) ifelse(x, linear.part.1, linear.part.0))
  true_dhat = lapply(linear.part.dhat, function(x) apply(GenSurvCurv(x, lambda, gamma), 2, mean))
  
  # misclassification
  misclass = sapply(dhat, function(x) mean(true_treat_rule != x))
  
  # return results
  return(mget(c('true_0', 'true_1', 'true_d',
                'true_dhat', 'misclass')))
  
}


# convert dictionary
conv.dict <- function(vec, dict) {
  new_vec = dict[vec]
  new_vec[is.na(new_vec)] = vec[is.na(new_vec)]
  return(new_vec)
}


format_results <- function(results, method_name, desc_names, append=NULL) {
  
  # convert to data frame
  results = as.data.frame(results)
  
  # add col names
  colnames(results) = method_name
  
  # pivot to long format and change names using desc_names
  results = results %>%
    mutate(sim_num = 1:n()) %>%
    pivot_longer(cols = -sim_num,
                 names_to = 'method',
                 values_to = 'value') %>%
    mutate(method = conv.dict(method, desc_names))
  
  # append a row
  if (!is.null(append)) {
    results = bind_rows(results,
                        data.frame(sim_num = 1:num_sims,
                                   method = 'Maximum',
                                   value = append)) %>%
      arrange(sim_num)
  }
  
  # return
  return(results)
  
}


combine_results <- function(results_lst) {
  
  # combine all results into single data frame
  results = bind_rows(results_lst, .id = 'scenario')
  
  # format scenario column
  results$scenario = as.character(as.numeric(results$scenario)-1)
  results$scenario = ifelse(results$scenario=='0',
                            'Scenario 0',
                            results$scenario)
  results$scenario = factor(results$scenario,
                            levels = c('Scenario 0', '1', '2', '3'))
  
  # return
  return(results)
  
}


plot_results <- function(results, zero_line=FALSE) {
  
  gg = ggplot(results, aes(x=factor(method, levels=rev(desc_names)), 
                      y=value)) +
    geom_boxplot() +
    coord_flip() +
    labs(x = '', y='') +
  facet_wrap(vars(scenario))
  
  if (zero_line) {
    gg = gg +
      geom_hline(yintercept=0, linetype='dashed', alpha=0.5)
  }
  
  return(gg)
  
}


# Other -------------------------------------------------------------------


do_rist <- function(data, nmin=6, M=50, L=2, tao=60) {
  
  # data: data must be ordered X, delta, times
  # nmin: minimum number of observed data in each node
  # M: number of trees in each fold
  # L: number of folds
  # tao: length of study
  
  # below file from: https://sites.google.com/site/teazrq/software
  source('RISTfunctions.r')
  
  # set other hyperparameters
  P = ncol(rist_data)-2   # number of dimension for X
  K = round(sqrt(P))
  
  # split dataset into not censored
  data_c = data[data$delta==0, ]
  
  # RIST
  R_Muti_ERT_build = Muti_ERT_fit(as.matrix(data[, 1:P]), M, K, L, nmin, SupLogRank=1, tao=tao, impute="random")
  R_Muti_ERT_predict = Muti_ERT_Predict(as.matrix(data_c), R_Muti_ERT_build$Forest_seq[[L]], 
                                        R_Muti_ERT_build$SurvMat_seq[[L]], R_Muti_ERT_build$time_intrest)
  
  # predictions
  times_pred_c = R_Muti_ERT_predict$TMT_predict
  times_pred = data$obs_times
  times_pred[data$delta==0] = times_pred_c
  
  # return observed times with censored times imputed
  return(times_pred)
  
}


Genetic.IPWE <- function(fn, X, u, z, delta, prop_scores, censor, smooth, nvars, standardize,
                         max.generations, pop.size, aug_terms=NULL, init_vals=NULL) {  
  
  # set initial values to 0 if not specified
  if (is.null(init_vals)) (init_vals = rep(0,nvars))
  
  # make sure number of initial values equals number of variables
  stopifnot(length(init_vals) == nvars)
  
  # standardize X?
  if (standardize) {
    X = scale(X)
  }
  
  # augmented version
  if (!is.null(aug_terms)) {
    
    failure.surv.0 = aug_terms$failure.surv.0
    failure.surv.1 = aug_terms$failure.surv.1
    lambda.surv.0 = aug_terms$lambda.surv.0
    lambda.surv.1 = aug_terms$lambda.surv.1
    
    temp <- genoud(fn=fn, X=X, obs_times=u, z=z, delta=delta, prop_scores=prop_scores, censor=censor, 
                   failure.surv.0=failure.surv.0, failure.surv.1=failure.surv.1, 
                   lambda.surv.0=lambda.surv.0, lambda.surv.1=lambda.surv.1, smooth=smooth,
                   nvars=nvars, 
                   Domains=cbind(rep(-5,nvars),rep(5,nvars)),
                   starting.values=init_vals,
                   max=TRUE,
                   pop.size = pop.size,
                   
                   print.level=1,
                   BFGS=FALSE, 
                   optim.method="Nelder-Mead",
                   P9=0, 
                   
                   max.generations=max.generations,
                   hard.generation.limit=TRUE)
  }
  
  # non-augmented version
  else {
    temp <- genoud(fn=fn, X=X, obs_times=u, z=z, delta=delta, prop_scores=prop_scores, censor=censor, smooth=smooth,
                   nvars=nvars, 
                   Domains=cbind(rep(-5,nvars),rep(5,nvars)),
                   starting.values=init_vals,
                   max=TRUE,
                   pop.size = pop.size,
                   
                   print.level=1,
                   BFGS=FALSE, 
                   optim.method="Nelder-Mead",
                   P9=0, 
                   
                   max.generations=max.generations,
                   hard.generation.limit=TRUE)
  }
  
  # extract estimates
  eta.est <- temp$par
  if (sum(eta.est^2)==0) eta.est <- c(1,rep(0,nvars-1))
  valhat.etahat <- temp$value
  
  return(c(eta.est, valhat.etahat))
}


scale_2sd <- function(vec) {
  vec = vec - mean(vec)
  sdx = sd(vec)
  vec = vec / (2*sdx)
  return(vec)
}
