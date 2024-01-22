# estimators.R
# Compares different estimators for coarsened and un-coarsened data
# Goal is unbiasedness


# Set-up ------------------------------------------------------------------

source('code/utils_sim.R')

library(survival)
library(tidyverse)

# set ggplot theme
theme_set(theme_bw())



# Scenario ----------------------------------------------------------------

# scenario
scenario = 2

# number of sims
num_sims = 1000

# sample size
n = 2500

# coarsen
coarsen = TRUE

# censor function
censor.func = function(data) -4 + 0.25*data$X16 + 0.2*data$X18
censor.lambda = 0.01  # scale
censor.gamma = 0.03  # shape

# generate survival times
lambda = 0.002  # scale
gamma = 0.02  # shape

# linear part under different treatment rules
event.func.0 = function(data) -1.5 - 0.2*data$X16 + 0.2*data$X17 - 0.2*data$X18 + 0.5*data$X19 + 1.1*rowSums(data %>% select(X1:X15))
event.func.1 = function(data) -1.6 - 0.1*data$X16 + 0.1*data$X17 - 0.1*data$X18 + 0.4*data$X19 + 1.4*rowSums(data %>% select(X1:X6, X13:X15)) - 0.1*rowSums(data %>% select(X7:X12))
event.func = function(data) ifelse(data$treat==0, event.func.0(data), event.func.1(data))

# true treatment rule function
true.treat.func = function(data) as.numeric(-0.1 + 0.1*data$X16 - 0.1*data$X17 + 0.1*data$X18 - 0.1*data$X19 + 0.3*rowSums(data %>% select(X1:X6, X13:X15)) - 1.2*rowSums(data %>% select(X7:X12)) < 0)

# create sims
set.seed(7894117)
simlist = lapply(1:num_sims, function(b) GenSims(n, censor.func, censor.lambda, censor.gamma,
                                                 event.func, lambda, gamma,
                                                 end=61, coarsen=coarsen, corr=TRUE))


# Simulate ----------------------------------------------------------------


# initialize matrices to store results
km0 = km1 = kmd = ht0 = ht1 = htd = hj0 = hj1 = hjd =
  curv0_truth = curv0_truth = curv1_truth = curvd_truth = 
  matrix(NA, nrow=num_sims, ncol=60)

# run for each sim
for (sim_num in 1:num_sims) {
  
  if (sim_num==1) (print(Sys.time()))
  
  data = simlist[[sim_num]]
  
  # prop scores
  prop_model = glm(treat ~ ., data=data %>% select(treat, starts_with('X')),
                   family = binomial(link='logit'))
  prop_scores = prop_model$fitted.values
  
  # censoring weights
  censor_model = coxph(Surv(obs_times, ltfu) ~ ., data=data %>% select(obs_times, ltfu, starts_with('X')))
  tmp = survfit(censor_model, newdata = data)
  cens_weights = summary(tmp, times=1:61, extend=TRUE)$surv %>% t()
  
  which_u = as.integer(factor(data$obs_times))
  cens_weights2 = sapply(1:n, function(i) cens_weights[i, which_u[i]])
  
  
  censor_model = coxph(Surv(obs_times, ltfu) ~ ., data=data %>% select(obs_times, ltfu, starts_with('X')))
  tmp = survfit(censor_model, newdata = data)
  cens_weights_km = summary(tmp, times=1:60, extend=TRUE)$surv %>% t()
  
  # true treatment rule
  true_treat_rule = true.treat.func(data)
  
  # estimated curves ipw-km
  km0[sim_num,] = ipwe_km(data$obs_times, data$treat, data$delta, prop_scores, cens_weights_km, 0, 1:60)
  km1[sim_num,] = ipwe_km(data$obs_times, data$treat, data$delta, prop_scores, cens_weights_km, 1, 1:60)
  kmd[sim_num,] = ipwe_km(data$obs_times, data$treat, data$delta, prop_scores, cens_weights_km, true_treat_rule, 1:60)
  
  # estimated curves ipw-ht
  ht0[sim_num,] = sapply(1:60, function(tpt) ipw_estimator(data$obs_times, data$treat, 1-data$ltfu, prop_scores, cens_weights2, 0, tpt))
  ht1[sim_num,] = sapply(1:60, function(tpt) ipw_estimator(data$obs_times, data$treat, 1-data$ltfu, prop_scores, cens_weights2, 1, tpt))
  htd[sim_num,] = sapply(1:60, function(tpt) ipw_estimator(data$obs_times, data$treat, 1-data$ltfu, prop_scores, cens_weights2, true_treat_rule, tpt))
  
  # estimated curves ipw-hj
  hj0[sim_num,] = sapply(1:60, function(tpt) ipw_estimator(data$obs_times, data$treat, 1-data$ltfu, prop_scores, cens_weights2, 0, tpt, hajek=TRUE))
  hj1[sim_num,] = sapply(1:60, function(tpt) ipw_estimator(data$obs_times, data$treat, 1-data$ltfu, prop_scores, cens_weights2, 1, tpt, hajek=TRUE))
  hjd[sim_num,] = sapply(1:60, function(tpt) ipw_estimator(data$obs_times, data$treat, 1-data$ltfu, prop_scores, cens_weights2, true_treat_rule, tpt, hajek=TRUE))
  
  
  # true curves
  linear.part.0 = event.func.0(data)
  linear.part.1 = event.func.1(data)
  linear.part.d = pmin(linear.part.0, linear.part.1)
  
  true_surv_0 = GenSurvCurv(linear.part.0, lambda, gamma, 1:60)
  true_surv_1 = GenSurvCurv(linear.part.1, lambda, gamma, 1:60)
  true_surv_d = GenSurvCurv(linear.part.d, lambda, gamma, 1:60)
  
  curv0_truth[sim_num,] = apply(true_surv_0, 2, mean)
  curv1_truth[sim_num,] = apply(true_surv_1, 2, mean)
  curvd_truth[sim_num,] = apply(true_surv_d, 2, mean)
  
  # check progress
  if (sim_num %% 10 == 0) {
    print(sim_num)
    print(Sys.time())
  }
  
}



# Evaluate ----------------------------------------------------------------


# boxplots

# treat 0
boxplot(km0 - curv0_truth, main='treat 0, km', ylim=c(-0.02,0.02))
abline(h=0, col='red')

boxplot(ht0 - curv0_truth, main='treat 0, ht', ylim=c(-0.02,0.02))
abline(h=0, col='red')

boxplot(hj0 - curv0_truth, main='treat 0, hajek', ylim=c(-0.02,0.02))
abline(h=0, col='red')

# treat 1
boxplot(km1 - curv1_truth, main='treat 1, km', ylim=c(-0.02,0.02))
abline(h=0, col='red')

boxplot(ht1 - curv1_truth, main='treat 1, ht', ylim=c(-0.02,0.02))
abline(h=0, col='red')

boxplot(hj1 - curv1_truth, main='treat 1, hajek', ylim=c(-0.02,0.02))
abline(h=0, col='red')

# treat d
boxplot(kmd - curvd_truth, main='treat dtr, km', ylim=c(-0.02,0.02))
abline(h=0, col='red')

boxplot(htd - curvd_truth, main='treat dtr, ht', ylim=c(-0.02,0.02))
abline(h=0, col='red')

boxplot(hjd - curvd_truth, main='treat dtr, hajek', ylim=c(-0.02,0.02))
abline(h=0, col='red')


# make nice plot
v1 = data.frame(kmd-curvd_truth)
v2 = data.frame(htd-curvd_truth)
v3 = data.frame(hjd-curvd_truth)
plot.data = bind_rows(v1, v2, v3, .id='estimator')
plot.data = plot.data %>%
  pivot_longer(cols = !estimator,
               names_prefix = 'X',
               names_to = 'time',
               values_to = 'diff')
plot.data = plot.data %>%
  mutate(time = factor(time, levels=1:60, ordered=TRUE))

ggplot(plot.data, aes(x=time, y=diff)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0) +
  facet_grid(rows = vars(estimator)) +
  labs(x = 'Time',
       y = '') +
  lims(y = c(-0.03, 0.03)) +
  scale_x_discrete(breaks = c(1, seq(5,60,by=5)),
                   labels = c(1, seq(5,60,by=5)))


ggplot(plot.data %>% filter(time==60), aes(x=estimator, y=diff)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0) +
  labs(x = '',
       y = '') +
  scale_x_discrete(labels = c('Kaplan-Meier', 'Horvitz-Thompson', 'Hajek'))
ggsave('results/unbias_d.png', plot=last_plot(),
       width = 5, height = 4, units = 'in')


colMeans(km0 - curv0_truth) %>% round(5)
colMeans(km1 - curv1_truth) %>% round(5)
colMeans(kmd - curvd_truth) %>% round(5)

apply(km0 - curv0_truth, 2, median) %>% round(5)
apply(km1 - curv1_truth, 2, median) %>% round(5)
apply(kmd - curvd_truth, 2, median) %>% round(5)


# statistics to report
plot.data %>% 
  filter(time==60) %>% 
  group_by(estimator) %>% 
  summarize(mean = mean(diff) %>% round(3),
            sd = sd(diff) %>% round(3))
