# summarizes simulation scenarios given in `sims_create.R`

source('code/utils_sim.R')

library(survival)
library(haven)
library(readr)
library(ggpubr)
library(tidyverse)

theme_set(theme_bw())


# Covariate stats ---------------------------------------------------------

# i.e. Table 1

n = 2500

data = GenCovs(n, corr=TRUE)

# check if variables match original means
colMeans(data)

# continuous and discrete vars
c_vars = c('X16', 'X17', 'X18', 'X19')
d_vars = paste0('X', 1:15)

# correlations
round(cor(data), 2)

library(ggcorrplot)
corr_c = cor(data[,c_vars])
ggcorrplot(corr_c, method='square', type='upper',
           colors = c('red', 'white', 'blue')) +
  scale_x_discrete(labels = c(expression(X[16]),
                              expression(X[17]),
                              expression(X[18]))) +
  scale_y_discrete(labels = c(expression(X[17]),
                              expression(X[18]),
                              expression(X[19]))) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0),
        axis.text.y = element_text(angle = 0, hjust = 0))
ggsave('results/corr_plot_c.png', plot=last_plot(),
       width = 4, height = 4, units = 'in', dpi = 300)

# png('results/corr_plot_c.png',
#     width = 4, height = 4, units = 'in', res = 300)
# corrplot::corrplot(corr_c, type='upper', diag=FALSE, tl.col='black')
# dev.off()


corr_d = cor(data[,d_vars])
ggcorrplot(corr_d, method='square', type='upper',
           colors = c('red', 'white', 'blue'),
           show.legend = FALSE) +
  scale_x_discrete(labels = c(expression(X[1]),
                              expression(X[2]),
                              expression(X[3]),
                              expression(X[4]),
                              expression(X[5]),
                              expression(X[6]),
                              expression(X[7]),
                              expression(X[8]),
                              expression(X[9]),
                              expression(X[10]),
                              expression(X[11]),
                              expression(X[12]),
                              expression(X[13]),
                              expression(X[14]))) +
  scale_y_discrete(labels = c(expression(X[2]),
                              expression(X[3]),
                              expression(X[4]),
                              expression(X[5]),
                              expression(X[6]),
                              expression(X[7]),
                              expression(X[8]),
                              expression(X[9]),
                              expression(X[10]),
                              expression(X[11]),
                              expression(X[12]),
                              expression(X[13]),
                              expression(X[14]),
                              expression(X[15]))) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0),
        axis.text.y = element_text(angle = 0, hjust = 0))
ggsave('results/corr_plot_d.png', plot=last_plot(),
       width = 4, height = 4, units = 'in', dpi = 300)


# means / medians
data %>% 
  select_at(all_of(c_vars)) %>% 
  mutate(X16 = exp(X16),
         X17 = exp(X17),
         X18 = exp(X18)) %>%
  summarize_all(quantile, probs=c(0.25, 0.5, 0.75)) %>%
  round(1)

data %>%
  select_at(all_of(d_vars)) %>%
  summarize_all(mean) %>%
  mutate_all(function(x) round(x*100, 1))

mean(data$treat)




# Censoring ---------------------------------------------------------------

# censoring parameters
censor.func = function(data) -4 + 0.25*data$X16 + 0.2*data$X18
censor.lambda = 0.01  # scale
censor.gamma = 0.03  # shape


# create covariate data
data = GenCovs(n, corr=TRUE)

# generate censor times from Gompertz distribution
censor.linear.part = censor.func(data)

# generate censoring probabilities
true_cens = GenSurvCurv(censor.linear.part, censor.lambda, censor.gamma)
true_cens_df = as.data.frame(true_cens) %>%
  mutate(id = 1:n()) %>%
  pivot_longer(cols = -id, names_to = 'time', values_to = 'prob',
               names_prefix = 'V') %>%
  mutate(time = as.numeric(time))

# plot
ggplot(true_cens_df, aes(x=time, y=prob, group=id)) +
  geom_line(alpha = 0.02) +
  lims(y = c(0,1)) +
  labs(x = 'Time',
       y = 'Probability')

ggsave('results/cens_plot.png', plot=last_plot(),
       width = 6, height = 4, units = 'in')



# True survival probability -----------------------------------------------

# event parameters
lambda = 0.002  # scale
gamma = 0.02  # shape
coarsen = TRUE


# scenario 0
event.func.0 = function(data) -1 - 2*data$X3 + 0.1*data$X18 + 0.5*data$X19
event.func.1 = function(data) -1.1 - 2*data$X3 + 0.1*data$X18 + 0.5*data$X19
event.func = function(data) ifelse(data$treat==0, event.func.0(data), event.func.1(data))

plot.0 = create_plot(n, censor.func, censor.lambda, censor.gamma,
                     event.func.0, event.func.1, lambda, gamma,
                     corr=TRUE, plot_all_cens = FALSE)


# scenario 1
event.func.0 = function(data) -1 - 2*data$X3 + 0.1*data$X18 + 0.5*data$X19
event.func.1 = function(data) -1 - 2*data$X3 - 1.9 - 0.1*data$X18 + 1.2*data$X19
event.func = function(data) ifelse(data$treat==0, event.func.0(data), event.func.1(data))

plot.1 = create_plot(n, censor.func, censor.lambda, censor.gamma,
                     event.func.0, event.func.1, lambda, gamma,
                     corr=TRUE, plot_all_cens = FALSE)


# scenario 2
event.func.0 = function(data) -1.5 - 0.2*data$X16 + 0.2*data$X17 - 0.2*data$X18 + 0.5*data$X19 + 1.1*rowSums(data %>% select(X1:X15))
event.func.1 = function(data) -1.6 - 0.1*data$X16 + 0.1*data$X17 - 0.1*data$X18 + 0.4*data$X19 + 1.4*rowSums(data %>% select(X1:X6, X13:X15)) - 0.1*rowSums(data %>% select(X7:X12))
event.func = function(data) ifelse(data$treat==0, event.func.0(data), event.func.1(data))

plot.2 = create_plot(n, censor.func, censor.lambda, censor.gamma,
                     event.func.0, event.func.1, lambda, gamma,
                     corr=TRUE, plot_all_cens = FALSE)


# scenario 3
event.func.0 = function(data) -0.4 - 0.3*data$X18 + 0.8*data$X19
event.func.1 = function(data) -0.4 - 0.3*data$X18 + 0.8*data$X19 - 0.95 + (data$X18-6)^2 + (data$X19-4)^2
event.func = function(data) ifelse(data$treat==0, event.func.0(data), event.func.1(data))

plot.3 = create_plot(n, censor.func, censor.lambda, censor.gamma,
                     event.func.0, event.func.1, lambda, gamma,
                     corr=TRUE, plot_all_cens = FALSE)


# combine plots
plot.data = bind_rows(plot.0, plot.1, plot.2, plot.3, .id='scenario')
plot.data$scenario = as.character(as.numeric(plot.data$scenario)-1)
plot.data$scenario = ifelse(plot.data$scenario=='0',
                            'Scenario 0',
                            plot.data$scenario) %>%
  factor(levels = c('Scenario 0', '1', '2', '3'), ordered = TRUE)

ggplot(plot.data, aes(x=time, y=surv_prob, color=dtr, linetype=dtr)) +
  geom_line() +
  lims(y = c(0,1)) +
  facet_wrap(vars(scenario)) +
  labs(x = 'Time',
       y = 'Survival probability',
       color = '',
       linetype = '') +
  scale_color_manual(labels = c('0', '1', 'Optimal'),
                     values = c('grey75', 'grey25', 'black')) +
  scale_linetype_manual(labels = c('0', '1', 'Optimal'),
                        values = c('dashed', 'dashed', 'solid')) +
  theme(legend.position = 'bottom')

# save plot
ggsave('results/true_surv.tiff', plot=last_plot(),
       width = 5, height = 6, dpi = 800,
       units = 'in')



# Visualize decision boundary ---------------------------------------------


# censoring parameters
censor.func = function(data) -4 + 0.25*data$X16 + 0.2*data$X18
censor.lambda = 0.01  # scale
censor.gamma = 0.03  # shape

coarsen = TRUE

# scenario 1
event.func.0 = function(data) -1 - 2*data$X3 + 0.1*data$X18 + 0.5*data$X19
event.func.1 = function(data) -1 - 2*data$X3 - 1.9 - 0.1*data$X18 + 1.2*data$X19
event.func = function(data) ifelse(data$treat==0, event.func.0(data), event.func.1(data))

data.1 = GenSims(n, censor.func, censor.lambda, censor.gamma,
                 event.func, lambda, gamma,
                 end=60, coarsen=coarsen)

gg_1 = ggplot() +
  geom_point(data=data.1, aes(x=X18, y=X19), alpha=0.1, shape=1, inherit.aes=FALSE) +
  geom_abline(intercept = 1.9/0.7, slope = 0.2/0.7) +
  scale_x_continuous(limits = c(0, 10)) +
  scale_y_continuous(limits = c(0, 10)) +
  labs(x = expression(X[18]), y = expression(X[19])) +
  theme(text = element_text(size = 15)) +
  ggtitle('Scenario 1')

# ggsave('results/treat_boundary_1.png', 
#        width = 4, height = 4, units = 'in')


# scenario 3
event.func.0 = function(data) -0.4 - 0.3*data$X18 + 0.8*data$X19
event.func.1 = function(data) -0.4 - 0.3*data$X18 + 0.8*data$X19 - 0.95 + (data$X18-6)^2 + (data$X19-4)^2
event.func = function(data) ifelse(data$treat==0, event.func.0(data), event.func.1(data))

data.3 = GenSims(n, censor.func, censor.lambda, censor.gamma,
                 event.func, lambda, gamma,
                 end=60, coarsen=coarsen)

grid_vals = seq(0, 10, length.out=100)
treat_grid = expand.grid(x=grid_vals, y=grid_vals)
treat_grid$z = -0.95 + (treat_grid$x-6)^2 + (treat_grid$y-4)^2
gg_2 = ggplot(treat_grid, aes(x=x, y=y, z=z)) +
  geom_contour(breaks = 0, color='black') +
  geom_point(data=data.3, aes(x=X18, y=X19), alpha=0.1, shape=1, inherit.aes=FALSE) +
  # geom_density_2d(data=data, aes(x=X18, y=X19), color='red', inherit.aes=FALSE) +
  scale_x_continuous(limits = c(0, 10)) +
  scale_y_continuous(limits = c(0, 10)) +
  labs(x = expression(X[18]), y = expression(X[19])) +
  theme(text = element_text(size = 15)) +
  ggtitle('Scenario 3')

# ggsave('results/treat_boundary_3.png', 
#        width = 4, height = 4, units = 'in')

ggarrange(gg_1, gg_2, nrow = 1, ncol = 2)

ggsave('results/treat_boundary.tiff', 
       width = 9, height = 4, dpi = 800, units = 'in')