devtools::install_github("emmekeaarts/mHMMbayes@one-state")
library(mHMMbayes)
library(ggplot2)
library(reshape2)

# Setting up ####
# Simulation settings 
settings4 <- read.csv("sim_scenarios_V4.csv")
true_m  <- 3
n_sim <- 100
J <- 2000
burn_in <- 500

# Reading in tables
Label_switch_proxy_3st <- readRDS("result_tables/Label_switch_proxy_3st.RDS")
Performance_emission_3st <- readRDS(file = "result_tables/Performance_emission_3st.RDS")
Model_selection_3st <- readRDS(file = "result_tables/Model_selection_3st.RDS")
Performance_gamma_3st <- readRDS(file = "result_tables/Performance_gamma_3st.RDS")



# MODEL SELECTION ####
## model selection, corrected for label switching ####
### Combining model selection table with label switching ####

# taking mean RMSE over each simulation run such that we combine RMSE over dependent variables 
aggr_Label_switch_proxy_3st <- aggregate(Label_switch_proxy_3st, by = list(Label_switch_proxy_3st$sim_iteration, 
                                                                           Label_switch_proxy_3st$KL_div,
                                                                           Label_switch_proxy_3st$n_dep, 
                                                                           Label_switch_proxy_3st$n_t, 
                                                                           Label_switch_proxy_3st$n), FUN = mean)


Model_selection_3st_AIC <- reshape(Model_selection_3st, direction = "wide", timevar = "K_candidate", 
               idvar = c("sim_iteration", "true_m", "n", "n_t", "n_dep", "KL_div"), drop = "AICc")

Model_selection_3st_AIC$selected <- apply(Model_selection_3st_AIC[,7:10], 1, which.min)
Model_selection_3st_AIC$correct_selected <- 0
Model_selection_3st_AIC$correct_selected[Model_selection_3st_AIC$selected == true_m] <- 1  





# Making sure ordering is same in performance table and label switching table, saving in new object 
noL_Model_selection_3st_AIC <- Model_selection_3st_AIC[order(
  Model_selection_3st_AIC$KL_div, 
  Model_selection_3st_AIC$n_dep, 
  Model_selection_3st_AIC$n_t, 
  Model_selection_3st_AIC$n,
  Model_selection_3st_AIC$sim_iteration
),]


# Bit double, but to make sure ordering is the same compared to performance table 
order_aggr_Label_switch_proxy_3st <- aggr_Label_switch_proxy_3st[order(aggr_Label_switch_proxy_3st$KL_div, 
                                                                       aggr_Label_switch_proxy_3st$n_dep, 
                                                                       aggr_Label_switch_proxy_3st$n_t, 
                                                                       aggr_Label_switch_proxy_3st$n,
                                                                       aggr_Label_switch_proxy_3st$sim_iteration
),]


# adding the RMSE to the performance table, RMSE needs to be repeated M * n_dep timees, where n_dep varies over the scenarios.. 
noL_Model_selection_3st_AIC$RMSE <- order_aggr_Label_switch_proxy_3st$RMSE

# removing 1) simulation runs that have RMSE < 0.20 AND 2) scenarios that pertain to KLD = 3 
noL_Model_selection_3st_AIC <- noL_Model_selection_3st_AIC[noL_Model_selection_3st_AIC$RMSE > 0.20,]


### Aggregating to MEAN within one simulation setting (= percentage of iterations resulting in correct selection)
aggr_noL_Model_selection_3st_AIC <- aggregate(noL_Model_selection_3st_AIC, by = list(noL_Model_selection_3st_AIC$KL_div, 
                                                                                     noL_Model_selection_3st_AIC$n_dep, 
                                                                                     noL_Model_selection_3st_AIC$n_t, 
                                                                                     noL_Model_selection_3st_AIC$n), FUN = mean)

aggr_noL_Model_selection_3st_AIC <- aggr_noL_Model_selection_3st_AIC[,c(6:10,16)]

# removing the scenarios that had less than 20% of the simulation runs in 
aggr_noL_Model_selection_3st_AIC <- aggr_noL_Model_selection_3st_AIC[-which((aggr_noL_Model_selection_3st_AIC$KL_div == 5 &
                                                                               aggr_noL_Model_selection_3st_AIC$n_dep == 2 &
                                                                               aggr_noL_Model_selection_3st_AIC$n == 120 &
                                                                               aggr_noL_Model_selection_3st_AIC$n_t > 100) | 
                                                                                          (aggr_noL_Model_selection_3st_AIC$KL_div == 5 &
                                                                                             aggr_noL_Model_selection_3st_AIC$n_dep == 2 &
                                                                                             aggr_noL_Model_selection_3st_AIC$n == 60 &
                                                                                             aggr_noL_Model_selection_3st_AIC$n_t > 100)),]

mean(aggr_noL_Model_selection_3st_AIC$correct_selected[aggr_noL_Model_selection_3st_AIC$KL_div > 5])
mean(aggr_noL_Model_selection_3st_AIC$correct_selected[aggr_noL_Model_selection_3st_AIC$KL_div == 5 & aggr_noL_Model_selection_3st_AIC$n_t > 50])
mean(aggr_noL_Model_selection_3st_AIC$correct_selected[aggr_noL_Model_selection_3st_AIC$KL_div == 5 & aggr_noL_Model_selection_3st_AIC$n_t == 50])

mean(aggr_bias_noL_Performance_emission_3st$abs_rel_bias[aggr_bias_noL_Performance_emission_3st$KL_div == 5])
mean(aggr_bias_noL_Performance_emission_3st$abs_rel_bias[aggr_bias_noL_Performance_emission_3st$KL_div == 7])

# BIAS ####
## Bias in emission means, no label switching correction ####

### Calculating absolute relative bias 
Performance_emission_3st$abs_rel_bias <- abs(Performance_emission_3st$mean_hat - Performance_emission_3st$mean_true) / abs(Performance_emission_3st$mean_true)
Performance_emission_3st$SD_rel_bias <- (Performance_emission_3st$SD_hat - Performance_emission_3st$SD_true) / abs(Performance_emission_3st$SD_true)

### Aggregating to MEDIAN within one simulation setting 
aggr_Performance_emission_3st <- aggregate(Performance_emission_3st, by = list(Performance_emission_3st$KL_div, 
                                                                               Performance_emission_3st$n_dep, 
                                                                               Performance_emission_3st$n_t, 
                                                                               Performance_emission_3st$n), FUN = median)

ggplot(data = aggr_Performance_emission_3st, mapping = aes(x = n_t, y = abs_rel_bias, group = as.factor(n), color = as.factor(n))) +
  geom_point() + 
  geom_line() +
  ylab("Bias emission means") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(n_dep)), cols = vars(as.factor(KL_div))) +
  xlab("number of observations per subject") +
  geom_hline(yintercept=c(.10), linetype="dashed", color = "grey") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))


## Bias in emission means, corrected for label switching ####
### Combining emission performance table with label switching ####

# taking mean RMSE over each simulation run such that we combine RMSE over dependent variables 
aggr_Label_switch_proxy_3st <- aggregate(Label_switch_proxy_3st, by = list(Label_switch_proxy_3st$sim_iteration, 
                                                                           Label_switch_proxy_3st$KL_div,
                                                                           Label_switch_proxy_3st$n_dep, 
                                                                           Label_switch_proxy_3st$n_t, 
                                                                           Label_switch_proxy_3st$n), FUN = mean)

# Making sure ordering is same in performance table and label switching table, saving in new object 
bias_noL_Performance_emission_3st <- Performance_emission_3st[order(
  Performance_emission_3st$KL_div, 
  Performance_emission_3st$n_dep, 
  Performance_emission_3st$n_t, 
  Performance_emission_3st$n,
  Performance_emission_3st$sim_iteration,
  Performance_emission_3st$dep,
  Performance_emission_3st$k
  
),]

# Bit double, but to make sure ordering is the same compared to performance table 
order_aggr_Label_switch_proxy_3st <- aggr_Label_switch_proxy_3st[order(aggr_Label_switch_proxy_3st$KL_div, 
                                                                       aggr_Label_switch_proxy_3st$n_dep, 
                                                                       aggr_Label_switch_proxy_3st$n_t, 
                                                                       aggr_Label_switch_proxy_3st$n,
                                                                       aggr_Label_switch_proxy_3st$sim_iteration
),]

# adding the RMSE to the performance table, RMSE needs to be repeated M * n_dep timees, where n_dep varies over the scenarios.. 
bias_noL_Performance_emission_3st$RMSE <- rep(order_aggr_Label_switch_proxy_3st$RMSE, 
                                              times = order_aggr_Label_switch_proxy_3st$n_dep * true_m)

# removing 1) simulation runs that have RMSE < 0.20 AND 2) scenarios that pertain to KLD = 3 
bias_noL_Performance_emission_3st <- bias_noL_Performance_emission_3st[bias_noL_Performance_emission_3st$RMSE > 0.20 & 
                                                                         bias_noL_Performance_emission_3st$KL_div > 3,]

### Inspecting absolute relative bias for Gaussian emission means ####
# to MEDIAN within one simulation setting 
aggr_bias_noL_Performance_emission_3st <- aggregate(bias_noL_Performance_emission_3st, by = list(bias_noL_Performance_emission_3st$KL_div, 
                                                                                                 bias_noL_Performance_emission_3st$n_dep, 
                                                                                                 bias_noL_Performance_emission_3st$n_t, 
                                                                                                 bias_noL_Performance_emission_3st$n), FUN = median)

# removing the scenarios that had less than 20% of the simulation runs in 
aggr_bias_noL_Performance_emission_3st <- aggr_bias_noL_Performance_emission_3st[-which((aggr_bias_noL_Performance_emission_3st$KL_div == 5 &
                                                                                          aggr_bias_noL_Performance_emission_3st$n_dep == 2 &
                                                                                          aggr_bias_noL_Performance_emission_3st$n == 120 &
                                                                                          aggr_bias_noL_Performance_emission_3st$n_t > 100) | 
                                                                                          (aggr_bias_noL_Performance_emission_3st$KL_div == 5 &
                                                                                              aggr_bias_noL_Performance_emission_3st$n_dep == 2 &
                                                                                              aggr_bias_noL_Performance_emission_3st$n == 60 &
                                                                                              aggr_bias_noL_Performance_emission_3st$n_t > 100)),]

# inspecting result
# comparing the figures, I have the feeling that scenarios affected more by label switching show larger differences. 
ggplot(data = aggr_bias_noL_Performance_emission_3st, mapping = aes(x = n_t, y = abs_rel_bias, group = as.factor(n), color = as.factor(n))) +
  geom_point() + 
  geom_line() +
  ylab("Bias emission means") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  # scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep))) +
  xlab("number of observations per subject") +
  geom_hline(yintercept=c(.10), linetype="dashed", color = "grey") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))


mean(aggr_bias_noL_Performance_emission_3st$abs_rel_bias)
mean(aggr_bias_noL_Performance_emission_3st$abs_rel_bias[aggr_bias_noL_Performance_emission_3st$KL_div == 5])
mean(aggr_bias_noL_Performance_emission_3st$abs_rel_bias[aggr_bias_noL_Performance_emission_3st$KL_div == 7])

mean(aggr_bias_noL_Performance_emission_3st$abs_rel_bias[aggr_bias_noL_Performance_emission_3st$n_dep == 2])
mean(aggr_bias_noL_Performance_emission_3st$abs_rel_bias[aggr_bias_noL_Performance_emission_3st$n_dep == 4])
mean(aggr_bias_noL_Performance_emission_3st$abs_rel_bias[aggr_bias_noL_Performance_emission_3st$n_dep == 8])

mean(aggr_bias_noL_Performance_emission_3st$abs_rel_bias[aggr_bias_noL_Performance_emission_3st$n == 15])
mean(aggr_bias_noL_Performance_emission_3st$abs_rel_bias[aggr_bias_noL_Performance_emission_3st$n == 30])
mean(aggr_bias_noL_Performance_emission_3st$abs_rel_bias[aggr_bias_noL_Performance_emission_3st$n == 60])
mean(aggr_bias_noL_Performance_emission_3st$abs_rel_bias[aggr_bias_noL_Performance_emission_3st$n == 120])



#### inspecting relative bias for gaussian emission SDs ####
mean(aggr_bias_noL_Performance_emission_3st$SD_rel_bias)
mean(aggr_bias_noL_Performance_emission_3st$SD_rel_bias[aggr_bias_noL_Performance_emission_3st$KL_div == 5])
mean(aggr_bias_noL_Performance_emission_3st$SD_rel_bias[aggr_bias_noL_Performance_emission_3st$KL_div == 7])

mean(aggr_bias_noL_Performance_emission_3st$SD_rel_bias[aggr_bias_noL_Performance_emission_3st$n_dep == 2])
mean(aggr_bias_noL_Performance_emission_3st$SD_rel_bias[aggr_bias_noL_Performance_emission_3st$n_dep == 4])
mean(aggr_bias_noL_Performance_emission_3st$SD_rel_bias[aggr_bias_noL_Performance_emission_3st$n_dep == 8])

mean(aggr_bias_noL_Performance_emission_3st$SD_rel_bias[aggr_bias_noL_Performance_emission_3st$n == 15])
mean(aggr_bias_noL_Performance_emission_3st$SD_rel_bias[aggr_bias_noL_Performance_emission_3st$n == 30])
mean(aggr_bias_noL_Performance_emission_3st$SD_rel_bias[aggr_bias_noL_Performance_emission_3st$n == 60])
mean(aggr_bias_noL_Performance_emission_3st$SD_rel_bias[aggr_bias_noL_Performance_emission_3st$n == 120])

mean(aggr_bias_noL_Performance_emission_3st$SD_rel_bias[aggr_bias_noL_Performance_emission_3st$n_t == 50])
mean(aggr_bias_noL_Performance_emission_3st$SD_rel_bias[aggr_bias_noL_Performance_emission_3st$n_t == 100])
mean(aggr_bias_noL_Performance_emission_3st$SD_rel_bias[aggr_bias_noL_Performance_emission_3st$n_t == 200])
mean(aggr_bias_noL_Performance_emission_3st$SD_rel_bias[aggr_bias_noL_Performance_emission_3st$n_t == 400])
mean(aggr_bias_noL_Performance_emission_3st$SD_rel_bias[aggr_bias_noL_Performance_emission_3st$n_t == 800])

mean(aggr_bias_noL_Performance_emission_3st$SD_rel_bias[aggr_bias_noL_Performance_emission_3st$KL_div == 5 & 
                                                          aggr_bias_noL_Performance_emission_3st$n_dep == 2 &
                                                          aggr_bias_noL_Performance_emission_3st$n != 15])

ggplot(data = aggr_bias_noL_Performance_emission_3st, mapping = aes(x = n_t, y = SD_rel_bias, group = as.factor(n), color = as.factor(n))) +
  geom_point() + 
  geom_line() +
  ylab("Bias emission SD") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  # scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep))) +
  xlab("number of observations per subject") +
  geom_hline(yintercept=c(.10), linetype="dashed", color = "grey") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

## Bias in transition probabilities, corrected for label switching ####
# First, calculate bias
Performance_gamma_3st$mean_bias <- (Performance_gamma_3st$gamma_ij_hat - Performance_gamma_3st$gamma_ij_true) 
# add indicator for self-transitions
Performance_gamma_3st$Diag <- Performance_gamma_3st$from_state_i == Performance_gamma_3st$to_state_j * 1

### Combining transition performance table with label switching ####
# making sure ordering is the same, save in new object 
bias_noL_Performance_gamma_3st <- Performance_gamma_3st[order(
  Performance_gamma_3st$KL_div, 
  Performance_gamma_3st$n_dep, 
  Performance_gamma_3st$n_t, 
  Performance_gamma_3st$n,
  Performance_gamma_3st$sim_iteration,
  Performance_gamma_3st$from_state_i,
  Performance_gamma_3st$to_state_j
  
),]

# adding the RMSE to the performance table, RMSE needs to be repeated M * M times 
bias_noL_Performance_gamma_3st$RMSE <- rep(order_aggr_Label_switch_proxy_3st$RMSE, 
                                           each = true_m * true_m)

# removing 1) simulation runs that have RMSE < 0.20 AND 2) scenarios that pertain to KLD = 3 
bias_noL_Performance_gamma_3st <- bias_noL_Performance_gamma_3st[bias_noL_Performance_gamma_3st$RMSE > 0.20 & 
                                                                   bias_noL_Performance_gamma_3st$KL_div > 3,]


### Inspecting mean bias for transition probs ####
# to MEAN within one simulation setting 
bias_noL_summary <- aggregate(bias_noL_Performance_gamma_3st, by = list(bias_noL_Performance_gamma_3st$Diag, 
                                                                        bias_noL_Performance_gamma_3st$n_dep, 
                                                                        bias_noL_Performance_gamma_3st$KL_div, 
                                                                        bias_noL_Performance_gamma_3st$n_t,
                                                                        bias_noL_Performance_gamma_3st$n), FUN = mean)

# removing the scenarios that had less than 10% of the simulation runs in 
bias_noL_summary <- bias_noL_summary[-which((bias_noL_summary$KL_div == 5 &
                                              bias_noL_summary$n_dep == 2 &
                                              bias_noL_summary$n == 120 &
                                              bias_noL_summary$n_t > 100) |
                                              (bias_noL_summary$KL_div == 5 &
                                                 bias_noL_summary$n_dep == 2 &
                                                 bias_noL_summary$n == 60 &
                                                 bias_noL_summary$n_t > 100)),]

bias_noL_summary$Diag <- factor(bias_noL_summary$Diag)

mean(bias_noL_summary$mean_bias[bias_noL_summary$Diag == 1])
mean(bias_noL_summary$mean_bias[bias_noL_summary$Diag == 0])
mean(bias_noL_summary$mean_bias[bias_noL_summary$Diag == 1 & bias_noL_summary$KL_div == 5 & bias_noL_summary$n_dep == 2])
mean(bias_noL_summary$mean_bias[bias_noL_summary$Diag == 0 & bias_noL_summary$KL_div == 5 & bias_noL_summary$n_dep == 2])

mean(bias_noL_summary$mean_bias[bias_noL_summary$Diag == 1 & bias_noL_summary$KL_div != 5 & bias_noL_summary$n_dep != 2])
mean(bias_noL_summary$mean_bias[bias_noL_summary$Diag == 0 & bias_noL_summary$KL_div != 5 & bias_noL_summary$n_dep != 2])

mean(bias_noL_summary$mean_bias[bias_noL_summary$Diag == 1 & bias_noL_summary$KL_div == 7])
mean(bias_noL_summary$mean_bias[bias_noL_summary$Diag == 0 & bias_noL_summary$KL_div == 7])

# Inspecting results, I corrected a mistake with the labeling of the self-transitions in the plot 

# Here, I out commented the log2 scale for ease of comparison, 
# I think the largest differences are again in scenarios most affected by label switching, that is KLD 5 and p = 2 (most), 
# and KLD = 5 and p = 4 (smaller differences)

ggplot(data = bias_noL_summary, mapping = aes(x = n_t, y = mean_bias, group = interaction(as.factor(n), Diag), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  ylab("median relative bias") +
  geom_line(aes(linetype = Diag)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep))) +
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_linetype_manual(name = "Transition type", 
                        values = c(2, 1),
                        labels = c("Off-diagonal", "Diagonal")) +
  # scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  geom_hline(yintercept=0, linetype="solid", color = "grey") +
  geom_hline(yintercept=c(-0.04, 0.04), linetype="dashed", color = "grey") +
  xlab("number of observations per subject") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))




# COVERAGE ####
## Emsission distribution ####
### creating object that contains RMSE values and coverage ####

# Here, I reuse the object created under inspecting bias for the emission distribution, that allready includes RMSE in the correct format
cov_noL_Performance_emission_3st <- bias_noL_Performance_emission_3st

# Coverage yes/no for means
cov_noL_Performance_emission_3st$cov_mean <- (cov_noL_Performance_emission_3st$mean_true > cov_noL_Performance_emission_3st$mean_95_lower &
                                                cov_noL_Performance_emission_3st$mean_true < cov_noL_Performance_emission_3st$mean_95_upper) * 1 
# Coverage yes/no for SDs
cov_noL_Performance_emission_3st$cov_SD <- (cov_noL_Performance_emission_3st$SD_true > cov_noL_Performance_emission_3st$SD_95_lower &
                                              cov_noL_Performance_emission_3st$SD_true < cov_noL_Performance_emission_3st$SD_95_upper) * 1 

# aggregate to MEAN (= proportion) within one simulation setting 
aggr_cov_noL_Performance_emission_3st <- aggregate(cov_noL_Performance_emission_3st, by = list(cov_noL_Performance_emission_3st$KL_div, 
                                                                                               cov_noL_Performance_emission_3st$n_dep, 
                                                                                               cov_noL_Performance_emission_3st$n_t, 
                                                                                               cov_noL_Performance_emission_3st$n), FUN = mean)

# removing the scenarios that had less than 10% of the simulation runs in 
aggr_cov_noL_Performance_emission_3st <- aggr_cov_noL_Performance_emission_3st[-which((aggr_cov_noL_Performance_emission_3st$KL_div == 5 &
                                                                                        aggr_cov_noL_Performance_emission_3st$n_dep == 2 &
                                                                                        aggr_cov_noL_Performance_emission_3st$n == 120 &
                                                                                        aggr_cov_noL_Performance_emission_3st$n_t > 100) |
                                                                                        (aggr_cov_noL_Performance_emission_3st$KL_div == 5 &
                                                                                           aggr_cov_noL_Performance_emission_3st$n_dep == 2 &
                                                                                           aggr_cov_noL_Performance_emission_3st$n == 60 &
                                                                                           aggr_cov_noL_Performance_emission_3st$n_t > 100) ),]
mean(aggr_cov_noL_Performance_emission_3st$cov_mean)
mean(aggr_cov_noL_Performance_emission_3st$cov_mean[aggr_cov_noL_Performance_emission_3st$KL_div == 5])
mean(aggr_cov_noL_Performance_emission_3st$cov_mean[aggr_bias_noL_Performance_emission_3st$KL_div == 7])

mean(aggr_cov_noL_Performance_emission_3st$cov_mean[aggr_cov_noL_Performance_emission_3st$n_dep == 2 &
                                                      aggr_cov_noL_Performance_emission_3st$KL_div == 5])
mean(aggr_cov_noL_Performance_emission_3st$cov_mean[aggr_cov_noL_Performance_emission_3st$n_dep == 4 &
                                                      aggr_cov_noL_Performance_emission_3st$KL_div == 5])
mean(aggr_cov_noL_Performance_emission_3st$cov_mean[aggr_cov_noL_Performance_emission_3st$n_dep == 8 &
                                                      aggr_cov_noL_Performance_emission_3st$KL_div == 5])

mean(aggr_cov_noL_Performance_emission_3st$cov_mean[aggr_cov_noL_Performance_emission_3st$n_dep == 2 &
                                                      aggr_cov_noL_Performance_emission_3st$KL_div == 7])
mean(aggr_cov_noL_Performance_emission_3st$cov_mean[aggr_cov_noL_Performance_emission_3st$n_dep == 4 &
                                                      aggr_cov_noL_Performance_emission_3st$KL_div == 7])
mean(aggr_cov_noL_Performance_emission_3st$cov_mean[aggr_cov_noL_Performance_emission_3st$n_dep == 8 &
                                                      aggr_cov_noL_Performance_emission_3st$KL_div == 7])

mean(aggr_cov_noL_Performance_emission_3st$cov_mean[aggr_cov_noL_Performance_emission_3st$n == 15])
mean(aggr_cov_noL_Performance_emission_3st$cov_mean[aggr_cov_noL_Performance_emission_3st$n == 30])
mean(aggr_cov_noL_Performance_emission_3st$cov_mean[aggr_cov_noL_Performance_emission_3st$n == 60])
mean(aggr_cov_noL_Performance_emission_3st$cov_mean[aggr_cov_noL_Performance_emission_3st$n == 120])

mean(aggr_cov_noL_Performance_emission_3st$cov_mean[aggr_cov_noL_Performance_emission_3st$n_t == 50 &
                                                      aggr_cov_noL_Performance_emission_3st$n_dep == 2])
mean(aggr_cov_noL_Performance_emission_3st$cov_mean[aggr_cov_noL_Performance_emission_3st$n_t == 100 &
                                                      aggr_cov_noL_Performance_emission_3st$n_dep == 2])
mean(aggr_cov_noL_Performance_emission_3st$cov_mean[aggr_cov_noL_Performance_emission_3st$n_t == 200 &
                                                      aggr_cov_noL_Performance_emission_3st$n_dep == 2])
mean(aggr_cov_noL_Performance_emission_3st$cov_mean[(aggr_cov_noL_Performance_emission_3st$n_t == 400 &
                                                       aggr_cov_noL_Performance_emission_3st$n_dep == 2)])
mean(aggr_cov_noL_Performance_emission_3st$cov_mean[(aggr_cov_noL_Performance_emission_3st$n_t == 800 &
                                                       aggr_cov_noL_Performance_emission_3st$n_dep == 2)])



### inspecting Coverage emission SDS ####
ggplot(data = aggr_cov_noL_Performance_emission_3st, mapping = aes(x = n_t, y = cov_SD, group = as.factor(n), color = as.factor(n))) +
  geom_point() + 
  geom_line() +
  ylab("Coverage emission SDs") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep))) +
  xlab("number of observations per subject") +
  geom_hline(yintercept=c(.90, .95), linetype="dashed", color = "grey") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

### inspecting Coverage emission means ####
ggplot(data = aggr_cov_noL_Performance_emission_3st, mapping = aes(x = n_t, y = cov_mean, group = as.factor(n), color = as.factor(n))) +
  geom_point() + 
  geom_line() +
  ylab("Coverage emission MEANS") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep))) +
  xlab("number of observations per subject") +
  geom_hline(yintercept=c(.90, .95), linetype="dashed", color = "grey") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))
