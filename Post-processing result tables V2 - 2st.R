devtools::install_github("emmekeaarts/mHMMbayes@one-state")
library(mHMMbayes)
library(ggplot2)
library(reshape2)

# Setting up ####
# Simulation settings 
settings4 <- read.csv("sim_scenarios_V4.csv")
true_m  <- 2
n_sim <- 100
J <- 2000
burn_in <- 500

# Reading in tables
Label_switch_proxy_2st <- readRDS("result_tables/Label_switch_proxy_2st.RDS")
Performance_emission_2st <- readRDS(file = "result_tables/Performance_emission_2st.RDS")
Model_selection_2st <- readRDS(file = "result_tables/Model_selection_2st.RDS")
Performance_gamma_2st <- readRDS(file = "result_tables/Performance_gamma_2st.RDS")



# MODEL SELECTION ####
## model selection, corrected for label switching ####
### Combining model selection table with label switching ####

# taking mean RMSE over each simulation run such that we combine RMSE over dependent variables 
aggr_Label_switch_proxy_2st <- aggregate(Label_switch_proxy_2st, by = list(Label_switch_proxy_2st$sim_iteration, 
                                                                           Label_switch_proxy_2st$KL_div,
                                                                           Label_switch_proxy_2st$n_dep, 
                                                                           Label_switch_proxy_2st$n_t, 
                                                                           Label_switch_proxy_2st$n), FUN = mean)


Model_selection_2st_AIC <- reshape(Model_selection_2st, direction = "wide", timevar = "K_candidate", 
                                   idvar = c("sim_iteration", "true_m", "n", "n_t", "n_dep", "KL_div"), drop = "AICc")

Model_selection_2st_AIC$selected <- apply(Model_selection_2st_AIC[,7:10], 1, which.min)
Model_selection_2st_AIC$correct_selected <- 0
Model_selection_2st_AIC$correct_selected[Model_selection_2st_AIC$selected == true_m] <- 1  





# Making sure ordering is same in performance table and label switching table, saving in new object 
noL_Model_selection_2st_AIC <- Model_selection_2st_AIC[order(
  Model_selection_2st_AIC$KL_div, 
  Model_selection_2st_AIC$n_dep, 
  Model_selection_2st_AIC$n_t, 
  Model_selection_2st_AIC$n,
  Model_selection_2st_AIC$sim_iteration
),]


# Bit double, but to make sure ordering is the same compared to performance table 
order_aggr_Label_switch_proxy_2st <- aggr_Label_switch_proxy_2st[order(aggr_Label_switch_proxy_2st$KL_div, 
                                                                       aggr_Label_switch_proxy_2st$n_dep, 
                                                                       aggr_Label_switch_proxy_2st$n_t, 
                                                                       aggr_Label_switch_proxy_2st$n,
                                                                       aggr_Label_switch_proxy_2st$sim_iteration
),]


# adding the RMSE to the performance table, RMSE needs to be repeated M * n_dep timees, where n_dep varies over the scenarios.. 
noL_Model_selection_2st_AIC$RMSE <- order_aggr_Label_switch_proxy_2st$RMSE

View(noL_Model_selection_2st_AIC)

# removing 1) simulation runs that have RMSE < 0.20 AND 
noL_Model_selection_2st_AIC <- noL_Model_selection_2st_AIC[noL_Model_selection_2st_AIC$RMSE > 0.20,]



### Aggregating to MEAN within one simulation setting (= percentage of iterations resulting in correct selection)
aggr_noL_Model_selection_2st_AIC <- aggregate(noL_Model_selection_2st_AIC, by = list(noL_Model_selection_2st_AIC$KL_div, 
                                                                                     noL_Model_selection_2st_AIC$n_dep, 
                                                                                     noL_Model_selection_2st_AIC$n_t, 
                                                                                     noL_Model_selection_2st_AIC$n), FUN = mean)

aggr_noL_Model_selection_2st_AIC <- aggr_noL_Model_selection_2st_AIC[,c(6:10,16)]

# checking that minimum number of simulation runs is 20, and aligning with selction of M = 3 for comparability 
table(noL_Model_selection_2st_AIC$n, noL_Model_selection_2st_AIC$n_t, noL_Model_selection_2st_AIC$n_dep, noL_Model_selection_2st_AIC$KL_div)


red_aggr_noL_Model_selection_2st_AIC <- aggr_noL_Model_selection_2st_AIC[which(aggr_noL_Model_selection_2st_AIC$KL_div == 3 &
                                                                                 ((aggr_noL_Model_selection_2st_AIC$n_t == 50) | 
                                                                                 (aggr_noL_Model_selection_2st_AIC$n_t == 100))),]

red2_aggr_noL_Model_selection_2st_AIC <- red_aggr_noL_Model_selection_2st_AIC[which((red_aggr_noL_Model_selection_2st_AIC$n_dep == 2 &
                                                                                      red_aggr_noL_Model_selection_2st_AIC$n ==15) |
                                                                                (red_aggr_noL_Model_selection_2st_AIC$n_dep == 4 &
                                                                                  ( red_aggr_noL_Model_selection_2st_AIC$n_t == 50 & 
                                                                                      red_aggr_noL_Model_selection_2st_AIC$n < 120) |
                                                                                   ( red_aggr_noL_Model_selection_2st_AIC$n_t == 100 & 
                                                                                       red_aggr_noL_Model_selection_2st_AIC$n < 30)) |
                                                                                   (red_aggr_noL_Model_selection_2st_AIC$n_dep == 8 &
                                                                                      ( red_aggr_noL_Model_selection_2st_AIC$n_t == 50) |
                                                                                      ( red_aggr_noL_Model_selection_2st_AIC$n_t == 100 & 
                                                                                          red_aggr_noL_Model_selection_2st_AIC$n < 60))),]


red3_aggr_noL_Model_selection_2st_AIC <- red2_aggr_noL_Model_selection_2st_AIC[-c(4,9:13),]

1 - mean(red3_aggr_noL_Model_selection_2st_AIC$correct_selected[red3_aggr_noL_Model_selection_2st_AIC$KL_div == 3 & red3_aggr_noL_Model_selection_2st_AIC$n_t == 50])
1 - mean(red3_aggr_noL_Model_selection_2st_AIC$correct_selected[red3_aggr_noL_Model_selection_2st_AIC$KL_div == 3 & red3_aggr_noL_Model_selection_2st_AIC$n_t == 100])


