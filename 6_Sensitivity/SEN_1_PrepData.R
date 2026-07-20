# jonashaslbeck@protonmail.com; July 20th, 2026

# Prepare all data objects needed for the prior sensitivity plots.

true_m <- 3
n_dep <- 4

add_block <- function(object) {
  object$block[object$var_tau_sig > 0] <- "BL3_var"
  object$block[object$trans_unif > 0] <- "BL2_transitions"
  object$block[object$mu_sc > 0] <- "BL1_emissions"
  object$block <- as.factor(object$block)
  object$K0_sc <- object$mu_K0_sc + object$trans_K0_sc
  object
}

add_prior_condition <- function(object) {
  object$prior_condition <- NA
  object$prior_condition[object$block == "BL1_emissions" & object$mu_sc == 1] <- "Equal emission means"
  object$prior_condition[object$block == "BL1_emissions" & object$mu_sc == 2] <- "Closer emission means"
  object$prior_condition[object$block == "BL2_transitions"] <- "Equal transition probabilities"
  object$prior_condition[object$block == "BL3_var" & object$var_tau_sig == 1] <- "Low heterogeneity"
  object$prior_condition[object$block == "BL3_var" & object$var_tau_sig == 2] <- "High heterogeneity"
  object
}

Label_switch_proxy_Sensitivity <- readRDS("6_Sensitivity/result_tables/Label_switch_proxy_Sensitivity.RDS")
Performance_emission_Sensitivity <- readRDS("6_Sensitivity/result_tables/Performance_emission_Sensitivity.RDS")
Performance_gamma_Sensitivity <- readRDS("6_Sensitivity/result_tables/Performance_gamma_Sensitivity.RDS")

# Label switching
aggr_Label_switch_proxy_Sensitivity <- aggregate(
  Label_switch_proxy_Sensitivity,
  by = list(Label_switch_proxy_Sensitivity$sim_iteration,
            Label_switch_proxy_Sensitivity$var_tau_sig,
            Label_switch_proxy_Sensitivity$trans_K0_sc,
            Label_switch_proxy_Sensitivity$trans_unif,
            Label_switch_proxy_Sensitivity$mu_K0_sc,
            Label_switch_proxy_Sensitivity$mu_sc,
            Label_switch_proxy_Sensitivity$KL_div,
            Label_switch_proxy_Sensitivity$n_t,
            Label_switch_proxy_Sensitivity$n),
  FUN = mean
)

aggr_Label_switch_proxy_Sensitivity$present <- (aggr_Label_switch_proxy_Sensitivity$RMSE < 0.20) * 1

aggr_label2 <- aggregate(
  aggr_Label_switch_proxy_Sensitivity,
  by = list(aggr_Label_switch_proxy_Sensitivity$var_tau_sig,
            aggr_Label_switch_proxy_Sensitivity$trans_K0_sc,
            aggr_Label_switch_proxy_Sensitivity$trans_unif,
            aggr_Label_switch_proxy_Sensitivity$mu_K0_sc,
            aggr_Label_switch_proxy_Sensitivity$mu_sc,
            aggr_Label_switch_proxy_Sensitivity$KL_div,
            aggr_Label_switch_proxy_Sensitivity$n_t,
            aggr_Label_switch_proxy_Sensitivity$n),
  FUN = mean
)

aggr_label2 <- add_block(aggr_label2)
aggr_label2 <- aggr_label2[, -c(1:17)]
aggr_label2 <- add_prior_condition(aggr_label2)

# Emission bias
Performance_emission_Sensitivity$abs_rel_bias <- abs(Performance_emission_Sensitivity$mean_hat - Performance_emission_Sensitivity$mean_true) /
  abs(Performance_emission_Sensitivity$mean_true)
Performance_emission_Sensitivity$SD_rel_bias <- abs(Performance_emission_Sensitivity$SD_hat - Performance_emission_Sensitivity$SD_true) /
  abs(Performance_emission_Sensitivity$SD_true)

bias_noL_Performance_emission_Sensitivity <- Performance_emission_Sensitivity[order(
  Performance_emission_Sensitivity$var_tau_sig,
  Performance_emission_Sensitivity$trans_K0_sc,
  Performance_emission_Sensitivity$trans_unif,
  Performance_emission_Sensitivity$mu_K0_sc,
  Performance_emission_Sensitivity$mu_sc,
  Performance_emission_Sensitivity$KL_div,
  Performance_emission_Sensitivity$n_t,
  Performance_emission_Sensitivity$n,
  Performance_emission_Sensitivity$sim_iteration,
  Performance_emission_Sensitivity$dep,
  Performance_emission_Sensitivity$k
), ]

order_aggr_Label_switch_proxy_Sensitivity <- aggr_Label_switch_proxy_Sensitivity[order(
  aggr_Label_switch_proxy_Sensitivity$var_tau_sig,
  aggr_Label_switch_proxy_Sensitivity$trans_K0_sc,
  aggr_Label_switch_proxy_Sensitivity$trans_unif,
  aggr_Label_switch_proxy_Sensitivity$mu_K0_sc,
  aggr_Label_switch_proxy_Sensitivity$mu_sc,
  aggr_Label_switch_proxy_Sensitivity$KL_div,
  aggr_Label_switch_proxy_Sensitivity$n_t,
  aggr_Label_switch_proxy_Sensitivity$n,
  aggr_Label_switch_proxy_Sensitivity$sim_iteration
), ]

bias_noL_Performance_emission_Sensitivity$RMSE <- rep(order_aggr_Label_switch_proxy_Sensitivity$RMSE,
                                                      each = n_dep * true_m)
bias_noL_Performance_emission_Sensitivity <- bias_noL_Performance_emission_Sensitivity[
  bias_noL_Performance_emission_Sensitivity$RMSE > 0.20, ]

aggr_bias_noL_Performance_emission_Sensitivity <- aggregate(
  bias_noL_Performance_emission_Sensitivity,
  by = list(bias_noL_Performance_emission_Sensitivity$var_tau_sig,
            bias_noL_Performance_emission_Sensitivity$trans_K0_sc,
            bias_noL_Performance_emission_Sensitivity$trans_unif,
            bias_noL_Performance_emission_Sensitivity$mu_K0_sc,
            bias_noL_Performance_emission_Sensitivity$mu_sc,
            bias_noL_Performance_emission_Sensitivity$KL_div,
            bias_noL_Performance_emission_Sensitivity$n_t,
            bias_noL_Performance_emission_Sensitivity$n),
  FUN = median
)

aggr_bias_noL_Performance_emission_Sensitivity <- add_prior_condition(add_block(aggr_bias_noL_Performance_emission_Sensitivity))

# Transition bias
Performance_gamma_Sensitivity$rel_bias <- (Performance_gamma_Sensitivity$gamma_ij_hat - Performance_gamma_Sensitivity$gamma_ij_true) /
  Performance_gamma_Sensitivity$gamma_ij_true
Performance_gamma_Sensitivity$mean_bias <- Performance_gamma_Sensitivity$gamma_ij_hat - Performance_gamma_Sensitivity$gamma_ij_true
Performance_gamma_Sensitivity$Diag <- Performance_gamma_Sensitivity$from_state_i == Performance_gamma_Sensitivity$to_state_j * 1

bias_noL_Performance_gamma_Sensitivity <- Performance_gamma_Sensitivity[order(
  Performance_gamma_Sensitivity$var_tau_sig,
  Performance_gamma_Sensitivity$trans_K0_sc,
  Performance_gamma_Sensitivity$trans_unif,
  Performance_gamma_Sensitivity$mu_K0_sc,
  Performance_gamma_Sensitivity$mu_sc,
  Performance_gamma_Sensitivity$KL_div,
  Performance_gamma_Sensitivity$n_t,
  Performance_gamma_Sensitivity$n,
  Performance_gamma_Sensitivity$sim_iteration,
  Performance_gamma_Sensitivity$from_state_i,
  Performance_gamma_Sensitivity$to_state_j
), ]

bias_noL_Performance_gamma_Sensitivity$RMSE <- rep(order_aggr_Label_switch_proxy_Sensitivity$RMSE,
                                                   each = true_m * true_m)
bias_noL_Performance_gamma_Sensitivity <- bias_noL_Performance_gamma_Sensitivity[
  bias_noL_Performance_gamma_Sensitivity$RMSE > 0.20, ]

bias_noL_summary <- aggregate(
  bias_noL_Performance_gamma_Sensitivity,
  by = list(bias_noL_Performance_gamma_Sensitivity$Diag,
            bias_noL_Performance_gamma_Sensitivity$var_tau_sig,
            bias_noL_Performance_gamma_Sensitivity$trans_K0_sc,
            bias_noL_Performance_gamma_Sensitivity$trans_unif,
            bias_noL_Performance_gamma_Sensitivity$mu_K0_sc,
            bias_noL_Performance_gamma_Sensitivity$mu_sc,
            bias_noL_Performance_gamma_Sensitivity$KL_div,
            bias_noL_Performance_gamma_Sensitivity$n_t,
            bias_noL_Performance_gamma_Sensitivity$n),
  FUN = mean
)

bias_noL_summary <- add_prior_condition(add_block(bias_noL_summary))

# Emission coverage
cov_noL_Performance_emission_Sensitivity <- bias_noL_Performance_emission_Sensitivity
cov_noL_Performance_emission_Sensitivity$cov_mean <- (
  cov_noL_Performance_emission_Sensitivity$mean_true > cov_noL_Performance_emission_Sensitivity$mean_95_lower &
    cov_noL_Performance_emission_Sensitivity$mean_true < cov_noL_Performance_emission_Sensitivity$mean_95_upper
) * 1
cov_noL_Performance_emission_Sensitivity$cov_SD <- (
  cov_noL_Performance_emission_Sensitivity$SD_true > cov_noL_Performance_emission_Sensitivity$SD_95_lower &
    cov_noL_Performance_emission_Sensitivity$SD_true < cov_noL_Performance_emission_Sensitivity$SD_95_upper
) * 1
cov_noL_Performance_emission_Sensitivity$cov_SD_width <- (
  cov_noL_Performance_emission_Sensitivity$SD_95_upper - cov_noL_Performance_emission_Sensitivity$SD_95_lower
)

aggr_cov_noL_Performance_emission_Sensitivity <- aggregate(
  cov_noL_Performance_emission_Sensitivity,
  by = list(cov_noL_Performance_emission_Sensitivity$KL_div,
            cov_noL_Performance_emission_Sensitivity$var_tau_sig,
            cov_noL_Performance_emission_Sensitivity$trans_K0_sc,
            cov_noL_Performance_emission_Sensitivity$trans_unif,
            cov_noL_Performance_emission_Sensitivity$mu_K0_sc,
            cov_noL_Performance_emission_Sensitivity$mu_sc,
            cov_noL_Performance_emission_Sensitivity$n_t,
            cov_noL_Performance_emission_Sensitivity$n),
  FUN = mean
)

aggr_cov_noL_Performance_emission_Sensitivity <- add_prior_condition(add_block(aggr_cov_noL_Performance_emission_Sensitivity))

# Transition coverage
cov_noL_Performance_gamma_Sensitivity <- bias_noL_Performance_gamma_Sensitivity
cov_noL_Performance_gamma_Sensitivity$cov_trans <- (
  cov_noL_Performance_gamma_Sensitivity$gamma_ij_true > cov_noL_Performance_gamma_Sensitivity$gamma_ij_95_lower &
    cov_noL_Performance_gamma_Sensitivity$gamma_ij_true < cov_noL_Performance_gamma_Sensitivity$gamma_ij_95_upper
) * 1

cov_noL_summary <- aggregate(
  cov_noL_Performance_gamma_Sensitivity,
  by = list(cov_noL_Performance_gamma_Sensitivity$Diag,
            cov_noL_Performance_gamma_Sensitivity$var_tau_sig,
            cov_noL_Performance_gamma_Sensitivity$trans_K0_sc,
            cov_noL_Performance_gamma_Sensitivity$trans_unif,
            cov_noL_Performance_gamma_Sensitivity$mu_K0_sc,
            cov_noL_Performance_gamma_Sensitivity$mu_sc,
            cov_noL_Performance_gamma_Sensitivity$KL_div,
            cov_noL_Performance_gamma_Sensitivity$n_t,
            cov_noL_Performance_gamma_Sensitivity$n),
  FUN = mean
)

cov_noL_summary <- add_prior_condition(add_block(cov_noL_summary))

plot_data <- list(
  aggr_label2 = aggr_label2,
  aggr_bias_noL_Performance_emission_Sensitivity = aggr_bias_noL_Performance_emission_Sensitivity,
  bias_noL_summary = bias_noL_summary,
  aggr_cov_noL_Performance_emission_Sensitivity = aggr_cov_noL_Performance_emission_Sensitivity,
  cov_noL_summary = cov_noL_summary
)

saveRDS(plot_data, "6_Sensitivity/result_tables/SEN_plot_data.RDS")
