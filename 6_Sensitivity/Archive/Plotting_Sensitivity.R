# jonashaslbeck@protonmail.com; June 19th, 2026

# --------------------------------------------------------
# ---------- What is happening here? ---------------------
# --------------------------------------------------------

# Taking Emmeke's output of the prior sensitivity analysis 
# and making plots in the same style as in the rest of the paper


# Collecting notes for Emmeke:
# 1) The postprocessing file is not reproducible; all the output files are not in the "data" folder, which does not exist
# 2) The Readme.rmd does not fully correspond to the R-code/files, for example it mentions a "emission_performance_Sensitivity.RDS" which does not exist
# 3) 


# Notes for myself
# What we are plotting is: the effect of the three different blocks of prior spec variations on
# 1) Relative Bias emission Means
# 2) Abs Mean Bias transition Means
# 3) 

# --------------------------------------------------------
# ---------- Load Packages -------------------------------
# --------------------------------------------------------

library(ggplot2)
library(dplyr)
library(colorspace)

dir.create("Figures/BaseR", recursive = TRUE, showWarnings = FALSE)
dir.create("Figures/ggplot2", recursive = TRUE, showWarnings = FALSE)



# --------------------------------------------------------
# ---------- Load Preprocessed Output --------------------
# --------------------------------------------------------

# Label Switching
Label_switch_proxy_Sensitivity <- readRDS("result_tables/Label_switch_proxy_Sensitivity.RDS")

# ----- Bias -----
# Emissions
Performance_emission_Sensitivity <- readRDS(file = "result_tables/Performance_emission_Sensitivity.RDS")
# Transition Probabilities
Performance_gamma_Sensitivity <- readRDS(file = "result_tables/Performance_gamma_Sensitivity.RDS")

# XX




# --------------------------------------------------------
# ---------- More Code of Emmeke -------------------------
# --------------------------------------------------------


# --------------------------------------------------------
# ---------- Plotting: Label Switching -------------------
# --------------------------------------------------------

aggr_Label_switch_proxy_Sensitivity <- aggregate(Label_switch_proxy_Sensitivity, by = list(Label_switch_proxy_Sensitivity$sim_iteration, Label_switch_proxy_Sensitivity$var_tau_sig,
                                                                                           Label_switch_proxy_Sensitivity$trans_K0_sc, Label_switch_proxy_Sensitivity$trans_unif,
                                                                                           Label_switch_proxy_Sensitivity$mu_K0_sc, Label_switch_proxy_Sensitivity$mu_sc,
                                                                                           Label_switch_proxy_Sensitivity$KL_div,
                                                                                           Label_switch_proxy_Sensitivity$n_t, Label_switch_proxy_Sensitivity$n), FUN = mean)


aggr_Label_switch_proxy_Sensitivity$present <- (aggr_Label_switch_proxy_Sensitivity$RMSE < 0.20 ) * 1

aggr_label2 <- aggregate(aggr_Label_switch_proxy_Sensitivity, by = list(aggr_Label_switch_proxy_Sensitivity$var_tau_sig,
                                                                        aggr_Label_switch_proxy_Sensitivity$trans_K0_sc, 
                                                                        aggr_Label_switch_proxy_Sensitivity$trans_unif,
                                                                        aggr_Label_switch_proxy_Sensitivity$mu_K0_sc, 
                                                                        aggr_Label_switch_proxy_Sensitivity$mu_sc,
                                                                        aggr_Label_switch_proxy_Sensitivity$KL_div, 
                                                                        aggr_Label_switch_proxy_Sensitivity$n_t, 
                                                                        aggr_Label_switch_proxy_Sensitivity$n), FUN = mean)

aggr_label2$block[aggr_label2$var_tau_sig > 0] <- "BL3_var"
aggr_label2$block[aggr_label2$trans_unif > 0] <- "BL2_transitions"
aggr_label2$block[aggr_label2$mu_sc > 0] <- "BL1_emissions"

aggr_label2$block <- as.factor(aggr_label2$block)  
aggr_label2$K0_sc <-aggr_label2$mu_K0_sc + aggr_label2$trans_K0_sc

aggr_label2 <- aggr_label2[,-c(1:17)]


# Figure 1: Label Switching
p_label_switch_gg <- ggplot(data = aggr_label2, mapping = aes(x = n_t, y = present, group = interaction(as.factor(n), as.factor(mu_sc), 
                                                                                   as.factor(K0_sc), as.factor(var_tau_sig)), shape = as.factor(K0_sc), color = interaction(as.factor(n)), linetype = interaction(as.factor(mu_sc), as.factor(var_tau_sig)))) +
  geom_point() + 
  geom_line() +
  ylab("Proportion of label switching") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(block)) +
  xlab("number of observations per subject") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

print(p_label_switch_gg)

ggsave(filename = "Figures/ggplot2/SEN_Figure_1_LabelSwitching_ggplot.pdf",
       plot = p_label_switch_gg,
       width = 8.5,
       height = 5.5)



# --------------------------------------------------------
# ---------- Plotting: Bias ------------------------------
# --------------------------------------------------------

true_m  <- 3
n_dep   <- 4

# ----- Figure 2: Relative Bias: Emission Means -----
Performance_emission_Sensitivity <- readRDS(file = "result_tables/Performance_emission_Sensitivity.RDS")

Performance_emission_Sensitivity$abs_rel_bias <- abs(Performance_emission_Sensitivity$mean_hat - Performance_emission_Sensitivity$mean_true) / abs(Performance_emission_Sensitivity$mean_true)
Performance_emission_Sensitivity$SD_rel_bias <- abs(Performance_emission_Sensitivity$SD_hat - Performance_emission_Sensitivity$SD_true) / abs(Performance_emission_Sensitivity$SD_true)
aggr_Performance_emission_Sensitivity <- aggregate(Performance_emission_Sensitivity, by = list(Performance_emission_Sensitivity$var_tau_sig,
                                                                                               Performance_emission_Sensitivity$trans_K0_sc,
                                                                                               Performance_emission_Sensitivity$trans_unif,
                                                                                               Performance_emission_Sensitivity$mu_K0_sc,
                                                                                               Performance_emission_Sensitivity$mu_sc,
                                                                                               Performance_emission_Sensitivity$KL_div, 
                                                                                               Performance_emission_Sensitivity$n_t, 
                                                                                               Performance_emission_Sensitivity$n), FUN = median)


aggr_Performance_emission_Sensitivity$block[aggr_Performance_emission_Sensitivity$var_tau_sig > 0] <- "BL3_var"
aggr_Performance_emission_Sensitivity$block[aggr_Performance_emission_Sensitivity$trans_unif > 0] <- "BL2_transitions"
aggr_Performance_emission_Sensitivity$block[aggr_Performance_emission_Sensitivity$mu_sc > 0] <- "BL1_emissions"

aggr_Performance_emission_Sensitivity$block <- as.factor(aggr_Performance_emission_Sensitivity$block)  
aggr_Performance_emission_Sensitivity$K0_sc <-aggr_Performance_emission_Sensitivity$mu_K0_sc + aggr_Performance_emission_Sensitivity$trans_K0_sc

# Plotting: With Label Switching
# ggplot(data = aggr_Performance_emission_Sensitivity, mapping = aes(x = n_t, y = abs_rel_bias, group = interaction(as.factor(n), as.factor(mu_sc), 
#                                                                                                                   as.factor(K0_sc), as.factor(var_tau_sig)), shape = as.factor(K0_sc), color = interaction(as.factor(n)), linetype = interaction(as.factor(mu_sc), as.factor(var_tau_sig)))) +
#   geom_point() + 
#   geom_line() +
#   ylab("Bias emission means") + 
#   scale_color_discrete(name = "Number of\nsubjects") +
#   scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
#   facet_grid(rows = vars(as.factor(KL_div)), cols = vars(block)) +
#   xlab("number of observations per subject") +
#   geom_hline(yintercept=c(.10), linetype="dashed", color = "grey") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

###########################
#### correcting for label switching ####
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
  
),]

order_aggr_Label_switch_proxy_Sensitivity <- aggr_Label_switch_proxy_Sensitivity[order( aggr_Label_switch_proxy_Sensitivity$var_tau_sig,
                                                                                        aggr_Label_switch_proxy_Sensitivity$trans_K0_sc,
                                                                                        aggr_Label_switch_proxy_Sensitivity$trans_unif,
                                                                                        aggr_Label_switch_proxy_Sensitivity$mu_K0_sc,
                                                                                        aggr_Label_switch_proxy_Sensitivity$mu_sc,
                                                                                        aggr_Label_switch_proxy_Sensitivity$KL_div, 
                                                                                        aggr_Label_switch_proxy_Sensitivity$n_t, 
                                                                                        aggr_Label_switch_proxy_Sensitivity$n,
                                                                                        aggr_Label_switch_proxy_Sensitivity$sim_iteration
),]



# View(bias_noL_Performance_emission_Sensitivity)

bias_noL_Performance_emission_Sensitivity$RMSE <- rep(order_aggr_Label_switch_proxy_Sensitivity$RMSE, 
                                                      each = n_dep * true_m)

bias_noL_Performance_emission_Sensitivity <- bias_noL_Performance_emission_Sensitivity[bias_noL_Performance_emission_Sensitivity$RMSE > 0.20,]

#### inspecting absolute relative bias for gaussian emission means ####
aggr_bias_noL_Performance_emission_Sensitivity <- aggregate(bias_noL_Performance_emission_Sensitivity, by = list(bias_noL_Performance_emission_Sensitivity$var_tau_sig,
                                                                                                                 bias_noL_Performance_emission_Sensitivity$trans_K0_sc,
                                                                                                                 bias_noL_Performance_emission_Sensitivity$trans_unif,
                                                                                                                 bias_noL_Performance_emission_Sensitivity$mu_K0_sc,
                                                                                                                 bias_noL_Performance_emission_Sensitivity$mu_sc,
                                                                                                                 bias_noL_Performance_emission_Sensitivity$KL_div, 
                                                                                                                 bias_noL_Performance_emission_Sensitivity$n_t, 
                                                                                                                 bias_noL_Performance_emission_Sensitivity$n), FUN = median)


aggr_bias_noL_Performance_emission_Sensitivity$block[aggr_bias_noL_Performance_emission_Sensitivity$var_tau_sig > 0] <- "BL3_var"
aggr_bias_noL_Performance_emission_Sensitivity$block[aggr_bias_noL_Performance_emission_Sensitivity$trans_unif > 0] <- "BL2_transitions"
aggr_bias_noL_Performance_emission_Sensitivity$block[aggr_bias_noL_Performance_emission_Sensitivity$mu_sc > 0] <- "BL1_emissions"

aggr_bias_noL_Performance_emission_Sensitivity$block <- as.factor(aggr_bias_noL_Performance_emission_Sensitivity$block)  
aggr_bias_noL_Performance_emission_Sensitivity$K0_sc <-aggr_bias_noL_Performance_emission_Sensitivity$mu_K0_sc + aggr_bias_noL_Performance_emission_Sensitivity$trans_K0_sc


p_bias_emission_means_gg <- ggplot(data = aggr_bias_noL_Performance_emission_Sensitivity, mapping = aes(x = n_t, y = abs_rel_bias, group = interaction(as.factor(n), as.factor(mu_sc), 
                                                                                                                           as.factor(K0_sc), as.factor(var_tau_sig)), shape = as.factor(K0_sc), color = interaction(as.factor(n)), linetype = interaction(as.factor(mu_sc), as.factor(var_tau_sig)))) +
  geom_point() + 
  geom_line() +
  ylab("Bias emission means") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(block)) +
  xlab("number of observations per subject") +
  geom_hline(yintercept=c(.10), linetype="dashed", color = "grey") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

print(p_bias_emission_means_gg)

ggsave(filename = "Figures/ggplot2/SEN_Figure_2_Bias_EmissionMeans_ggplot.pdf",
       plot = p_bias_emission_means_gg,
       width = 8.5,
       height = 5.5)



# ----- Figure 3: Bias: Emission SDs -----
#### inspecting relative bias for gaussian emission SDs ####

# Plotting
p_bias_emission_sds_gg <- ggplot(data = aggr_bias_noL_Performance_emission_Sensitivity, mapping = aes(x = n_t, y = SD_rel_bias, group = interaction(as.factor(n), as.factor(mu_sc), 
                                                                                                                          as.factor(K0_sc), as.factor(var_tau_sig)), shape = as.factor(K0_sc), color = interaction(as.factor(n)), linetype = interaction(as.factor(mu_sc), as.factor(var_tau_sig)))) +
  geom_point() + 
  geom_line() +
  ylab("Bias: Emission SDs") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(block)) +
  xlab("number of observations per subject") +
  geom_hline(yintercept=c(.10), linetype="dashed", color = "grey") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

print(p_bias_emission_sds_gg)

ggsave(filename = "Figures/ggplot2/SEN_Figure_3_Bias_EmissionSDs_ggplot.pdf",
       plot = p_bias_emission_sds_gg,
       width = 8.5,
       height = 5.5)




# ----- Figure 4: Mean Abs Bias: Transition Probabilities -----

Performance_gamma_Sensitivity$rel_bias <- (Performance_gamma_Sensitivity$gamma_ij_hat - Performance_gamma_Sensitivity$gamma_ij_true) / Performance_gamma_Sensitivity$gamma_ij_true
Performance_gamma_Sensitivity$mean_bias <- (Performance_gamma_Sensitivity$gamma_ij_hat - Performance_gamma_Sensitivity$gamma_ij_true) 
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
  
),]


bias_noL_Performance_gamma_Sensitivity$RMSE <- rep(order_aggr_Label_switch_proxy_Sensitivity$RMSE, 
                                                   each = true_m * true_m)
bias_noL_Performance_gamma_Sensitivity <- bias_noL_Performance_gamma_Sensitivity[bias_noL_Performance_gamma_Sensitivity$RMSE > 0.20,]


#### inspecting bias in transition probs ####
bias_noL_summary <- aggregate(bias_noL_Performance_gamma_Sensitivity, by = list(bias_noL_Performance_gamma_Sensitivity$Diag, 
                                                                                bias_noL_Performance_gamma_Sensitivity$var_tau_sig,
                                                                                bias_noL_Performance_gamma_Sensitivity$trans_K0_sc,
                                                                                bias_noL_Performance_gamma_Sensitivity$trans_unif,
                                                                                bias_noL_Performance_gamma_Sensitivity$mu_K0_sc,
                                                                                bias_noL_Performance_gamma_Sensitivity$mu_sc,
                                                                                bias_noL_Performance_gamma_Sensitivity$KL_div, 
                                                                                bias_noL_Performance_gamma_Sensitivity$n_t,
                                                                                bias_noL_Performance_gamma_Sensitivity$n), FUN = mean)


bias_noL_summary$block[bias_noL_summary$var_tau_sig > 0] <- "BL3_var"
bias_noL_summary$block[bias_noL_summary$trans_unif > 0] <- "BL2_transitions"
bias_noL_summary$block[bias_noL_summary$mu_sc > 0] <- "BL1_emissions"

bias_noL_summary$block <- as.factor(bias_noL_summary$block)  
bias_noL_summary$K0_sc <-bias_noL_summary$mu_K0_sc + bias_noL_summary$trans_K0_sc

# Make plotting function to plot diag/offdiag separately

plot_bias <- function(data, diag_value) {
  
  data %>%
    filter(Diag == diag_value) %>%
    ggplot(
      aes(
        x = n_t,
        y = mean_bias,
        group = interaction(
          as.factor(n),
          as.factor(mu_sc),
          as.factor(K0_sc),
          as.factor(var_tau_sig)
        ),
        shape = as.factor(K0_sc),
        color = as.factor(n),
        linetype = interaction(
          as.factor(mu_sc),
          as.factor(var_tau_sig)
        )
      )
    ) +
    geom_point() +
    geom_line() +
    geom_hline(
      yintercept = 0,
      linetype = "solid",
      color = "grey"
    ) +
    geom_hline(
      yintercept = c(-0.04, 0.04),
      linetype = "dashed",
      color = "grey"
    ) +
    scale_x_continuous(
      trans = "log2",
      breaks = c(50, 100, 200, 400, 800)
    ) +
    facet_grid(
      rows = vars(as.factor(KL_div)),
      cols = vars(block)
    ) +
    labs(
      title = paste("Diag =", diag_value),
      x = "Number of observations per subject",
      y = "Mean bias transition probabilities"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        vjust = 0.5,
        hjust = 0.5
      )
    )
}

diag_levels <- unique(bias_noL_summary$Diag)

# Figure 7a: Diagonal
p_bias_transitions_diag_gg <- plot_bias(bias_noL_summary, 1)
print(p_bias_transitions_diag_gg)

ggsave(filename = "Figures/ggplot2/SEN_Figure_4a_Bias_TransitionProbabilities_Diagonal_ggplot.pdf",
       plot = p_bias_transitions_diag_gg,
       width = 8.5,
       height = 5.5)

# Figure 7b: Off-Diagonal
p_bias_transitions_offdiag_gg <- plot_bias(bias_noL_summary, 0)
print(p_bias_transitions_offdiag_gg)

ggsave(filename = "Figures/ggplot2/SEN_Figure_4b_Bias_TransitionProbabilities_OffDiagonal_ggplot.pdf",
       plot = p_bias_transitions_offdiag_gg,
       width = 8.5,
       height = 5.5)


# --------------------------------------------------------
# ---------- Plotting: Coverage --------------------------
# --------------------------------------------------------

# ----- Figure 5: Coverage: Emission Means -----
cov_noL_Performance_emission_Sensitivity <- bias_noL_Performance_emission_Sensitivity
cov_noL_Performance_emission_Sensitivity$cov_mean <- (cov_noL_Performance_emission_Sensitivity$mean_true > cov_noL_Performance_emission_Sensitivity$mean_95_lower &
                                                        cov_noL_Performance_emission_Sensitivity$mean_true < cov_noL_Performance_emission_Sensitivity$mean_95_upper) * 1 

cov_noL_Performance_emission_Sensitivity$cov_SD <- (cov_noL_Performance_emission_Sensitivity$SD_true > cov_noL_Performance_emission_Sensitivity$SD_95_lower &
                                                      cov_noL_Performance_emission_Sensitivity$SD_true < cov_noL_Performance_emission_Sensitivity$SD_95_upper) * 1 

cov_noL_Performance_emission_Sensitivity$cov_SD_width <- (cov_noL_Performance_emission_Sensitivity$SD_95_upper - cov_noL_Performance_emission_Sensitivity$SD_95_lower)


aggr_cov_noL_Performance_emission_Sensitivity <- aggregate(cov_noL_Performance_emission_Sensitivity, by = list(cov_noL_Performance_emission_Sensitivity$KL_div, 
                                                                                                               cov_noL_Performance_emission_Sensitivity$var_tau_sig,
                                                                                                               cov_noL_Performance_emission_Sensitivity$trans_K0_sc,
                                                                                                               cov_noL_Performance_emission_Sensitivity$trans_unif,
                                                                                                               cov_noL_Performance_emission_Sensitivity$mu_K0_sc,
                                                                                                               cov_noL_Performance_emission_Sensitivity$mu_sc,
                                                                                                               cov_noL_Performance_emission_Sensitivity$n_t, 
                                                                                                               cov_noL_Performance_emission_Sensitivity$n), FUN = mean)



mean(aggr_cov_noL_Performance_emission_Sensitivity$cov_mean)
mean(aggr_cov_noL_Performance_emission_Sensitivity$cov_mean[aggr_cov_noL_Performance_emission_Sensitivity$KL_div == 5])
mean(aggr_cov_noL_Performance_emission_Sensitivity$cov_mean[aggr_cov_noL_Performance_emission_Sensitivity$KL_div == 7])

mean(aggr_cov_noL_Performance_emission_Sensitivity$cov_mean[aggr_cov_noL_Performance_emission_Sensitivity$n == 30])
mean(aggr_cov_noL_Performance_emission_Sensitivity$cov_mean[aggr_cov_noL_Performance_emission_Sensitivity$n == 120])

aggr_cov_noL_Performance_emission_Sensitivity$block[aggr_cov_noL_Performance_emission_Sensitivity$var_tau_sig > 0] <- "BL3_var"
aggr_cov_noL_Performance_emission_Sensitivity$block[aggr_cov_noL_Performance_emission_Sensitivity$trans_unif > 0] <- "BL2_transitions"
aggr_cov_noL_Performance_emission_Sensitivity$block[aggr_cov_noL_Performance_emission_Sensitivity$mu_sc > 0] <- "BL1_emissions"

aggr_cov_noL_Performance_emission_Sensitivity$block <- as.factor(aggr_cov_noL_Performance_emission_Sensitivity$block)  
aggr_cov_noL_Performance_emission_Sensitivity$K0_sc <-aggr_cov_noL_Performance_emission_Sensitivity$mu_K0_sc + aggr_cov_noL_Performance_emission_Sensitivity$trans_K0_sc


p_coverage_emission_means_gg <- ggplot(data = aggr_cov_noL_Performance_emission_Sensitivity, mapping = aes(x = n_t, y = cov_mean, group = interaction(as.factor(n), as.factor(mu_sc), 
                                                                                                                      as.factor(K0_sc), as.factor(var_tau_sig)), shape = as.factor(K0_sc), color = interaction(as.factor(n)), linetype = interaction(as.factor(mu_sc), as.factor(var_tau_sig)))) +
  geom_point() + 
  geom_line() +
  ylab("Coverage Emission Means") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(block)) +
  xlab("number of observations per subject") +
  geom_hline(yintercept=c(.90, .95), linetype="dashed", color = "grey") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

print(p_coverage_emission_means_gg)

ggsave(filename = "Figures/ggplot2/SEN_Figure_5_Coverage_EmissionMeans_ggplot.pdf",
       plot = p_coverage_emission_means_gg,
       width = 8.5,
       height = 5.5)


# ----- Figure 6: Coverage: Emission SDs -----

p_coverage_emission_sds_gg <- ggplot(data = aggr_cov_noL_Performance_emission_Sensitivity, mapping = aes(x = n_t, y = cov_SD, group = interaction(as.factor(n), as.factor(mu_sc), 
                                                                                                                    as.factor(K0_sc), as.factor(var_tau_sig)), shape = as.factor(K0_sc), color = interaction(as.factor(n)), linetype = interaction(as.factor(mu_sc), as.factor(var_tau_sig)))) +
  geom_point() + 
  geom_line() +
  ylab("Coverage emission distribution - SD") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(block)) +
  xlab("number of observations per subject") +
  geom_hline(yintercept=c(.90, .95), linetype="dashed", color = "grey") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

print(p_coverage_emission_sds_gg)

ggsave(filename = "Figures/ggplot2/SEN_Figure_6_Coverage_EmissionSDs_ggplot.pdf",
       plot = p_coverage_emission_sds_gg,
       width = 8.5,
       height = 5.5)


# ----- Figure 7: Coverage: Transition Probabilities -----
cov_noL_Performance_gamma_Sensitivity <- bias_noL_Performance_gamma_Sensitivity
cov_noL_Performance_gamma_Sensitivity$cov_trans <- (cov_noL_Performance_gamma_Sensitivity$gamma_ij_true > cov_noL_Performance_gamma_Sensitivity$gamma_ij_95_lower &
                                                      cov_noL_Performance_gamma_Sensitivity$gamma_ij_true < cov_noL_Performance_gamma_Sensitivity$gamma_ij_95_upper) * 1 


cov_noL_summary <- aggregate(cov_noL_Performance_gamma_Sensitivity, by = list(cov_noL_Performance_gamma_Sensitivity$Diag, 
                                                                              cov_noL_Performance_gamma_Sensitivity$var_tau_sig,
                                                                              cov_noL_Performance_gamma_Sensitivity$trans_K0_sc,
                                                                              cov_noL_Performance_gamma_Sensitivity$trans_unif,
                                                                              cov_noL_Performance_gamma_Sensitivity$mu_K0_sc,
                                                                              cov_noL_Performance_gamma_Sensitivity$mu_sc,
                                                                              cov_noL_Performance_gamma_Sensitivity$KL_div, 
                                                                              cov_noL_Performance_gamma_Sensitivity$n_t,
                                                                              cov_noL_Performance_gamma_Sensitivity$n), FUN = mean)


cov_noL_summary$block[cov_noL_summary$var_tau_sig > 0] <- "BL3_var"
cov_noL_summary$block[cov_noL_summary$trans_unif > 0] <- "BL2_transitions"
cov_noL_summary$block[cov_noL_summary$mu_sc > 0] <- "BL1_emissions"

cov_noL_summary$block <- as.factor(cov_noL_summary$block)  
cov_noL_summary$K0_sc <-cov_noL_summary$mu_K0_sc + cov_noL_summary$trans_K0_sc


# Make plotting function to plot diag/offdiag separately
plot_coverage <- function(data, diag_value) {
  
  data %>%
    filter(Diag == diag_value) %>%
    ggplot(
      aes(
        x = n_t,
        y = cov_trans,
        group = interaction(
          as.factor(n),
          as.factor(mu_sc),
          as.factor(K0_sc),
          as.factor(var_tau_sig)
        ),
        shape = as.factor(K0_sc),
        color = as.factor(n),
        linetype = interaction(
          as.factor(mu_sc),
          as.factor(var_tau_sig)
        )
      )
    ) +
    geom_point() +
    geom_line() +
    geom_hline(
      yintercept = 0,
      linetype = "solid",
      color = "grey"
    ) +
    geom_hline(
      yintercept = c(0.90, 0.95),
      linetype = "dashed",
      color = "grey"
    ) +
    scale_x_continuous(
      trans = "log2",
      breaks = c(50, 100, 200, 400, 800)
    ) +
    facet_grid(
      rows = vars(as.factor(KL_div)),
      cols = vars(block)
    ) +
    labs(
      title = paste("Diag =", diag_value),
      x = "Number of observations per subject",
      y = "Coverage transition probabilities"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        vjust = 0.5,
        hjust = 0.5
      )
    )
}

diag_levels <- unique(cov_noL_summary$Diag)
# 7a: Diagonal
p_coverage_transitions_diag_gg <- plot_coverage(cov_noL_summary, 1)
print(p_coverage_transitions_diag_gg)

ggsave(filename = "Figures/ggplot2/SEN_Figure_7a_Coverage_TransitionProbabilities_Diagonal_ggplot.pdf",
       plot = p_coverage_transitions_diag_gg,
       width = 8.5,
       height = 5.5)

# 7b: Off-diagonal
p_coverage_transitions_offdiag_gg <- plot_coverage(cov_noL_summary, 0)
print(p_coverage_transitions_offdiag_gg)

ggsave(filename = "Figures/ggplot2/SEN_Figure_7b_Coverage_TransitionProbabilities_OffDiagonal_ggplot.pdf",
       plot = p_coverage_transitions_offdiag_gg,
       width = 8.5,
       height = 5.5)



# --------------------------------------------------------
# ---------- Overall Settings ----------------------------
# --------------------------------------------------------

# Simulation design
Nvar <- c(50, 100, 200)
pvar <- c(30, 120)

# Colors for Subjects
cols <- qualitative_hcl(n = 4, palette = "Dark3")[c(2,4)] # To match colors in main results figures
# plot(1:4, col = cols, pch = 19, cex = 2)

names(cols) <- pvar

block_order <- c("BL1_emissions", "BL2_transitions", "BL3_var")
block_labels <- c("Emission Priors", "Transition Priors", "Heterogeneity Priors")
kld_levels <- c(5, 7)

AddPriorCondition <- function(object) {
  object$prior_condition <- NA
  object$prior_condition[object$block == "BL1_emissions" & object$mu_sc == 1] <- "Equal emission means"
  object$prior_condition[object$block == "BL1_emissions" & object$mu_sc == 2] <- "Closer emission means"
  object$prior_condition[object$block == "BL2_transitions"] <- "Equal transition probabilities"
  object$prior_condition[object$block == "BL3_var" & object$var_tau_sig == 1] <- "Low heterogeneity"
  object$prior_condition[object$block == "BL3_var" & object$var_tau_sig == 2] <- "High heterogeneity"
  object
}

aggr_label2 <- AddPriorCondition(aggr_label2)
aggr_bias_noL_Performance_emission_Sensitivity <- AddPriorCondition(aggr_bias_noL_Performance_emission_Sensitivity)
bias_noL_summary <- AddPriorCondition(bias_noL_summary)
aggr_cov_noL_Performance_emission_Sensitivity <- AddPriorCondition(aggr_cov_noL_Performance_emission_Sensitivity)
cov_noL_summary <- AddPriorCondition(cov_noL_summary)

pch_sensitivity <- c("Equal emission means" = 1,
                     "Closer emission means" = 5,
                     "Equal transition probabilities" = 20,
                     "Low heterogeneity" = 6,
                     "High heterogeneity" = 2)

lty_K0 <- c("0" = 1,
            "1" = 2)



# --------------------------------------------------------
# ---------- Overall Plotting Function -------------------
# --------------------------------------------------------

# Plotting Labels
plotLabel <- function(text, cex=1.4, srt=0) {
  par(mar=rep(0,4))
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  text(0.5, 0.5, text, cex=cex, srt=srt)
}

plotColumnLabel <- function(text, cex=1.4, panel_mar=c(4,4.6,0,0.5)) {
  par(mar=c(0,panel_mar[2],0,panel_mar[4]))
  plot.new()
  plot.window(xlim=c(1,3), ylim=c(0,1))
  text(2, 0.5, text, cex=cex)
}

plotLegendTitle <- function(text, x, y) {
  text(x, y, text, adj = c(0, 1), font = 2, cex = 1)
}

# Main plotting function
PlotSensitivity <- function(object,
                            yvar,
                            ylim,
                            ylab = NULL,
                            h_ab = NULL,
                            h_solid = NULL,
                            main = NULL,
                            panel_mar = c(4,4.6,0,0.5),
                            ylab_line = 3.3,
                            column_cex = 1.4,
                            leg_pos="topright") {
  
  
  # ----- Define Layout -----
  lmat <- rbind(c(0, 1:3, 12),
                c(4, 6:8, 12),
                c(5, 9:11, 12))
  
  lo <- layout(mat = lmat,
               widths = c(.15, 1, 1, 1, .9),
               heights = c(.15, 1, 1))
  
  # ----- Plot Labels -----
  # Cols: KLDs
  plotColumnLabel(expression("Emission Priors"), cex=column_cex, panel_mar=panel_mar)
  plotColumnLabel(expression("Transition Priors"), cex=column_cex, panel_mar=panel_mar)
  plotColumnLabel(expression("Heterogeneity Priors"), cex=column_cex, panel_mar=panel_mar)
  
  # Rows: True K
  plotLabel(expression("         D"["KL"]*" = 5"), srt=90)
  plotLabel(expression("         D"["KL"]*" = 7"), srt=90)
  
  
  # ----- Loop in Data -----
  for(kld in kld_levels) {
    for(b in seq_along(block_order)) {
      # Canvas
      par(mar=panel_mar, pty="s")
      plot.new()
      plot.window(xlim=c(1, 3), ylim=ylim)
      grid()
      axis(1, labels=Nvar, at=1:3, las=1)
      axis(2, las=2)
      if(kld==max(kld_levels)) title(xlab = expression(N[t]), line=2.4)
      if(b==1) title(ylab=ylab, line=ylab_line)
      if(!is.null(main) & kld==kld_levels[1] & b==2) title(main=main, line=0.6, cex.main=.9)
      
      plot_data <- object[object$KL_div == kld & object$block == block_order[b], ]
      if(nrow(plot_data) == 0) next
      plot_data$x_pos <- match(plot_data$n_t, Nvar)
      plot_data$line_group <- interaction(plot_data$n,
                                          plot_data$prior_condition,
                                          plot_data$K0_sc,
                                          drop = TRUE)
      
      for(g in levels(plot_data$line_group)) {
        line_data <- plot_data[plot_data$line_group == g, ]
        line_data <- line_data[order(line_data$x_pos), ]
        line_col <- cols[as.character(line_data$n[1])]
        line_lty <- lty_K0[as.character(line_data$K0_sc[1])]
        line_pch <- pch_sensitivity[line_data$prior_condition[1]]
        
        lines(x = line_data$x_pos,
              y = line_data[[yvar]],
              col = line_col,
              lty = line_lty,
              lwd = 1.2)
        
        points(x = line_data$x_pos,
               y = line_data[[yvar]],
               col = line_col,
               pch = line_pch,
               cex = 1)
      }
      
      if(!is.null(h_solid)) abline(h = h_solid, col = "grey", lty = 1)
      if(!is.null(h_ab)) abline(h = h_ab, col = "grey", lty = 2)
      
    } #end for: p
  } # end for: kld
  
  par(mar=rep(0, 4), pty="m")
  plot.new()
  plot.window(xlim=c(0, 1), ylim=c(0, 1))
  plotLegendTitle(expression(bold(N)[pers]), .08, .94)
  legend(.08, .91, legend=paste0(pvar, " Persons"),
         text.col = cols,
         col = cols,
         lwd = 1.4,
         bty="n",
         cex=1)
  plotLegendTitle(expression(bold(K)[0]), .08, .76)
  legend(.08, .73, legend=c(expression(K[0]*" = 1"),
                            expression(K[0]*" = 0.2" %.% N[pers])),
         lty = lty_K0,
         lwd = 1.4,
         bty = "n",
         cex = 1)
  plotLegendTitle("Prior variation", .08, .57)
  legend(.08, .54, legend=names(pch_sensitivity),
         pch = pch_sensitivity,
         bty = "n",
         cex = .9)
  
} # eoF


# TESTING
pdf("Figures/BaseR/SEN_Figure_1_LabelSwitching_baseR.pdf",
    width = 7.42,
    height = 3.85,
    bg = "white")
PlotSensitivity(aggr_label2, 
                yvar="present",
                ylim=c(0, 1), 
                ylab="Prop. Label Switching",
                panel_mar=c(4,3.4,0,0.5),
                ylab_line=2.45,
                column_cex=1.15)
dev.off()

pdf("Figures/BaseR/SEN_Figure_2_Bias_EmissionMeans_baseR.pdf",
    width = 7.42,
    height = 3.85,
    bg = "white")
PlotSensitivity(aggr_bias_noL_Performance_emission_Sensitivity,
                yvar="abs_rel_bias",
                ylim=c(0, .40),
                ylab="Bias emission means",
                h_ab=.10,
                panel_mar=c(4,3.4,0,0.5),
                ylab_line=2.45,
                column_cex=1.15)
dev.off()

pdf("Figures/BaseR/SEN_Figure_3_Bias_EmissionSDs_baseR.pdf",
    width = 7.42,
    height = 3.85,
    bg = "white")
PlotSensitivity(aggr_bias_noL_Performance_emission_Sensitivity,
                yvar="SD_rel_bias",
                ylim=c(0, .20),
                ylab="Bias: Emission SDs",
                h_ab=.10,
                panel_mar=c(4,3.4,0,0.5),
                ylab_line=2.45,
                column_cex=1.15)
dev.off()

pdf("Figures/BaseR/SEN_Figure_4a_Bias_TransitionProbabilities_Diagonal_baseR.pdf",
    width = 7.42,
    height = 3.85,
    bg = "white")
PlotSensitivity(bias_noL_summary[bias_noL_summary$Diag == 1, ],
                yvar="mean_bias",
                ylim=c(-.10, .05),
                ylab="Mean bias transition probabilities",
                h_ab=c(-.04, .04),
                h_solid=0,
                main="Diagonal",
                panel_mar=c(4,3.4,0,0.5),
                ylab_line=2.45,
                column_cex=1.15)
dev.off()

pdf("Figures/BaseR/SEN_Figure_4b_Bias_TransitionProbabilities_OffDiagonal_baseR.pdf",
    width = 7.42,
    height = 3.85,
    bg = "white")
PlotSensitivity(bias_noL_summary[bias_noL_summary$Diag == 0, ],
                yvar="mean_bias",
                ylim=c(-.05, .05),
                ylab="Mean bias transition probabilities",
                h_ab=c(-.04, .04),
                h_solid=0,
                main="Off-diagonal",
                panel_mar=c(4,3.4,0,0.5),
                ylab_line=2.45,
                column_cex=1.15)
dev.off()

pdf("Figures/BaseR/SEN_Figure_5_Coverage_EmissionMeans_baseR.pdf",
    width = 7.42,
    height = 3.85,
    bg = "white")
PlotSensitivity(aggr_cov_noL_Performance_emission_Sensitivity,
                yvar="cov_mean",
                ylim=c(0, 1),
                ylab="Coverage emission means",
                h_ab=c(.90, .95),
                panel_mar=c(4,3.4,0,0.5),
                ylab_line=2.45,
                column_cex=1.15)
dev.off()

pdf("Figures/BaseR/SEN_Figure_6_Coverage_EmissionSDs_baseR.pdf",
    width = 7.42,
    height = 3.85,
    bg = "white")
PlotSensitivity(aggr_cov_noL_Performance_emission_Sensitivity,
                yvar="cov_SD",
                ylim=c(0, 1),
                ylab="Coverage emission SDs",
                h_ab=c(.90, .95),
                panel_mar=c(4,3.4,0,0.5),
                ylab_line=2.45,
                column_cex=1.15)
dev.off()

pdf("Figures/BaseR/SEN_Figure_7a_Coverage_TransitionProbabilities_Diagonal_baseR.pdf",
    width = 7.42,
    height = 3.85,
    bg = "white")
PlotSensitivity(cov_noL_summary[cov_noL_summary$Diag == 1, ],
                yvar="cov_trans",
                ylim=c(0, 1),
                ylab="Coverage transition probabilities",
                h_ab=c(.90, .95),
                h_solid=0,
                main="Diagonal",
                panel_mar=c(4,3.4,0,0.5),
                ylab_line=2.45,
                column_cex=1.15)
dev.off()

pdf("Figures/BaseR/SEN_Figure_7b_Coverage_TransitionProbabilities_OffDiagonal_baseR.pdf",
    width = 7.42,
    height = 3.85,
    bg = "white")
PlotSensitivity(cov_noL_summary[cov_noL_summary$Diag == 0, ],
                yvar="cov_trans",
                ylim=c(0, 1),
                ylab="Coverage transition probabilities",
                h_ab=c(.90, .95),
                h_solid=0,
                main="Off-diagonal",
                panel_mar=c(4,3.4,0,0.5),
                ylab_line=2.45,
                column_cex=1.15)
dev.off()
