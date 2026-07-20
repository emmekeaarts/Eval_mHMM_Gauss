# jonashaslbeck@protonmail.com; July 20th, 2026

# Plot prior sensitivity figures with ggplot2.

library(ggplot2)
library(dplyr)

dir.create("Figures/ggplot2", recursive = TRUE, showWarnings = FALSE)

plot_data <- readRDS("6_Sensitivity/result_tables/SEN_plot_data.RDS")

aggr_label2 <- plot_data$aggr_label2
aggr_bias_noL_Performance_emission_Sensitivity <- plot_data$aggr_bias_noL_Performance_emission_Sensitivity
bias_noL_summary <- plot_data$bias_noL_summary
aggr_cov_noL_Performance_emission_Sensitivity <- plot_data$aggr_cov_noL_Performance_emission_Sensitivity
cov_noL_summary <- plot_data$cov_noL_summary

plot_sensitivity_gg <- function(data, yvar, ylab, hline = NULL) {
  p <- ggplot(data = data,
              mapping = aes(x = n_t,
                            y = .data[[yvar]],
                            group = interaction(as.factor(n), as.factor(mu_sc),
                                                as.factor(K0_sc), as.factor(var_tau_sig)),
                            shape = as.factor(K0_sc),
                            color = interaction(as.factor(n)),
                            linetype = interaction(as.factor(mu_sc), as.factor(var_tau_sig)))) +
    geom_point() +
    geom_line() +
    ylab(ylab) +
    scale_color_discrete(name = "Number of\nsubjects") +
    scale_x_continuous(trans = "log2", breaks = c(50, 100, 200, 400, 800)) +
    facet_grid(rows = vars(as.factor(KL_div)), cols = vars(block)) +
    xlab("number of observations per subject") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

  if(!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = "dashed", color = "grey")
  }

  p
}

plot_transition_gg <- function(data, diag_value, yvar, ylab, hline = NULL, hline_solid = NULL) {
  p <- data %>%
    filter(Diag == diag_value) %>%
    ggplot(aes(x = n_t,
               y = .data[[yvar]],
               group = interaction(as.factor(n), as.factor(mu_sc),
                                   as.factor(K0_sc), as.factor(var_tau_sig)),
               shape = as.factor(K0_sc),
               color = as.factor(n),
               linetype = interaction(as.factor(mu_sc), as.factor(var_tau_sig)))) +
    geom_point() +
    geom_line() +
    scale_x_continuous(trans = "log2", breaks = c(50, 100, 200, 400, 800)) +
    facet_grid(rows = vars(as.factor(KL_div)), cols = vars(block)) +
    labs(title = paste("Diag =", diag_value),
         x = "Number of observations per subject",
         y = ylab) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

  if(!is.null(hline_solid)) {
    p <- p + geom_hline(yintercept = hline_solid, linetype = "solid", color = "grey")
  }
  if(!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = "dashed", color = "grey")
  }

  p
}

p_label_switch_gg <- plot_sensitivity_gg(aggr_label2,
                                         yvar = "present",
                                         ylab = "Proportion of label switching")
ggsave("Figures/ggplot2/SEN_Figure_1_LabelSwitching_ggplot.pdf",
       plot = p_label_switch_gg, width = 8.5, height = 5.5)

p_bias_emission_means_gg <- plot_sensitivity_gg(aggr_bias_noL_Performance_emission_Sensitivity,
                                                yvar = "abs_rel_bias",
                                                ylab = "Bias emission means",
                                                hline = .10)
ggsave("Figures/ggplot2/SEN_Figure_2_Bias_EmissionMeans_ggplot.pdf",
       plot = p_bias_emission_means_gg, width = 8.5, height = 5.5)

p_bias_emission_sds_gg <- plot_sensitivity_gg(aggr_bias_noL_Performance_emission_Sensitivity,
                                              yvar = "SD_rel_bias",
                                              ylab = "Bias: Emission SDs",
                                              hline = .10)
ggsave("Figures/ggplot2/SEN_Figure_3_Bias_EmissionSDs_ggplot.pdf",
       plot = p_bias_emission_sds_gg, width = 8.5, height = 5.5)

p_bias_transitions_diag_gg <- plot_transition_gg(bias_noL_summary,
                                                diag_value = 1,
                                                yvar = "mean_bias",
                                                ylab = "Mean bias transition probabilities",
                                                hline = c(-.04, .04),
                                                hline_solid = 0)
ggsave("Figures/ggplot2/SEN_Figure_4a_Bias_TransitionProbabilities_Diagonal_ggplot.pdf",
       plot = p_bias_transitions_diag_gg, width = 8.5, height = 5.5)

p_bias_transitions_offdiag_gg <- plot_transition_gg(bias_noL_summary,
                                                   diag_value = 0,
                                                   yvar = "mean_bias",
                                                   ylab = "Mean bias transition probabilities",
                                                   hline = c(-.04, .04),
                                                   hline_solid = 0)
ggsave("Figures/ggplot2/SEN_Figure_4b_Bias_TransitionProbabilities_OffDiagonal_ggplot.pdf",
       plot = p_bias_transitions_offdiag_gg, width = 8.5, height = 5.5)

p_coverage_emission_means_gg <- plot_sensitivity_gg(aggr_cov_noL_Performance_emission_Sensitivity,
                                                    yvar = "cov_mean",
                                                    ylab = "Coverage Emission Means",
                                                    hline = c(.90, .95))
ggsave("Figures/ggplot2/SEN_Figure_5_Coverage_EmissionMeans_ggplot.pdf",
       plot = p_coverage_emission_means_gg, width = 8.5, height = 5.5)

p_coverage_emission_sds_gg <- plot_sensitivity_gg(aggr_cov_noL_Performance_emission_Sensitivity,
                                                  yvar = "cov_SD",
                                                  ylab = "Coverage emission distribution - SD",
                                                  hline = c(.90, .95))
ggsave("Figures/ggplot2/SEN_Figure_6_Coverage_EmissionSDs_ggplot.pdf",
       plot = p_coverage_emission_sds_gg, width = 8.5, height = 5.5)

p_coverage_transitions_diag_gg <- plot_transition_gg(cov_noL_summary,
                                                    diag_value = 1,
                                                    yvar = "cov_trans",
                                                    ylab = "Coverage transition probabilities",
                                                    hline = c(.90, .95),
                                                    hline_solid = 0)
ggsave("Figures/ggplot2/SEN_Figure_7a_Coverage_TransitionProbabilities_Diagonal_ggplot.pdf",
       plot = p_coverage_transitions_diag_gg, width = 8.5, height = 5.5)

p_coverage_transitions_offdiag_gg <- plot_transition_gg(cov_noL_summary,
                                                       diag_value = 0,
                                                       yvar = "cov_trans",
                                                       ylab = "Coverage transition probabilities",
                                                       hline = c(.90, .95),
                                                       hline_solid = 0)
ggsave("Figures/ggplot2/SEN_Figure_7b_Coverage_TransitionProbabilities_OffDiagonal_ggplot.pdf",
       plot = p_coverage_transitions_offdiag_gg, width = 8.5, height = 5.5)
