# jonashaslbeck@protonmail.com; July 20th, 2026

# Plot prior sensitivity figures with base R.

library(colorspace)

dir.create("Figures/BaseR", recursive = TRUE, showWarnings = FALSE)

plot_data <- readRDS("result_tables/SEN_plot_data.RDS")

aggr_label2 <- plot_data$aggr_label2
aggr_bias_noL_Performance_emission_Sensitivity <- plot_data$aggr_bias_noL_Performance_emission_Sensitivity
bias_noL_summary <- plot_data$bias_noL_summary
aggr_cov_noL_Performance_emission_Sensitivity <- plot_data$aggr_cov_noL_Performance_emission_Sensitivity
cov_noL_summary <- plot_data$cov_noL_summary

Nvar <- c(50, 100, 200)
pvar <- c(30, 120)
cols <- qualitative_hcl(n = 4, palette = "Dark3")[c(2, 4)]
names(cols) <- pvar

block_order <- c("BL1_emissions", "BL2_transitions", "BL3_var")
kld_levels <- c(5, 7)

pch_sensitivity <- c("Equal emission means" = 1,
                     "Closer emission means" = 5,
                     "Equal transition probabilities" = 20,
                     "Low heterogeneity" = 6,
                     "High heterogeneity" = 2)

lty_K0 <- c("0" = 1,
            "1" = 2)

plotLabel <- function(text, cex = 1.4, srt = 0) {
  par(mar = rep(0, 4))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0.5, 0.5, text, cex = cex, srt = srt)
}

plotColumnLabel <- function(text, cex = 1.4, panel_mar = c(4, 4.6, 0, 0.5)) {
  par(mar = c(0, panel_mar[2], 0, panel_mar[4]))
  plot.new()
  plot.window(xlim = c(1, 3), ylim = c(0, 1))
  text(2, 0.5, text, cex = cex)
}

plotLegendTitle <- function(text, x, y) {
  text(x, y, text, adj = c(0, 1), font = 2, cex = 1)
}

PlotSensitivity <- function(object,
                            yvar,
                            ylim,
                            ylab = NULL,
                            h_ab = NULL,
                            h_solid = NULL,
                            main = NULL,
                            panel_mar = c(4, 4.6, 0, 0.5),
                            ylab_line = 3.3,
                            column_cex = 1.4) {
  lmat <- rbind(c(0, 1:3, 12),
                c(4, 6:8, 12),
                c(5, 9:11, 12))

  layout(mat = lmat,
         widths = c(.15, 1, 1, 1, .9),
         heights = c(.15, 1, 1))

  plotColumnLabel(expression("Emission Priors"), cex = column_cex, panel_mar = panel_mar)
  plotColumnLabel(expression("Transition Priors"), cex = column_cex, panel_mar = panel_mar)
  plotColumnLabel(expression("Heterogeneity Priors"), cex = column_cex, panel_mar = panel_mar)

  plotLabel(expression("         D"["KL"]*" = 5"), srt = 90)
  plotLabel(expression("         D"["KL"]*" = 7"), srt = 90)

  for(kld in kld_levels) {
    for(b in seq_along(block_order)) {
      par(mar = panel_mar, pty = "s")
      plot.new()
      plot.window(xlim = c(1, 3), ylim = ylim)
      grid()
      axis(1, labels = Nvar, at = 1:3, las = 1)
      axis(2, las = 2)
      if(kld == max(kld_levels)) title(xlab = expression(N[t]), line = 2.4)
      if(b == 1) title(ylab = ylab, line = ylab_line)
      if(!is.null(main) & kld == kld_levels[1] & b == 2) title(main = main, line = 0.6, cex.main = .9)

      plot_data_panel <- object[object$KL_div == kld & object$block == block_order[b], ]
      if(nrow(plot_data_panel) == 0) next

      plot_data_panel$x_pos <- match(plot_data_panel$n_t, Nvar)
      plot_data_panel$line_group <- interaction(plot_data_panel$n,
                                                plot_data_panel$prior_condition,
                                                plot_data_panel$K0_sc,
                                                drop = TRUE)

      for(g in levels(plot_data_panel$line_group)) {
        line_data <- plot_data_panel[plot_data_panel$line_group == g, ]
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
    }
  }

  par(mar = rep(0, 4), pty = "m")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  plotLegendTitle(expression(bold(N)[pers]), .08, .94)
  legend(.08, .91, legend = paste0(pvar, " Persons"),
         text.col = cols,
         col = cols,
         lwd = 1.4,
         bty = "n",
         cex = 1)
  plotLegendTitle(expression(bold(K)[0]), .08, .76)
  legend(.08, .73, legend = c(expression(K[0]*" = 1"),
                              expression(K[0]*" = 0.2" %.% N[pers])),
         lty = lty_K0,
         lwd = 1.4,
         bty = "n",
         cex = 1)
  plotLegendTitle("Prior variation", .08, .57)
  legend(.08, .54, legend = names(pch_sensitivity),
         pch = pch_sensitivity,
         bty = "n",
         cex = .9)
}

PlotSensitivityCompact <- function(file, object, yvar, ylim, ylab,
                                   h_ab = NULL, h_solid = NULL, main = NULL) {
  pdf(file, width = 7.42, height = 3.85, bg = "white")
  PlotSensitivity(object,
                  yvar = yvar,
                  ylim = ylim,
                  ylab = ylab,
                  h_ab = h_ab,
                  h_solid = h_solid,
                  main = main,
                  panel_mar = c(4, 3.4, 0, 0.5),
                  ylab_line = 2.45,
                  column_cex = 1.15)
  dev.off()
}

PlotSensitivityCompact("Figures/BaseR/SEN_Figure_1_LabelSwitching_baseR.pdf",
                       aggr_label2,
                       yvar = "present",
                       ylim = c(0, 1),
                       ylab = "Prop. Label Switching")

PlotSensitivityCompact("Figures/BaseR/SEN_Figure_2_Bias_EmissionMeans_baseR.pdf",
                       aggr_bias_noL_Performance_emission_Sensitivity,
                       yvar = "abs_rel_bias",
                       ylim = c(0, .40),
                       ylab = "Bias emission means",
                       h_ab = .10)

PlotSensitivityCompact("Figures/BaseR/SEN_Figure_3_Bias_EmissionSDs_baseR.pdf",
                       aggr_bias_noL_Performance_emission_Sensitivity,
                       yvar = "SD_rel_bias",
                       ylim = c(0, .20),
                       ylab = "Bias: Emission SDs",
                       h_ab = .10)

PlotSensitivityCompact("Figures/BaseR/SEN_Figure_4a_Bias_TransitionProbabilities_Diagonal_baseR.pdf",
                       bias_noL_summary[bias_noL_summary$Diag == 1, ],
                       yvar = "mean_bias",
                       ylim = c(-.10, .05),
                       ylab = "Mean bias transition prob.",
                       h_ab = c(-.04, .04),
                       h_solid = 0,
                       main = "Diagonal")

PlotSensitivityCompact("Figures/BaseR/SEN_Figure_4b_Bias_TransitionProbabilities_OffDiagonal_baseR.pdf",
                       bias_noL_summary[bias_noL_summary$Diag == 0, ],
                       yvar = "mean_bias",
                       ylim = c(-.05, .05),
                       ylab = "Mean bias transition prob.",
                       h_ab = c(-.04, .04),
                       h_solid = 0,
                       main = "Off-diagonal")

PlotSensitivityCompact("Figures/BaseR/SEN_Figure_5_Coverage_EmissionMeans_baseR.pdf",
                       aggr_cov_noL_Performance_emission_Sensitivity,
                       yvar = "cov_mean",
                       ylim = c(0, 1),
                       ylab = "Cov. emission means",
                       h_ab = c(.90, .95))

PlotSensitivityCompact("Figures/BaseR/SEN_Figure_6_Coverage_EmissionSDs_baseR.pdf",
                       aggr_cov_noL_Performance_emission_Sensitivity,
                       yvar = "cov_SD",
                       ylim = c(0, 1),
                       ylab = "Cov. emission SDs",
                       h_ab = c(.90, .95))

PlotSensitivityCompact("Figures/BaseR/SEN_Figure_7a_Coverage_TransitionProbabilities_Diagonal_baseR.pdf",
                       cov_noL_summary[cov_noL_summary$Diag == 1, ],
                       yvar = "cov_trans",
                       ylim = c(0, 1),
                       ylab = "Cov. transition prob.",
                       h_ab = c(.90, .95),
                       h_solid = 0,
                       main = "Diagonal")

PlotSensitivityCompact("Figures/BaseR/SEN_Figure_7b_Coverage_TransitionProbabilities_OffDiagonal_baseR.pdf",
                       cov_noL_summary[cov_noL_summary$Diag == 0, ],
                       yvar = "cov_trans",
                       ylim = c(0, 1),
                       ylab = "Cov. transition prob.",
                       h_ab = c(.90, .95),
                       h_solid = 0,
                       main = "Off-diagonal")
