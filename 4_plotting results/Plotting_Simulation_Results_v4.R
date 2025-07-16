# jonashaslbeck@protonmail.com; Feb 16th, 2024

# --------------------------------------------------------
# ---------- What is happening here? ---------------------
# --------------------------------------------------------

# Plotting all results figures based on Emmeke's (aggregate) simulation results

# --------------------------------------------------------
# ---------- Load Packages -------------------------------
# --------------------------------------------------------


library(reshape2)
library(viridis)
library(alluvial)
library(mHMMbayes)

library(coda) # to compute Gelman-Rubin statistic


### Packages Jonas
# Wrangling
library(abind)


# Plotting
library(RColorBrewer)
library(xtable)
library(scales)
library(colorspace) # color-blind

source("0_Helpers.R")


# --------------------------------------------------------
# ---------- Some Global Stuff ---------------------------
# --------------------------------------------------------

Ntvar <- c(50, 100, 200, 400, 800)
n_Ntvar <- length(Ntvar)
KLDvar <- c(3, 5, 7)

# --------------------------------------------------------
# ---------- Figure: LabelSwitching Proportion -----------
# --------------------------------------------------------

# -----------------------
# ------- For K=2 -------
# -----------------------
# We need this to subset on label switching also for the model selection results (acc/bias)

# ----- Load Data -----
res_lblsw_2st <- readRDS("result_tables/Label_switch_proxy_2st.RDS")
head(res_lblsw_2st)

# ----- Process -----
lblsw_2st <- xtabs(RMSE ~ KL_div + n_dep + n_t + n + dep + sim_iteration, data = res_lblsw_2st)
lblsw_2st[lblsw_2st==0] <- NA
dim(lblsw_2st)

# Get Proportion
lblsw_bin_2st <- lblsw_2st < 0.2

# Aggregate this across dependent variables
dim(lblsw_bin_2st)
lblsw_bin_hom_2st <- apply(lblsw_bin_2st, c(1:4, 6), function(x) {
  log_out <- any(x, na.rm=TRUE)
  N_noNA <- sum(!is.na(x))
  c(rep(log_out, N_noNA), rep(NA, 8-N_noNA))
})
dim(lblsw_bin_hom_2st)
lblsw_bin_hom_2st <- aperm(lblsw_bin_hom_2st, c(2:5, 1, 6))
dim(lblsw_bin_hom_2st)

# ------ Compute Indicator: Proportion violated across conditions? -----
# First step: Across iterations
lblsw_bin_hom_90sat_2st <- apply(lblsw_bin_hom_2st[, , , , 1, ], 1:4, function(x) mean(x, na.rm=TRUE)>0.8)
dim(lblsw_bin_hom_90sat_2st)


# -----------------------
# ------- For K=3 -------
# -----------------------

# ----- Load Data -----
res_lblsw <- readRDS("result_tables/Label_switch_proxy_3st.RDS")
head(res_lblsw)

# ----- Process -----
lblsw <- xtabs(RMSE ~ KL_div + n_dep + n_t + n + dep + sim_iteration, data = res_lblsw)
lblsw[lblsw==0] <- NA
dim(lblsw)

# # Get Proportion
# lblsw_bin <- lblsw < 0.2
#
# # Aggregate this across dependent variables
# dim(lblsw_bin)
# lblsw_bin_hom_old <- apply(lblsw_bin, c(1:4, 6), function(x) {
#   log_out <- any(x, na.rm=TRUE)
#   N_noNA <- sum(!is.na(x))
#   c(rep(log_out, N_noNA), rep(NA, 8-N_noNA))
# })
# dim(lblsw_bin_hom_old)
# lblsw_bin_hom_old <- aperm(lblsw_bin_hom_old, c(2:5, 1, 6))
# dim(lblsw_bin_hom_old)

# Apply rule: across the 2/4/8 dependent variable, take average, then apply the 0.2 threshold on those; then copy the decision into the 2/4/8 non-NA entries
dim(lblsw)
lblsw_bin_hom <- apply(lblsw, c(1:4, 6), function(x) {
  x_mean <- mean(x, na.rm=TRUE)
  log_out <- x_mean < 0.2
  N_noNA <- sum(!is.na(x))
  c(rep(log_out, N_noNA), rep(NA, 8-N_noNA))
})
dim(lblsw_bin_hom)
lblsw_bin_hom <- aperm(lblsw_bin_hom, c(2:5, 1, 6))
dim(lblsw_bin_hom)


# ----- Plot -----

cols <- qualitative_hcl(n = 4, palette = "Dark3")

# Aggregate for Plot
lblsw_bin_hom_agg_iter <- apply(lblsw_bin_hom, 1:4, function(x) mean(x, na.rm=TRUE))

dim(lblsw_bin_hom_agg_iter)

lblsw_bin_hom_agg_iter[3, 1, 1, 1]
lblsw_bin_hom_agg_iter[3, 3, 1, 1]

sc <- 0.72
pdf("Figures/Fig_Res_LabelSwitching_prop_k=3.pdf", width = 8*sc, height = 8*sc*0.95)

PlotSimEmiss(object = lblsw_bin_hom_agg_iter,
             ylim = c(0, 1),
             ylab="Prop Label Switching",
             h_ab=NULL,
             leg_pos="topright")

dev.off()

# ------ Compute Indicator: Proportion violated across conditions? -----
# First step: Across iterations
# If there is more than 80% label switching then NA; that is: we require 20% or more
# Note that we subset on dep=1, because all of them contain the same values (see above calulation of lblsw_bin_hom)
lblsw_bin_hom_90sat <- apply(lblsw_bin_hom[, , , , 1, ], 1:4, function(x) mean(x, na.rm=TRUE)>0.8)
dim(lblsw_bin_hom_90sat)



# --------------------------------------------------------
# ---------- Load Results: Convergence Emissions ---------
# --------------------------------------------------------

# Load Data
res_conv_em <- readRDS("result_tables/Convergence_emission_3st_corrected_V2.RDS")
# unique(res_conv$k) # This is the state; of course we have n_dep means PER state
unique(res_conv_em$sim_rep) # Four estimates of the GR-statistic, from four sim iterations (which themselves are based on I think 4 runs, of which we requite ??? to be not affected by label switching)

dim(res_conv_em)
table(res_conv_em$valid_runs)

# Subset: Need 3 or 4 / 4 iterations, to compute
res_conv_em_ss <- res_conv_em[res_conv_em$valid_runs > 2, ]

# ----- GR-Statistic -----
# Rework
conv_GR_em <- xtabs(GR_stat ~ KL_div + n_dep + n_t + n + true_m + sim_rep + dep + k + valid_runs, data = res_conv_em_ss)
conv_GR_em[conv_GR_em==0] <- NA

# Aggregate
conv_GR_em_agg <- apply(conv_GR_em, 1:4, mean, na.rm=TRUE)

# Compute Proportion Converged
conv_GR_em_prop <- conv_GR_em
conv_GR_em_prop[conv_GR_em_prop < 1.1] <- TRUE
conv_GR_em_prop[conv_GR_em_prop > 1.1] <- FALSE
# Aggregate
conv_GR_em_prop_agg <- apply(conv_GR_em_prop, 1:4, mean, na.rm=TRUE)


# ----- Convergence: Transition Probabilities -----
# Load Data
res_conv_trns <- readRDS("result_tables/Convergence_gamma_3st_corrected_V2.RDS")
head(res_conv_trns)

# Subset: Need 3 or 4 / 4 iterations, to compute
res_conv_trns_ss <- res_conv_trns[res_conv_trns$valid_runs > 2, ]


# Rework
conv_GR_trn <- xtabs(GR_stat ~ KL_div + n_dep + n_t  + n + from_state_i + to_state_j + true_m + sim_rep, data = res_conv_trns_ss)
conv_GR_trn[conv_GR_trn==0] <- NA
dim(conv_GR_trn)
# Compute Proportion Converged
conv_GR_trn_prop <- conv_GR_trn
conv_GR_trn_prop[conv_GR_trn_prop < 1.1] <- TRUE
conv_GR_trn_prop[conv_GR_trn_prop > 1.1] <- FALSE

# Aggregate
conv_GR_trns_agg <- apply(conv_GR_trn, 1:6, mean, na.rm=TRUE)
conv_GR_trn_prop_agg <- apply(conv_GR_trn_prop, 1:6, mean, na.rm=TRUE)
dim(conv_GR_trns_agg)


# --------------------------------------------------------
# ---------- Figure: Convergence: Emission (GR) -----------
# --------------------------------------------------------

sc <- 0.85
pdf("Figures/Fig_Res_Conv_k=3_Emission_GR_NEW_V2.pdf", width = 8*sc, height = 8*sc)

PlotSimEmiss(object = conv_GR_em_agg,
             ylim = c(0.8, 1.6),
             ylab = "Gelman-Ruben Statistic",
             h_ab=c(1, 1.1),
             leg_pos="topright")

dev.off()


# ----- KLD = 2,3 ------
sc <- 0.72
pdf("Figures/Fig_Res_Conv_k=3_Emission_GR_NEW_KLD23_V2.pdf", width = 8*sc, height = 8*sc*(2/3)*0.95)

PlotSimEmiss_KLD23(object = conv_GR_em_agg,
                   ylim = c(0.8, 1.6),
                   ylab = "Gelman-Ruben Statistic",
                   h_ab=c(1, 1.1),
                   leg_pos="topright")

dev.off()




# --------------------------------------------------------
# ---------- Figure: Convergence: Emission (Prop) --------
# --------------------------------------------------------

sc <- 0.85
pdf("Figures/Fig_Res_Conv_k=3_Emission_Prop_NEW.pdf", width = 8*sc, height = 8*sc)

PlotSimEmiss(object = conv_GR_em_prop_agg,
             ylim = c(0, 1),
             ylab = "Proportion Converged",
             # h_ab=c(1, 1.1),
             leg_pos="topright")

dev.off()

# ----- KLD = 2,3 ------
sc <- 0.72
pdf("Figures/Fig_Res_Conv_k=3_Emission_Prop_NEW_KLD23_V2.pdf", width = 8*sc, height = 8*sc*(2/3)*0.95)

PlotSimEmiss_KLD23(object = conv_GR_em_prop_agg,
                   ylim = c(0, 1),
                   ylab = "Proportion Convergence",
                   h_ab=c(1, 1.1),
                   leg_pos="topright",
                   leg_box = c(1,2))

dev.off()



# --------------------------------------------------------
# ---------- Figure: Convergence: Transition Prob  (GR) --
# --------------------------------------------------------

sc <- 0.85
pdf("Figures/Fig_Res_Conv_k=3_Transition_GR_NEW.pdf", width = 8*sc, height = 8*sc)

PlotSimTrans(object = conv_GR_trns_agg,
             ylim = c(0.8, 1.6),
             h_ab=c(1, 1.1),
             ylab = "Gelman-Ruben Statistic",
             leg_pos="topright",
             leg_pos2 = "left")

dev.off()


# ----- KLD = 2,3 ------
sc <- 0.72
pdf("Figures/Fig_Res_Conv_k=3_Transition_GR_NEW_KLD23_V2.pdf", width = 8*sc, height = 8*sc*(2/3)*0.95)

PlotSimTrans_KLD23(object = conv_GR_trns_agg,
                   ylim = c(0.8, 1.6),
                   h_ab=c(1, 1.1),
                   ylab = "Gelman-Ruben Statistic",
                   leg_pos="topright",
                   leg_pos2 = "left",
                   leg2_box=c(2,3))

dev.off()


# --------------------------------------------------------
# ---------- Figure: Convergence: Transition Prob (Prop) -
# --------------------------------------------------------

sc <- 0.85
pdf("Figures/Fig_Res_Conv_k=3_Transition_Prop_NEW.pdf", width = 8*sc, height = 8*sc)

PlotSimTrans(object = conv_GR_trn_prop_agg,
             ylim = c(0, 1),
             # h_ab=c(1, 1.1),
             ylab = "Proportion Converged",
             leg_pos="right",
             leg_pos2 = "bottomleft")

dev.off()


# ----- KLD = 2,3 ------
sc <- 0.72
pdf("Figures/Fig_Res_Conv_k=3_Transition_Prop_NEW_KLD23_V2.pdf", width = 8*sc, height = 8*sc*(2/3)*0.95)

PlotSimTrans_KLD23(object = conv_GR_trn_prop_agg,
                   ylim = c(0, 1),
                   h_ab=c(1, 1.1),
                   ylab = "Proportion Converged",
                   leg_pos="right",
                   leg_pos2 = "bottomleft")

dev.off()



# --------------------------------------------------------
# ---------- Load Results: Model Selection ---------------
# --------------------------------------------------------

results_1st_new <- readRDS("result_tables/Model_selection_1st.RDS")
results_2st_new <- readRDS("result_tables/Model_selection_2st.RDS")
results_3st_new <- readRDS("result_tables/Model_selection_3st.RDS")

l_res <- list(results_1st_new,
              results_2st_new,
              results_3st_new)

# --------------------------------------------------------
# ---------- Preprocessing: Model Selection: AICc ---------
# --------------------------------------------------------

# ----- Collapse into one Array -----

# Each flat table into array
l_ms <- list()
for(k in 1:3) {
  l_ms[[k]] <- xtabs(AICc ~ K_candidate + KL_div + n_dep + n_t + n + sim_iteration, data = l_res[[k]]) # Doesn't work because I have different dependent variables
  l_ms[[k]][l_ms[[k]] == 0] <- NA # because xtabs outrageously inserts 0s when cmb are missing
} # end for: k
D <- length(dim(l_ms[[k]]))

# ----- Make Predictions -----
l_ms_pred <- list()
for(k in 1:3) l_ms_pred[[k]] <- apply(l_ms[[k]] , 2:D, function(x) {
  if(any(is.na(x))) NA else which.min(x)
})

# ----- Compute Accuracy -----
l_ms_acc <- list()
for(k in 1:3) l_ms_acc[[k]] <- l_ms_pred[[k]] == k
# Combine in our array
a_ms_acc <- abind(l_ms_acc[[1]],
                  l_ms_acc[[2]],
                  l_ms_acc[[3]],
                  along = length(c(3, dim(l_ms_acc[[k]]))))

### Subset on "no Label Switching"
a_ms_acc_noL <- a_ms_acc
a_ms_acc_noL[, , , , , 2][lblsw_bin_hom_2st[, , , , 1, ]] <- NA # for K=2 (for K=1 not needed, because no label switching possible)
a_ms_acc_noL[, , , , , 3][lblsw_bin_hom[, , , , 1, ]] <- NA # for K=3
### Aggregate
a_ms_acc_agg_it <- apply(a_ms_acc, c(1:4, 6), mean, na.rm=TRUE)
a_ms_acc_noL_agg_it <- apply(a_ms_acc_noL, c(1:4, 6), mean, na.rm=TRUE)
### Apply 90% rule
a_ms_acc_noL_agg_it[, , , , 2][lblsw_bin_hom_90sat_2st] <- NA
a_ms_acc_noL_agg_it[, , , , 3][lblsw_bin_hom_90sat] <- NA
### Aggregate across p and Nsubj
a_ms_acc_agg_it_p_Nsubj <- apply(a_ms_acc_agg_it, c(1, 3, 5), mean, na.rm=TRUE)
a_ms_acc_noL_agg_it_p_Nsubj <- apply(a_ms_acc_noL_agg_it, c(1, 3, 5), mean, na.rm=TRUE)
### Cmb in list
l_ms_acc_agg_it_p_Nsubj_AICc <- list(a_ms_acc_agg_it_p_Nsubj,
                                     a_ms_acc_noL_agg_it_p_Nsubj)

# ----- Compute Bias -----
l_ms_bias <- list()
for(k in 1:3) l_ms_bias[[k]] <- l_ms_pred[[k]] - k
# Combine in our array
a_ms_bias <- abind(l_ms_bias[[1]],
                   l_ms_bias[[2]],
                   l_ms_bias[[3]],
                   along = length(c(3, dim(l_ms_bias[[k]]))))

### Subset on "no Label Switching"
a_ms_bias_noL <- a_ms_bias
a_ms_bias_noL[, , , , , 2][lblsw_bin_hom_2st[, , , , 1, ]] <- NA # for K=2 (for K=1 not needed, because no label switching possible)
a_ms_bias_noL[, , , , , 3][lblsw_bin_hom[, , , , 1, ]] <- NA # for K=3
### Aggregate
a_ms_bias_agg_it <- apply(a_ms_bias, c(1:4, 6), mean, na.rm=TRUE)
a_ms_bias_noL_agg_it <- apply(a_ms_bias_noL, c(1:4, 6), mean, na.rm=TRUE)
### Apply 90% rule
a_ms_bias_noL_agg_it[, , , , 2][lblsw_bin_hom_90sat_2st] <- NA
a_ms_bias_noL_agg_it[, , , , 3][lblsw_bin_hom_90sat] <- NA
### Aggregate across p and Nsubj
a_ms_bias_agg_it_p_Nsubj <- apply(a_ms_bias_agg_it, c(1, 3, 5), mean, na.rm=TRUE)
a_ms_bias_agg_it_p_Nsubj <- apply(a_ms_bias_noL_agg_it, c(1, 3, 5), mean, na.rm=TRUE)
### Cmb in list
l_ms_bias_agg_it_p_Nsubj_AICc <- list(a_ms_bias_agg_it_p_Nsubj,
                                      a_ms_bias_agg_it_p_Nsubj)



# --------------------------------------------------------
# ---------- Preprocessing: Model Selection: AIC ---------
# --------------------------------------------------------

# ----- Collapse into one Array -----

# Each flat table into array
l_ms <- list()
for(k in 1:3) {
  l_ms[[k]] <- xtabs(AIC ~ K_candidate + KL_div + n_dep + n_t + n + sim_iteration, data = l_res[[k]]) # Doesn't work because I have different dependent variables
  l_ms[[k]][l_ms[[k]] == 0] <- NA # because xtabs outrageously inserts 0s when cmb are missing
} # end for: k
D <- length(dim(l_ms[[k]]))

# ----- Make Predictions -----
l_ms_pred <- list()
for(k in 1:3) l_ms_pred[[k]] <- apply(l_ms[[k]] , 2:D, function(x) {
  if(any(is.na(x))) NA else which.min(x)
})

# ----- Compute Accuracy -----
l_ms_acc <- list()
for(k in 1:3) l_ms_acc[[k]] <- l_ms_pred[[k]] == k
# Combine in our array
a_ms_acc <- abind(l_ms_acc[[1]],
                  l_ms_acc[[2]],
                  l_ms_acc[[3]],
                  along = length(c(3, dim(l_ms_acc[[k]]))))

### Subset on "no Label Switching"
a_ms_acc_noL <- a_ms_acc
a_ms_acc_noL[, , , , , 2][lblsw_bin_hom_2st[, , , , 1, ]] <- NA # for K=2 (for K=1 not needed, because no label switching possible)
a_ms_acc_noL[, , , , , 3][lblsw_bin_hom[, , , , 1, ]] <- NA # for K=3
### Aggregate
a_ms_acc_agg_it <- apply(a_ms_acc, c(1:4, 6), mean, na.rm=TRUE)
a_ms_acc_noL_agg_it <- apply(a_ms_acc_noL, c(1:4, 6), mean, na.rm=TRUE)
### Apply 90% rule
a_ms_acc_noL_agg_it[, , , , 2][lblsw_bin_hom_90sat_2st] <- NA
a_ms_acc_noL_agg_it[, , , , 3][lblsw_bin_hom_90sat] <- NA
## Cmb in list
l_ms_acc_noL_agg_it <- list(a_ms_acc_agg_it,
                            a_ms_acc_noL_agg_it)
### Aggregate across p and Nsubj
a_ms_acc_agg_it_p_Nsubj <- apply(a_ms_acc_agg_it, c(1, 3, 5), mean, na.rm=TRUE)
a_ms_acc_noL_agg_it_p_Nsubj <- apply(a_ms_acc_noL_agg_it, c(1, 3, 5), mean, na.rm=TRUE)
### Cmb in list
l_ms_acc_agg_it_p_Nsubj <- list(a_ms_acc_agg_it_p_Nsubj,
                                a_ms_acc_noL_agg_it_p_Nsubj)

# ----- Compute Bias -----
l_ms_bias <- list()
for(k in 1:3) l_ms_bias[[k]] <- l_ms_pred[[k]] - k
# Combine in our array
a_ms_bias <- abind(l_ms_bias[[1]],
                   l_ms_bias[[2]],
                   l_ms_bias[[3]],
                   along = length(c(3, dim(l_ms_bias[[k]]))))

### Subset on "no Label Switching"
a_ms_bias_noL <- a_ms_bias
a_ms_bias_noL[, , , , , 2][lblsw_bin_hom_2st[, , , , 1, ]] <- NA # for K=2 (for K=1 not needed, because no label switching possible)
a_ms_bias_noL[, , , , , 3][lblsw_bin_hom[, , , , 1, ]] <- NA # for K=3
### Aggregate
a_ms_bias_agg_it <- apply(a_ms_bias, c(1:4, 6), mean, na.rm=TRUE)
a_ms_bias_noL_agg_it <- apply(a_ms_bias_noL, c(1:4, 6), mean, na.rm=TRUE)
### Apply 90% rule
a_ms_bias_noL_agg_it[, , , , 2][lblsw_bin_hom_90sat_2st] <- NA
a_ms_bias_noL_agg_it[, , , , 3][lblsw_bin_hom_90sat] <- NA
## Cmb in list
l_ms_bias_noL_agg_it <- list(a_ms_bias_agg_it,
                             a_ms_bias_noL_agg_it)
### Aggregate across p and Nsubj
a_ms_bias_agg_it_p_Nsubj <- apply(a_ms_bias_agg_it, c(1, 3, 5), mean, na.rm=TRUE)
a_ms_bias_agg_it_p_Nsubj <- apply(a_ms_bias_noL_agg_it, c(1, 3, 5), mean, na.rm=TRUE)
### Cmb in list
l_ms_bias_agg_it_p_Nsubj <- list(a_ms_bias_agg_it_p_Nsubj,
                                 a_ms_bias_agg_it_p_Nsubj)


# -------------------------------------------------
# ----- Figure: Model Selection 2x3 Main Text -----
# -------------------------------------------------

cols_MS <- qualitative_hcl(n = 3, palette = "Set2")
# plot(1:3, col=cols_MS, pch=20, cex=8)

LS <- c("All", "noL")

# Loop over: AIC and AICc
IC <- c("AIC", "AICc")
l_ms_acc_agg_it_p_Nsubj_LCs <- list(l_ms_acc_agg_it_p_Nsubj,
                                    l_ms_acc_agg_it_p_Nsubj_AICc)
l_ms_bias_agg_it_p_Nsubj_LCs <- list(l_ms_bias_agg_it_p_Nsubj,
                                     l_ms_bias_agg_it_p_Nsubj_AICc)
for(lc in 1:2) {
  # Loop over: with/without label switching
  for(v in 1:2) {
    sc <- 0.72
    pdf(paste0("Figures/Fig_Res_ModelSelect_Accuracy_Bias_Agg_", IC[lc], "_", LS[v], ".pdf"), width = 8*sc, height = 8*sc*(2/3)*1)

    # ----- Define Layout -----
    lmat <- rbind(c(0, 1:3),
                  c(4, 6:8),
                  c(5, 9:11))

    lo <- layout(mat = lmat,
                 widths = c(.1, 1, 1, 1),
                 heights = c(.2, 1, 1, 1))
    # layout.show(lo)

    # ----- Plot Labels -----
    # Rows: True K
    plotLabel("           M = 1")
    plotLabel("           M = 2")
    plotLabel("           M = 3")

    # Rows: True K
    plotLabel("           Accuracy", srt=90)
    plotLabel("           Bias", srt=90)

    # ----- Accuracy -----
    for(k in 1:3) {
      # par(mar=c(3.5,3,0,0.5))
      par(mar=c(4,3.75,0,0.5))
      plot.new()
      plot.window(xlim=c(1, 5), ylim=c(0,1))
      grid()
      axis(1, labels=Ntvar, at=1:5)
      axis(2, las=2)
      for(kld in 1:3)  lines(1:5, l_ms_acc_agg_it_p_Nsubj_LCs[[lc]][[v]][kld, , k], col=cols_MS[kld], lwd=2, type="o", lty=kld)

      # Legend
      if(k==1) legend("center", legend=paste0("KLD = ", KLDvar),
                      text.col = cols_MS, bty="n", cex=1.2)
    }

    # ----- Bias -----
    for(k in 1:3) {
      # par(mar=c(3.5,3,0,0.5))
      par(mar=c(4,3.75,0,0.5))
      plot.new()
      plot.window(xlim=c(1, 5), ylim=c(-3,3))
      grid()
      axis(1, labels=Ntvar, at=1:5)
      axis(2, las=2)
      abline(h=0, lty=2, col="lightgrey")
      title(xlab=expression(N[t]), line=2.4)
      for(kld in 1:3)  lines(1:5, l_ms_bias_agg_it_p_Nsubj_LCs[[lc]][[v]][kld, , k], col=cols_MS[kld], lwd=2, type="o", lty=kld)
    }

    dev.off()
  } # end for: version
} # end for: lc



# -------------------------------------------------
# ----- Figure: Model Sel: Accuracy 3x3 Appendix --
# -------------------------------------------------

# Pick 3 colors (for # dependent variables)
cols <- brewer.pal(4, "Set2")

# Loop over: with/without label switching
for(v in 1:2) {
  sc <- 0.72
  pdf(paste0("Figures/Fig_Res_ModelSelect_Accuracy_", LS[v], ".pdf"), width = 8*sc, height = 8*sc*0.95)

  # ----- Define Layout -----
  lmat <- rbind(c(0, 1:3),
                c(4, 7:9),
                c(5, 10:12),
                c(6, 13:15))

  lo <- layout(mat = lmat,
               widths = c(.15, 1, 1, 1),
               heights = c(.15, 1, 1, 1))
  # layout.show(lo)

  # ----- Plot Labels -----
  # Cols: KLDs
  plotLabel(expression("           D"["KL"]*" = 3"))
  plotLabel(expression("           D"["KL"]*" = 5"))
  plotLabel(expression("           D"["KL"]*" = 7"))

  # Rows: True K
  plotLabel("           M = 1", srt=90)
  plotLabel("           M = 2", srt=90)
  plotLabel("           M = 3", srt=90)

  # ----- Loop in Data -----
  for(k in 1:3) {
    for(kld in 1:3) {
      par(mar=c(4,3.75,0,0.5))
      plot.new()
      plot.window(xlim=c(1, 5), ylim=c(0,1))
      grid()
      axis(1, labels=Ntvar, at=1:5)
      axis(2, las=2)
      if(kld==1) title(ylab="Accuracy", line=2.4)
      if(k==3) title(xlab=expression(N[t]), line=2.4)
      for(p in 1:3) for(s in 1:4) lines(1:5, l_ms_acc_noL_agg_it[[v]][kld, p, , s, k],
                                        col=cols_MS[p], lwd=1.5, type="l", lty=s)

      # Legend
      if(k==1 & kld==1) legend("center", legend=paste0(rep(c(2,4,8), each=4), " Variables, ", rep(c(15, 30, 60, 120), 3), " Subjecs"),
                               text.col = rep(cols_MS, each=4),
                               col = rep(cols_MS, each=4),
                               bty="n", cex=0.67, lty=rep(1:4, 3))

    }
  }

  dev.off()
} # For: versions


# -------------------------------------------------
# ----- Figure: Model Sel: Bias 3x3 Appendix --
# -------------------------------------------------

# Pick 3 colors (for # dependent variables)
cols <- brewer.pal(4, "Set2")

# Loop over: with/without label switching
for(v in 1:2) {
  sc <- 0.72
  pdf(paste0("Figures/Fig_Res_ModelSelect_Bias_", LS[v], ".pdf"), width = 8*sc, height = 8*sc*0.95)

  # ----- Define Layout -----
  lmat <- rbind(c(0, 1:3),
                c(4, 7:9),
                c(5, 10:12),
                c(6, 13:15))

  lo <- layout(mat = lmat,
               widths = c(.15, 1, 1, 1),
               heights = c(.15, 1, 1, 1))
  # layout.show(lo)

  # ----- Plot Labels -----
  # Cols: KLDs
  plotLabel(expression("           D"["KL"]*" = 3"))
  plotLabel(expression("           D"["KL"]*" = 5"))
  plotLabel(expression("           D"["KL"]*" = 7"))

  # Rows: True K
  plotLabel("           M = 1", srt=90)
  plotLabel("           M = 2", srt=90)
  plotLabel("           M = 3", srt=90)

  # ----- Loop in Data -----
  for(k in 1:3) {
    for(kld in 1:3) {
      par(mar=c(4,3.75,0,0.5))
      plot.new()
      plot.window(xlim=c(1, 5), ylim=c(-3,3))
      grid()
      axis(1, at=1:5, labels=Ntvar)
      axis(2, las=2)
      if(kld==1) title(ylab="Abs Relative Bias", line=2.4)
      if(k==3) title(xlab=expression(N[t]), line=2.4)
      for(p in 1:3) for(s in 1:4) lines(1:5, l_ms_bias_noL_agg_it[[v]][kld, p, , s, k],
                                        col=cols_MS[p], lwd=1.5, type="l", lty=s)

      # Legend
      if(k==1 & kld==1) legend("center", legend=paste0(rep(c(2,4,8), each=4), " Variables, ", rep(c(15, 30, 60, 120), 3), " Subjecs"),
                               text.col = rep(cols_MS, each=4),
                               col = rep(cols_MS, each=4),
                               bty="n", cex=0.67, lty=rep(1:4, 3))

    }
  }

  dev.off()
} # for: versions


# ---------------------------------------------------------------
# ---------- Load Results: Emission Distributions ---------------
# ---------------------------------------------------------------

res_emiss <- readRDS("result_tables/Performance_emission_3st.RDS")
# NOTE: this actually contains the data for all K=1,2,3

head(res_emiss)
dim(res_emiss)


# --------------------------------------------------------
# ---------- Preprocessing: Emission Distributions -------
# --------------------------------------------------------

# ----- Emission Mean: Bias -----
emiss_mean_true <- xtabs(mean_true ~ KL_div + n_dep + n_t + n + dep + sim_iteration + k, data = res_emiss)
emiss_mean_true[emiss_mean_true==0] <- NA
emiss_mean_est <- xtabs(mean_hat ~ KL_div + n_dep + n_t + n + dep + sim_iteration + k, data = res_emiss)
emiss_mean_est[emiss_mean_est==0] <- NA
dim(emiss_mean_est)
# Calculate Bias
emiss_mean_bias <- abs( (emiss_mean_est - emiss_mean_true )  / emiss_mean_true )
# Subset on: "no label switching"
emiss_mean_bias_noL <- emiss_mean_bias
emiss_mean_bias_noL[lblsw_bin_hom] <- NA
# Aggregate
emiss_mean_bias_agg_iter <- apply(emiss_mean_bias, 1:4, median, na.rm=TRUE)
emiss_mean_bias_noL_agg_iter <- apply(emiss_mean_bias_noL, 1:4, median, na.rm=TRUE)
# Apply 90%-rule
emiss_mean_bias_noL_agg_iter[lblsw_bin_hom_90sat] <- NA
# Combine to list
l_emiss_mean_bias_agg_iter <- list(emiss_mean_bias_agg_iter,
                                   emiss_mean_bias_noL_agg_iter)


# ----- Emission Mean: Coverage -----
emiss_mean_est_lower <- xtabs(mean_95_lower ~ KL_div + n_dep + n_t + n + dep + sim_iteration + k, data = res_emiss)
emiss_mean_est_lower[emiss_mean_est_lower==0] <- NA
emiss_mean_est_upper <- xtabs(mean_95_upper ~ KL_div + n_dep + n_t + n + dep + sim_iteration + k, data = res_emiss)
emiss_mean_est_upper[emiss_mean_est_upper==0] <- NA
dim(emiss_mean_est_upper)
# Again subset on K=3
emiss_mean_cvrg <- (emiss_mean_est_lower < emiss_mean_true) & (emiss_mean_true < emiss_mean_est_upper)
dim(emiss_mean_cvrg)
# Subset on: "no label switching"
emiss_mean_cvrg_noL <- emiss_mean_cvrg
emiss_mean_cvrg_noL[lblsw_bin_hom] <- NA
# Aggregate
emiss_mean_cvrg_agg_iter <- apply(emiss_mean_cvrg, 1:4, mean, na.rm=TRUE)
emiss_mean_cvrg_noL_agg_iter <- apply(emiss_mean_cvrg_noL, 1:4, mean, na.rm=TRUE)
# Apply 90%-rule
emiss_mean_cvrg_noL_agg_iter[lblsw_bin_hom_90sat] <- NA
# Combine to list
l_emiss_mean_cvrg_agg_iter <- list(emiss_mean_cvrg_agg_iter,
                                   emiss_mean_cvrg_noL_agg_iter)

# ----- Emission SD: Bias -----
emiss_SD_true <- xtabs(SD_true ~ KL_div + n_dep + n_t + n + dep + sim_iteration + k, data = res_emiss)
emiss_SD_true[emiss_SD_true==0] <- NA
emiss_SD_est <- xtabs(SD_hat ~ KL_div + n_dep + n_t + n + dep + sim_iteration + k, data = res_emiss)
emiss_SD_est[emiss_SD_est==0] <- NA
# Calculate Bias
emiss_SD_bias <-  (emiss_SD_est - emiss_SD_true )  / emiss_SD_true
# Subset on: "no label switching"
emiss_SD_bias_noL <- emiss_SD_bias
emiss_SD_bias_noL[lblsw_bin_hom] <- NA
# Aggregate
emiss_SD_bias_agg_iter <- apply(emiss_SD_bias, 1:4, median, na.rm=TRUE)
emiss_SD_bias_noL_agg_iter <- apply(emiss_SD_bias_noL, 1:4, median, na.rm=TRUE)
# Apply 90%-rule
emiss_SD_bias_noL_agg_iter[lblsw_bin_hom_90sat] <- NA
# Combine to list
l_emiss_SD_bias_agg_iter <- list(emiss_SD_bias_agg_iter,
                                 emiss_SD_bias_noL_agg_iter)

# ----- Emission SD: Coverage -----
emiss_SD_est_lower <- xtabs(SD_95_lower ~ KL_div + n_dep + n_t + n + dep + sim_iteration + k, data = res_emiss)
emiss_SD_est_lower[emiss_SD_est_lower==0] <- NA
emiss_SD_est_upper <- xtabs(SD_95_upper ~ KL_div + n_dep + n_t + n + dep + sim_iteration + k, data = res_emiss)
emiss_SD_est_upper[emiss_SD_est_upper==0] <- NA
# Calculate Coverage
emiss_SD_cvrg <- (emiss_SD_est_lower < emiss_SD_true) & (emiss_SD_true < emiss_SD_est_upper)
# Subset on: "no label switching"
emiss_SD_cvrg_noL <- emiss_SD_cvrg
emiss_SD_cvrg_noL[lblsw_bin_hom] <- NA
# Aggregate
emiss_SD_cvrg_agg_iter <- apply(emiss_SD_cvrg, 1:4, mean, na.rm=TRUE)
emiss_SD_cvrg_noL_agg_iter <- apply(emiss_SD_cvrg_noL, 1:4, mean, na.rm=TRUE)
# Apply 90%-rule
emiss_SD_cvrg_noL_agg_iter[lblsw_bin_hom_90sat] <- NA
# Combine to list
l_emiss_SD_cvrg_agg_iter<- list(emiss_SD_cvrg_agg_iter,
                                emiss_SD_cvrg_noL_agg_iter)

# ----- Emission SD: Prediction / SD of Estimate -----
# Subset on: "no label switching"
emiss_SD_est_noL <- emiss_SD_est
emiss_SD_est_noL[lblsw_bin_hom] <- NA
# Aggregate
# Average over dependent variables and iteration
emiss_SD_est_agg_iter <- apply(emiss_SD_est, 1:4, sd, na.rm=TRUE)
emiss_SD_est_noL_agg_iter <- apply(emiss_SD_est_noL, 1:4, sd, na.rm=TRUE)
# Apply 90%-rule
emiss_SD_est_noL_agg_iter[lblsw_bin_hom_90sat] <- NA
# Cmb to list
l_emiss_SD_est_agg_iter <- list(emiss_SD_est_agg_iter,
                                emiss_SD_est_noL_agg_iter)


# --------------------------------------------------------
# ---------- Figure: Mean: Bias --------------------------
# --------------------------------------------------------
# median absolute relative bias

# Using "colorspace" package
cols <- qualitative_hcl(n = 4, palette = "Dark3")

LS <- c("All", "noL")

# ------ All Three KLDs -----

# Loop over: with/without label switching
for(v in 1:2) {
  sc <- 0.85
  pdf(paste0("Figures/Fig_Res_Estimation_Mean_AbsRelBias_", LS[v] ,".pdf"), width = 8*sc, height = 8*sc)

  PlotSimEmiss(object = l_emiss_mean_bias_agg_iter[[v]],
               ylim = c(0, .8),
               h_ab=0.1,
               leg_pos="topright")

  dev.off()
}

# ------ Version without KLD=1 -----

# Loop over: with/without label switching
for(v in 1:2) {
  sc <- 0.72
  pdf(paste0("Figures/Fig_Res_Estimation_Mean_AbsRelBias_KLD23_", LS[v] ,".pdf"), width = 8*sc, height = 8*sc*(2/3)*0.95)

  PlotSimEmiss_KLD23(object = l_emiss_mean_bias_agg_iter[[v]],
                     ylim = c(0, .15),
                     h_ab=0.1,
                     ylab="Absolute Relative Bias",
                     leg_pos=c(1, .148))

  dev.off()
}


# --------------------------------------------------------
# ---------- Figure: Mean: Coverage ----------------------
# --------------------------------------------------------
# How often is true within the CI

# ------ All Three KLDs -----

# Loop over: with/without label switching
for(v in 1:2) {
  sc <- 0.85
  pdf(paste0("Figures/Fig_Res_Estimation_Mean_Coverage_", LS[v] ,".pdf"), width = 8*sc, height = 8*sc)

  PlotSimEmiss(object = l_emiss_mean_cvrg_agg_iter[[v]],
               ylim = c(0, 1),
               h_ab=0.95,
               leg_pos="bottomright")

  dev.off()
}

# ------ Version without KLD=1 -----

for(v in 1:2) {
  sc <- 0.72
  pdf(paste0("Figures/Fig_Res_Estimation_Mean_Coverage__KLD23_", LS[v] ,".pdf"), width = 8*sc, height = 8*sc*(2/3)*0.95)

  PlotSimEmiss_KLD23(object = l_emiss_mean_cvrg_agg_iter[[v]],
                     ylim = c(0, 1),
                     h_ab=0.95,
                     leg_pos="right",
                     leg_box = c(1,2),
                     ylab="Coverage")

  dev.off()
}





# --------------------------------------------------------
# ---------- Figure: SD: Bias ----------------------------
# --------------------------------------------------------
# median relative bias (bc always the same)


# ------ All Three KLDs -----

for(v in 1:2) {
  sc <- 0.85
  pdf(paste0("Figures/Fig_Res_Estimation_SD_RelBias_", LS[v] ,".pdf"), width = 8*sc, height = 8*sc)

  PlotSimEmiss(object = l_emiss_SD_bias_agg_iter[[v]],
               ylim = c(0, 0.4),
               h_ab=0.1,
               leg_pos="topright")

  dev.off()
}

# ------ Version without KLD=1 -----

# Loop over: with/without label switching
for(v in 1:2) {
  sc <- 0.72
  pdf(paste0("Figures/Fig_Res_Estimation_SD_RelBias_KLD23_", LS[v] ,".pdf"), width = 8*sc, height = 8*sc*(2/3)*0.95)

  PlotSimEmiss_KLD23(object = l_emiss_SD_bias_agg_iter[[v]],
                     ylim = c(0, 0.15),
                     leg_pos="bottom",
                     leg_box=c(2,1),
                     ylab="Relative Bias")

  dev.off()
}


# --------------------------------------------------------
# ---------- Figure: SD: Precision -----------------------
# --------------------------------------------------------
# SD of estimates


# ------ All Three KLDs -----
# Loop over: with/without label switching
for(v in 1:2) {
  sc <- 0.85
  pdf(paste0("Figures/Fig_Res_Estimation_SD_Precision_", LS[v] ,".pdf"), width = 8*sc, height = 8*sc)

  PlotSimEmiss(object = l_emiss_SD_est_agg_iter[[v]],
               ylim = c(0, 0.2),
               h_ab=NULL,
               leg_pos="topright")

  dev.off()
}

# ------ Version without KLD=1 -----
# Loop over: with/without label switching
for(v in 1:2) {
  sc <- 0.72
  pdf(paste0("Figures/Fig_Res_Estimation_SD_Precision_KLD23_", LS[v] ,".pdf"), width = 8*sc, height = 8*sc*(2/3)*0.95)

  PlotSimEmiss_KLD23(object = l_emiss_SD_est_agg_iter[[v]],
                     ylim = c(0, 0.15),
                     h_ab=NULL,
                     leg_pos="topright", leg_box = c(2,2),
                     ylab="Precision")

  dev.off()
}



# --------------------------------------------------------
# ---------- Figure: SD: Coverage ------------------------
# --------------------------------------------------------
# How often is true within the CI

# ------ All Three KLDs -----
for(v in 1:2) {
  sc <- 0.85
  pdf(paste0("Figures/Fig_Res_Estimation_SD_Coverage_", LS[v] ,".pdf"), width = 8*sc, height = 8*sc)

  PlotSimEmiss(object = l_emiss_SD_cvrg_agg_iter[[v]],
               ylim = c(0, 1),
               h_ab=0.95,
               leg_pos="bottomright")

  dev.off()
}

# ------ Version without KLD=1 -----
for(v in 1:2) {
  sc <- 0.72
  pdf(paste0("Figures/Fig_Res_Estimation_SD_Coverage__KLD23_", LS[v] ,".pdf"), width = 8*sc, height = 8*sc*(2/3)*0.95)

  PlotSimEmiss_KLD23(object = l_emiss_SD_cvrg_agg_iter[[v]],
                     ylim = c(0, 1),
                     h_ab=0.95,
                     leg_pos="right",
                     leg_box = c(1,2),
                     ylab="Coverage")

  dev.off()
}



# ---------------------------------------------------------------
# ---------- Load Results: Transition Probabilities -------------
# ---------------------------------------------------------------

res_gam <- readRDS("result_tables/Performance_gamma_3st.RDS")


# --------------------------------------------------------
# ---------- Preprocessing: Transition Probabilities -----
# --------------------------------------------------------
#
# # ----- Transition: Bias -----
# # Issue because data structure is different than I imagined, so will do manually
# a_trans_true <- a_trans_est <- array(NA, dim=c(3, 3, 3, 3, 5, 4, 99))
# for(kld in 1:3) for(p in 1:3) for(nt in 1:5) for (ns in 1:4) for(fi in 1:3) for(ti in 1:3) {
#   a_trans_true[fi, ti, kld, p, nt, ns, ] <- try(res_gam$gamma_ij_true[res_gam$n==c(15, 30, 60, 120)[ns] &
#                                                                         res_gam$n_t==c(50, 100, 200, 400, 800)[n] &
#                                                                         res_gam$n_dep==c(2,4,8)[p] &
#                                                                         res_gam$KL_div == c(3,5,7)[kld] &
#                                                                         res_gam$from_state_i == fi &
#                                                                         res_gam$to_state_j == ti])
# } # end: for stack

# ----- Get this into Arrays ------
trans_true <- xtabs(gamma_ij_true ~ KL_div + n_dep + n_t + n + from_state_i + to_state_j + sim_iteration + true_m, data = res_gam)
trans_true[trans_true==0] <- NA
trans_est <- xtabs(gamma_ij_hat ~ KL_div + n_dep + n_t + n + from_state_i + to_state_j + sim_iteration + true_m, data = res_gam)
trans_est[trans_est==0] <- NA
# Dim
dim(trans_true)

# ----- Adapt Label Switching array to match transition probability array -----
lblsw_bin_tnsm <- array(NA, dim=c(3, 3, 5, 4, 3, 3, 100, 1))
for(kld in 1:3) for(p in 1:3) for(nt in 1:5) for (ns in 1:4) for(fi in 1:3) for(ti in 1:3) for(i in 1:100) {
  # We use the same aggregation rule as above, requiring TRUE for
  lblsw_bin_tnsm[kld, p, nt, ns, fi, ti, i, 1] <- lblsw_bin_hom[kld, p, nt, ns, 1, i] # We take first dep variable, because above we already homogenized this
}
dim(lblsw_bin_tnsm)

# ----- Same for the 90% rule array -----
rule90_tnsm <- array(NA, dim=c(3, 3, 5, 4, 3, 3))
for(kld in 1:3) for(p in 1:3) for(nt in 1:5) for (ns in 1:4) for(fi in 1:3) for(ti in 1:3) {
  # We use the same aggregation rule as above, requiring TRUE for
  rule90_tnsm[kld, p, nt, ns, fi, ti] <- lblsw_bin_hom_90sat[kld, p, nt, ns] # We take first dep variable, because above we already homogenized this
}
dim(rule90_tnsm)


# ----- Calculate Raw Bias ------
trans_bias_raw <- trans_est - trans_true
# Subset on: "no label switching"
trans_bias_raw_noL <- trans_bias_raw
dim(trans_bias_raw_noL)
trans_bias_raw_noL[lblsw_bin_tnsm] <- NA
# Aggregate
trans_bias_raw_agg_iter <- apply(trans_bias_raw, 1:6, mean, na.rm=TRUE)
trans_bias_raw_noL_agg_iter <- apply(trans_bias_raw_noL, 1:6, mean, na.rm=TRUE)
# Apply 90%-rule
trans_bias_raw_noL_agg_iter[rule90_tnsm] <- NA
# Cmb in list
l_trans_bias_raw_agg_iter <- list(trans_bias_raw_agg_iter,
                                  trans_bias_raw_noL_agg_iter)

# ----- Calculate Absolute Bias ------
trans_bias_abs <- abs(trans_est - trans_true)
# Subset on: "no label switching"
trans_bias_abs_noL <- trans_bias_abs
trans_bias_abs_noL[lblsw_bin_tnsm] <- NA
# Aggregate
trans_bias_abs_agg_iter <- apply(trans_bias_abs, 1:6, mean, na.rm=TRUE)
trans_bias_abs_noL_agg_iter <- apply(trans_bias_abs_noL, 1:6, mean, na.rm=TRUE)
# Apply 90%-rule
trans_bias_abs_noL_agg_iter[rule90_tnsm] <- NA
# Cmb in list
l_trans_bias_abs_agg_iter <- list(trans_bias_abs_agg_iter,
                                  trans_bias_abs_noL_agg_iter)

# ----- Transition: Coverage -----
trans_est_lower <- xtabs(gamma_ij_95_lower ~ KL_div + n_dep + n_t + n + from_state_i + to_state_j + sim_iteration + true_m, data = res_gam)
trans_est_lower[trans_est_lower==0] <- NA
trans_est_upper <- xtabs(gamma_ij_95_upper ~ KL_div + n_dep + n_t + n + from_state_i + to_state_j + sim_iteration + true_m, data = res_gam)
trans_est_upper[trans_est_upper==0] <- NA
# Get Coverage
trans_est_cvrg <- (trans_est_lower < trans_true) & (trans_true < trans_est_upper)
dim(trans_est_cvrg)
# Subset on: "no label switching"
trans_est_noL_cvrg <- trans_est_cvrg
trans_est_noL_cvrg[lblsw_bin_tnsm] <- NA
# Average over dependent variables and iteration
emiss_mean_cvrg_agg_iter <- apply(trans_est_cvrg, 1:6, mean, na.rm=TRUE)
emiss_mean_cvrg_noL_agg_iter <- apply(trans_est_noL_cvrg, 1:6, mean, na.rm=TRUE)
dim(emiss_mean_cvrg_agg_iter)
# Apply 90%-rule
emiss_mean_cvrg_noL_agg_iter[rule90_tnsm] <- NA
# Cmb in list
l_emiss_mean_cvrg_agg_iter <- list(emiss_mean_cvrg_agg_iter,
                                   emiss_mean_cvrg_noL_agg_iter)

# ----- Transition SD: Prediction / SD of Estimate -----
trans_est_forSD <- trans_est
# Subset on: "no label switching"
trans_est_forSD_noL <- trans_est_forSD
trans_est_forSD_noL[lblsw_bin_tnsm] <- NA
# Compute SD over dependent variables and iteration
trans_est_forSD_agg_iter <- apply(trans_est_forSD, 1:6, sd, na.rm=TRUE)
trans_est_forSD_noL_agg_iter <- apply(trans_est_forSD_noL, 1:6, sd, na.rm=TRUE)
# Apply 90%-rule
trans_est_forSD_noL_agg_iter[rule90_tnsm] <- NA
# Cmb in list
l_trans_est_forSD_agg_iter <- list(trans_est_forSD_agg_iter,
                                   trans_est_forSD_noL_agg_iter)


# # Average over dependent variables and iteration
# emiss_SD_Precision <- apply(emiss_SD_est, 1:4, function(x) sd(x, na.rm=TRUE))



# --------------------------------------------------------
# ---------- Figure: Transition: Bias Raw ----------------
# --------------------------------------------------------
# Plot Separately: diagonal and off-diagonal

cols <- qualitative_hcl(n = 4, palette = "Dark3")

# ------ All Three KLDs -----
for(v in 1:2) {
  sc <- 0.85
  pdf(paste0("Figures/Fig_Res_Transition_Bias_Raw_k=3_", LS[v] ,".pdf"), width = 8*sc, height = 8*sc)

  PlotSimTrans(object = l_trans_bias_raw_agg_iter[[v]],
               ylim = c(-0.2, 0.2),
               h_ab=c(-0.1, 0.1),
               leg_pos="topright",
               leg_pos2 = "bottomleft")

  dev.off()
}


# ------ Version without KLD=1 -----
for(v in 1:2) {
  sc <- 0.72
  pdf(paste0("Figures/Fig_Res_Transition_Bias_Raw_KLD23_k=3_", LS[v] ,".pdf"), width = 8*sc, height = 8*sc*(2/3))

  PlotSimTrans_KLD23(object = l_trans_bias_raw_agg_iter[[v]],
                     ylim = c(-0.05, 0.05),
                     # h_ab=c(-0.1, 0.1),
                     ylab="    Mean Bias",
                     ylab_line =2.9,
                     leg_pos="topright",
                     leg_pos2 = "bottomleft")

  dev.off()
}



# --------------------------------------------------------
# ---------- Figure: Transition: Bias Abs ----------------
# --------------------------------------------------------
# Plot Separately: diagonal and off-diagonal

for(v in 1:2) {
  sc <- 0.85
  pdf(paste0("Figures/Fig_Res_Transition_Bias_Abs_k=3_", LS[v] ,".pdf"), width = 8*sc, height = 8*sc)

  PlotSimTrans(object = l_trans_bias_abs_agg_iter[[v]],
               ylim = c(0, 0.3),
               h_ab=c(0.1),
               leg_pos="topright",
               leg_pos2 = "left")

  dev.off()
}


# --------------------------------------------------------
# ---------- Figure: Transition: Coverage ----------------
# --------------------------------------------------------
# Plot Separately: diagonal and off-diagonal

# ------ All Three KLDs -----
for(v in 1:2) {
  sc <- 0.85
  pdf(paste0("Figures/Fig_Res_Transition_Coverage_k=3_", LS[v] ,".pdf"), width = 8*sc, height = 8*sc)

  PlotSimTrans(object = l_emiss_mean_cvrg_agg_iter[[v]],
               ylim = c(0, 1),
               h_ab=0.95,
               leg_pos="right",
               leg_pos2 = "bottomleft")

  dev.off()
}

# ------ Version without KLD=1 -----
for(v in 1:2) {
  sc <- 0.72
  pdf(paste0("Figures/Fig_Res_Transition_Coverage_KLD23_k=3_", LS[v] ,".pdf"), width = 8*sc, height = 8*sc*(2/3))

  PlotSimTrans_KLD23(object = l_emiss_mean_cvrg_agg_iter[[v]],
                     ylim = c(0.8, 1),
                     ylab_line=2.8,
                     h_ab=0.95,
                     leg_pos="right",
                     leg_pos2 = "bottomleft",
                     ylab="Coverage")

  dev.off()
}




# --------------------------------------------------------
# ---------- Figure: Transition: Precision ---------------
# --------------------------------------------------------
# Plot Separately: diagonal and off-diagonal


# ------ Version without KLD=1 -----

for(v in 1:2) {
  sc <- 0.72
  pdf(paste0("Figures/Fig_Res_Transition_Precision_KLD23_k=3_", LS[v] ,".pdf"), width = 8*sc, height = 8*sc*(2/3))

  PlotSimTrans_KLD23(object = l_trans_est_forSD_agg_iter[[v]],
                     ylim = c(0, 0.04),
                     h_ab=0.95,
                     leg_pos="topright", leg2_box = c(1,2),
                     ylab="Precision",
                     ylab_line=2.8,
                     leg_pos2 = "bottomleft")

  dev.off()
}


# --------------------------------------------------------
# ---------- State Decoding: Load Data & Preprocess ------
# --------------------------------------------------------

res_decod <- readRDS("result_tables/State_decoding_3st.RDS")
head(res_decod)

# ----- State Decoding: Accuracy -----
decod_acc <- xtabs(mean_prop_correct ~ KL_div + n_dep + n_t + n + true_m + sim_iteration, data = res_decod)
decod_acc[decod_acc==0] <- NA
decod_acc_drK <- decod_acc[, , , , 1, ] # Drop useless k-dimension
# Subset on "No label switching"
lblsw_bin_hom_drK <- lblsw_bin_hom[, , , , 1, ] # get matching array
decod_acc_drK_noL <- decod_acc_drK
decod_acc_drK_noL[lblsw_bin_hom_drK] <- NA
# Aggregate
decod_acc_drK_agg_ter <- apply(decod_acc, 1:4, mean, na.rm=TRUE)
decod_acc_drK_noL_agg_ter <- apply(decod_acc_drK_noL, 1:4, mean, na.rm=TRUE)
# Apply 90% Rule
decod_acc_drK_noL_agg_ter[lblsw_bin_hom_90sat] <- NA
# Cmb to List
l_decod_acc_drK_agg_ter <- list(decod_acc_drK_agg_ter,
                                decod_acc_drK_noL_agg_ter)

# ----- State Decoding: Cohen's Kappa -----
decod_kappa <- xtabs(mean_kappa ~ KL_div + n_dep + n_t + n + true_m + sim_iteration, data = res_decod)
decod_kappa[decod_kappa==0] <- NA
decod_kappa_drK <- decod_kappa[, , , , 1, ] # Drop useless k-dimension
# Subset on "No label switching"
lblsw_bin_hom_drK <- lblsw_bin_hom[, , , , 1, ] # get matching array
decod_kappa_drK_noL <- decod_kappa_drK
decod_kappa_drK_noL[lblsw_bin_hom_drK] <- NA
# Aggregate
decod_kappa_drK_agg_ter <- apply(decod_kappa_drK, 1:4, mean, na.rm=TRUE)
decod_kappa_drK_noL_agg_ter <- apply(decod_kappa_drK_noL, 1:4, mean, na.rm=TRUE)
# Apply 90% Rule
decod_kappa_drK_noL_agg_ter[lblsw_bin_hom_90sat] <- NA
# Cmb to List
l_decod_kappa_drK_agg_ter <- list(decod_kappa_drK_agg_ter,
                                  decod_kappa_drK_noL_agg_ter)


# --------------------------------------------------------
# ---------- Figure: State Decoding: Accuracy ------------
# --------------------------------------------------------

# ------ All Three KLDs -----

for(v in 1:2) {
  sc <- 0.85
  pdf(paste0("Figures/Fig_Res_Decoding_k=3_Accuracy_", LS[v] ,".pdf"), width = 8*sc, height = 8*sc)

  PlotSimEmiss(object = l_decod_acc_drK_agg_ter[[v]],
               ylim = c(0, 1),
               ylab = "Accuracy",
               h_ab=NULL,
               leg_pos="bottomleft")

  dev.off()
}

# ------ Version without KLD=1 -----

for(v in 1:2) {
  sc <- 0.72
  pdf(paste0("Figures/Fig_Res_Decoding_k=3_Accuracy_KLD23_", LS[v] ,".pdf"), width = 8*sc, height = 8*sc*(2/3))

  PlotSimEmiss_KLD23(object = l_decod_acc_drK_agg_ter[[v]],
                     ylim = c(0.85, 1),
                     ylab = "Accuracy",
                     h_ab=NULL,
                     leg_pos="bottomright",
                     leg_box = c(1,3))

  dev.off()
}



# --------------------------------------------------------
# ---------- Figure: State Decoding: Cohen's Kappa -------
# --------------------------------------------------------


# ------ All Three KLDs -----

for(v in 1:2) {
  sc <- 0.85
  pdf(paste0("Figures/Fig_Res_Decoding_k=3_Kappa_", LS[v] ,".pdf"), width = 8*sc, height = 8*sc)

  PlotSimEmiss(object = l_decod_kappa_drK_agg_ter[[v]],
               ylim = c(0, 1),
               ylab = "Cohen's Kappa",
               h_ab=NULL,
               leg_pos="bottomleft")

  dev.off()
}


# ------ Version without KLD=1 -----

for(v in 1:2) {
  sc <- 0.72
  pdf(paste0("Figures/Fig_Res_Decoding_k=3_Kappa_KLD23_", LS[v] ,".pdf"), width = 8*sc, height = 8*sc*(2/3))

  PlotSimEmiss_KLD23(object = l_decod_kappa_drK_agg_ter[[v]],
                     ylim = c(0.8, 1),
                     ylab = "Cohen's Kappa",
                     h_ab=NULL,
                     leg_pos="top",
                     leg_box = c(1,2))

  dev.off()
}









