# jonashaslbeck@protonmail.com; June 28th, 2024

# --------------------------------------------------------
# ---------- What is happening here? ---------------------
# --------------------------------------------------------

# Helper Functions used throughout the repo


# --------------------------------------------------------
# ---------- Load Packages -------------------------------
# --------------------------------------------------------

library(ggplot2)
library(reshape2)
library(viridis)
library(alluvial)
library(mHMMbayes)

library(coda) # to compute Gelman-Rubin statistic
library(xtable)
library(scales)


# --------------------------------------------------------
# ---------- Compute Gelman-Ruben Statistic --------------
# --------------------------------------------------------

# m <- 3
# model_list <- list(l_Models_1[[m]], l_Models_2[[m]])

f_GR <- function(model_list, digits=2, burnin=500) {

  # Get m
  nChain <- length(model_list) # number of Chains
  m <- model_list[[1]]$input$m

  # Transition Matrix
  iter <- m*(m-1)
  v_RG_gamma <- rep(NA, m*(m-1))
  if(m>1) {
    for(j in 1:iter) {
      l_obj <- lapply(model_list, function(x) mcmc(x$gamma_int_bar[-(1:burnin), j]))
      mcmc_chains <- mcmc.list(l_obj) #make this a list input (must be possible)
      gelman_rubin_stat <- gelman.diag(mcmc_chains)
      v_RG_gamma[j] <- gelman_rubin_stat$psrf[1,1]
    }
  }
  m_RG_gamma <- matrix(v_RG_gamma, m, m-1, byrow = TRUE)

  # Emission Distribution
  RG_emiss <- matrix(NA, n_dep, m)
  for(p in 1:n_dep) {
    for(ms in 1:m) {
      l_obj <- lapply(model_list, function(x) mcmc(x$emiss_mu_bar[[p]][-(1:burnin), ms]))
      mcmc_chains <- mcmc.list(l_obj)
      gelman_rubin_stat <- gelman.diag(mcmc_chains)
      RG_emiss[p, ms] <- gelman_rubin_stat$psrf[1,1]
    }
  }

  # Return
  outlist <- list("m_RG_gamma" = round(m_RG_gamma, digits=digits),
                  "RG_emiss" = round(RG_emiss, digits=digits))
  return(outlist)

} # eoF



# --------------------------------------------------------
# ---------- Plotting Labels in Canvas -------------------
# --------------------------------------------------------

plotLabel <- function(text, cex=1.4, srt=0) {
  par(mar=rep(0,4))
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  text(0.5, 0.5, text, cex=cex, srt=srt)
}


# --------------------------------------------------------
# ---------- Functions to Compute Pseudo-Residuals -------
# --------------------------------------------------------

GetResid <- function(data, # empirical data
                     model, # fitted mHMM model object
                     j, # Desired variable
                     i # Desired subject
) {

  # Aux
  v_subj_id <- unique(data$subj_id)

  # --- Emission Distributions ---
  emiss_subject <- obtain_emiss(model, level = "subject") # For all
  mu_ji <- emiss_subject[[labels[j]]][[i]][, 1] # For specified variable & subject; col=1 for means

  # Get most likely state sequence for each subject
  if(model$input$m > 1) {
    model_m_stateSeq <- suppressMessages(vit_mHMM(model, emotion_mHMM))
    State_i <- model_m_stateSeq$state[model_m_stateSeq$subj==v_subj_id[i]] # state seq for person i
  } else {
    State_i <- rep(1, 240+40*8) # if m=1, we know that state is 1 always
  }

  # Compute predictions for fixed j & i
  Pred_ji <- as.numeric(mu_ji[State_i])

  # Get Empirical data
  X_ji <- data[[labels[j]]][data$subj_id==v_subj_id[i]]

  # Compute residual
  res_ji <- X_ji - Pred_ji

  # Compute MSE
  RMSE_ji <- sqrt(mean(na.omit(res_ji)^2))

  # Return
  outlist <- list("emp" = X_ji,
                  "model" = Pred_ji,
                  "resid" = res_ji,
                  "RMSE" = RMSE_ji)
  return(outlist)

} # eoF


# --------------------------------------------------------
# ---------- Functions Visualize Pseudo-Residuals --------
# --------------------------------------------------------


PlotRes <- function(res_ji, j, i, layout=TRUE, title=TRUE) {

  # Compute AR
  ar_ji <- round(cor(res_ji$resid[-1], res_ji$resid[-(240+8*40)], use="complete.obs"), 3)
  # Fit linear trend
  time <- 1:(240+8*40)
  trend_model <- lm(res_ji$resid ~ time)
  trend_sum <- summary(trend_model)
  trend_coef <- round(trend_sum$coefficients[2, 1], 3)
  trend_pval <- round(trend_sum$coefficients[2, 4], 3)

  # Layout
  if(layout) layout(matrix(1:2, ncol=2), widths = c(1, .5))
  # LinePlot
  par(mar=c(4,3,2,1))
  plot.new()
  plot.window(xlim=c(1, 560), ylim=c(-110, 110))
  axis(1)
  axis(2, las=2)
  abline(h=0, lty=1, col="lightgrey", lwd=2)
  points(res_ji$resid, pch=20)
  abline(trend_model, lwd=2, col="steelblue")
  text(5,-70, paste0("Lin Trend: Slope = ", trend_coef, "; pval = ", trend_pval), col="steelblue", adj=0)
  text(5,-90, paste0("Lag-1 AR = ", ar_ji), col="tomato", adj=0)
  text(4, 80, paste0("RMSE = ", round(res_ji$RMSE, 2)), adj=0)
  if(title) title(main=paste0(labels[j], " [subj = ", v_subj_id[i], "]"), font.main=1)
  # Marginal
  par(mar=c(4,0,2,1))
  hist_data <- hist(res_ji$resid, plot = FALSE, breaks=seq(-110, 110, length=20))
  barplot(hist_data$counts,
          horiz = TRUE,  # Horizontal bars
          names.arg = NULL,
          axes=FALSE)
  x_seq <- seq(-100, 100, length=1000)
  gauss_den <- dnorm(x_seq,
                     mean = mean(res_ji$resid, na.rm = TRUE),
                     sd = sd(res_ji$resid, na.rm = TRUE))
  scaled_den <- (gauss_den/max(gauss_den) ) * max(hist_data$counts)
  lines(scaled_den, seq(0, 24, length=1000), col="grey") # Not entiresure where those 22 come from

} # eoF

# PlotRes(res_ji, j=j, i=i, layout=TRUE)


# --------------------------------------------------------
# ---------- Heatplot for Transition Matrix --------------
# --------------------------------------------------------

plotHeat <- function(phi,
                     k,
                     main="",
                     labels=NULL,
                     las.x=1,
                     cex.axis=1.5,
                     cex.ft = 1.5,
                     cex.val=1) {


  # Save the current par settings
  old_par <- par(no.readonly = TRUE)

  # -- Aux Variables --
  p <- ncol(phi)

  # -- Make color gradient --
  color.gradient <- function(x, colors=c("#E41A1C", "white", "#377EB8"), colsteps=201) {
    return( grDevices::colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
  }
  x <- 1:201
  grad <- color.gradient(x)

  # Make canvas
  par(mar=c(2,2.5,2.5,1)*2)
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0, 1))

  # Auxiliary plotting variables
  sfm <- 1/(p*2)
  seq_mp_x <- seq(0, 1, length=p+1)[-(p+1)] + sfm

  # Plot Axes & Axis labels
  # xy_labels <- paste0("S ", 1:p)
  # Adjust mgp for custom axis labels
  par(mgp = c(3, 0.35, 0))  # Move tick labels closer to the plot

  y_labels <- sapply(p:1, function(i) as.expression(bquote(S[.(i)])))
  x_labels <- sapply(1:p, function(i) as.expression(bquote(S[.(i)])))

  axis(3, labels = x_labels, at=seq_mp_x, cex.axis=cex.axis, tick=FALSE)
  axis(2, labels = y_labels, at=seq_mp_x, las=2, cex.axis=cex.axis, tick=FALSE)
  title(main, font.main=1)

  title(ylab="From", cex.lab=cex.ft, line=2)
  mtext("To", side=3, cex=cex.ft, line=2)

  phi_col <- matrix(NA, p, p)

  # Plot Data
  for(i in 1:p) {
    for(j in 1:p) {

      # Get color
      phi_ij <- phi[p:1, ][j, i]
      if(phi_ij < -1) {
        col_ij <- grad[1]
      } else if(phi_ij > 1 ) {
        col_ij <- grad[201]
      } else {
        col_ij <- grad[phi[p:1, ][j, i] * 100 + 101]
        phi_col[j,i] <- col_ij
      }

      # Plot box
      rect(xleft = seq_mp_x[i]-sfm,
           ybottom = seq_mp_x[j]-sfm,
           xright = seq_mp_x[i]+sfm,
           ytop = seq_mp_x[j]+sfm,
           col = col_ij)
      # Plot text
      text(seq_mp_x[i], seq_mp_x[j], round(phi_ij , 2), cex=cex.val, col="black")
    }
  }

  # Reset par settings back to original
  par(mgp = c(3, 1, 0))


  # Return colors
  return(phi_col)

} # eoF


# ----------------------------------------------------------------------
# ----- Generate GMMs with given K, p and pw KLD -----------------------
# ----------------------------------------------------------------------
# Taken from here: https://github.com/jmbh/OrdinalGMMSim_reparchive/blob/master/aux_functions.R

makeGGM_KLD <- function(K,
                        p,
                        target_KLD,
                        maxIter=20,
                        init_range = c(-1, 1),
                        sig = 0.25,
                        method = "Nelder-Mead",
                        tol = .001,
                        verbose = FALSE) {

  # --- Error Function ---

  fn <- function(par, target_KLD, K, p, sig) {

    # Create data structure for mixtures
    Sigma <- list()
    for(k in 1:K) Sigma[[k]] <- diag(p)*sig

    m_mu <- matrix(par, K, p, byrow=TRUE)
    mu <- list()
    for(k in 1:K) mu[[k]] <- m_mu[k, ]

    # Compute all pairwise KLD
    v_KLD <- rep(NA, K*(K-1)/2)
    count <- 1
    for(i in 1:K) {
      for(j in i:K) {
        if(j!=i) {
          v_KLD[count] <-  KLD(mu1 = mu[[i]],
                               mu2 = mu[[j]],
                               sig1 = Sigma[[i]],
                               sig2 = Sigma[[j]])
          count <- count + 1
        }
      }
    }

    # Compute Error
    v_error <- (v_KLD - target_KLD)^2

    return(sum(v_error))

  } # eof

  # --- Call Optim ---

  conv_err <- 20
  counter <- 1
  while(conv_err > tol) {

    n_pars <- K*p # number of mean-parameters
    par <- runif(n_pars, min = init_range[1], max = init_range[2])

    out <- optim(par = par,
                 fn = fn,
                 target_KLD = target_KLD,
                 K = K,
                 p = p,
                 sig = sig,
                 method = method)


    conv_err <- out$value
    if(verbose) print(out$value)
    counter <- counter + 1

    if(counter > maxIter) stop(paste0("No convergence after ", maxIter," tries."))

  }

  # --- Put Parameters in List Structure ---

  # Check results:
  m_mu_optim <- matrix(out$par, K, p, byrow=TRUE)
  mu <- list()
  for(k in 1:K) mu[[k]] <- m_mu_optim[k, ]
  Sigma <- list()
  for(k in 1:K) Sigma[[k]] <- diag(p)*sig

  outlist <- list("mu" = mu,
                  "Sigma" = Sigma,
                  "fin.error" = out$value,
                  "restarts" = counter)

  return(outlist)


} # eoF


# ----------------------------------------------------------------------
# ----- Make Contour plots of bivariate Gaussian Mixtures --------------
# ----------------------------------------------------------------------
# Adapted from: https://github.com/jmbh/OrdinalGMMSim_reparchive/blob/master/aux_functions.R


PlotCont <- function(mu, Sigma, drawlabels=FALSE) {

  # Number of components
  K <- length(mu)

  # Select colors
  n <- ifelse(K<3, 3, K)
  # cols <- brewer.pal(n = n, name = "Set1")[1:K]
  cols <- c("#F15A22", "#85AB8F", "#008EB0")[1:K] # Colors from Emmeke, to match with systems Figure in "Theoretical Background"-section

  # Make Plot area
  plot.new()
  scale_max <- 3
  plot.window(xlim=c(-scale_max, scale_max), ylim=c(-scale_max, scale_max))
  axis(1, seq(-scale_max, scale_max, length=5))
  axis(2, seq(-scale_max, scale_max, length=5), las=2)

  # Plot contours
  x.points <- seq(-scale_max,scale_max,length.out=100)
  y.points <- x.points
  z <- matrix(0,nrow=100,ncol=100)
  for(k in 1:K) {
    mu_k <- mu[[k]]
    sigma_k <- Sigma[[k]]

    for (i in 1:100) {
      for (j in 1:100) {
        z[i,j] <- mvtnorm::dmvnorm(c(x.points[i],y.points[j]),
                                   mean=mu_k,
                                   sigma=sigma_k)
      }
    }

    # browser()

    # Flatten the z-values
    z_flat <- sort(as.vector(z))
    # Calculate the cumulative density
    cumulative_density <- cumsum(z_flat) / sum(z_flat)
    # Set the desired quantiles (e.g., 50%, 90%)
    quantile_levels <- exp(seq(-.2, -5, length=8))
    # round(quantile_levels, 2)

    # Find the z-values corresponding to these quantiles
    contour_levels <- sapply(quantile_levels, function(q) {
      z_flat[which.min(abs(cumulative_density - q))]
    })

    # Now plot using `contour()` from the graphics package
    contour(x.points, y.points, z, levels = round(contour_levels, 3),
            drawlabels = drawlabels, lwd = 2, cex, add=TRUE,  col=cols[k], )

    # contour(x.points,y.points,z, add=T, col=cols[k], drawlabels = drawlabels)

  } # end for: k

} # eoF



# ----------------------------------------------------------------------
# ----- Compute KL-Divergence for Multivariate Gaussian ----------------
# ----------------------------------------------------------------------
# Taken from here: https://github.com/jmbh/OrdinalGMMSim_reparchive/blob/master/aux_functions.R


KLD <- function(mu1, mu2, sig1, sig2) {

  d <- length(mu1)
  sig1 <- matrix(sig1, d, d)
  sig2 <- matrix(sig2, d, d)
  mu1 <- matrix(mu1, nrow=d)
  mu2 <- matrix(mu2, nrow=d)

  kld <- 0.5 * (    log( det(sig2)/det(sig1) )
                    - d
                    + sum(diag(solve(sig2) %*% sig1))
                    + t(mu1 - mu2) %*% solve(sig2) %*% (mu1 - mu2)
  )

  as.numeric(kld)
}



# ---------------------------------------------------------------------
# ---------- Make Sim Res Plots for Emission Parameters ---------------
# ---------------------------------------------------------------------

PlotSimEmiss <- function(object,
                         ylim,
                         ylab = NULL,
                         h_ab = NULL,
                         leg_pos="topright") {


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
  plotLabel(expression("         D"["KL"]*" = 3"))
  plotLabel(expression("         D"["KL"]*" = 5"))
  plotLabel(expression("         D"["KL"]*" = 7"))

  # Rows: True K
  plotLabel("         p = 2", srt=90)
  plotLabel("         p = 4", srt=90)
  plotLabel("         p = 8", srt=90)

  # ----- Loop in Data -----
  for(p in 1:3) {
    for(kld in 1:3) {
      par(mar=c(4,3.75,0,0.5))
      plot.new()
      plot.window(xlim=c(1, 5), ylim=ylim)
      grid()
      axis(1, labels=Ntvar, at=1:5, las=1)
      axis(2, las=2)
      if(p==3) title(xlab = expression(N[t]), line=2.4)
      if(kld==1) title(ylab=ylab, line=2.5)
      abline(h=h_ab, col="grey")
      for(s in 1:4) lines(1:5, object[kld, p, , s],
                          col=cols[s], lwd=2, type="l")

      # Legend
      if(p==1 & kld==3) legend(leg_pos, legend=paste0(c(15, 30, 60, 120), " Subjects"),
                               text.col = cols,
                               col = cols,
                               bty="n", cex=1)

    } #end for: p
  } # end for: kld

} # eoF


# ---------------------------------------------------------------------
# ---------- Make Sim Res Plots for Emission Parameters [KLD=2,3 only] -
# ---------------------------------------------------------------------

PlotSimEmiss_KLD23 <- function(object,
                               ylim,
                               ylab = NULL,
                               h_ab = NULL,
                               leg_pos="topright",
                               leg_box=c(2,2)) {


  # ----- Define Layout -----
  lmat <- rbind(c(0, 1:3),
                c(4, 6:8),
                c(5, 9:11))

  lo <- layout(mat = lmat,
               widths = c(.15, 1, 1, 1),
               heights = c(.15, 1, 1, 1))
  # layout.show(lo)

  # ----- Plot Labels -----
  # Cols: p
  plotLabel("         p = 2")
  plotLabel("         p = 4")
  plotLabel("         p = 8")

  # Rows: KLD
  plotLabel(expression("            D"["KL"]*" = 5"), srt=90)
  plotLabel(expression("            D"["KL"]*" = 7"), srt=90)

  # # KLD=2, Nsubj = 120, p=2, N=100
  # kld <- 2
  # p <- 1
  # s <- 4
  # object[kld, p, , s] # 0.168

  # ----- Loop in Data -----
  for(kld in 2:3) {
    for(p in 1:3) {
      par(mar=c(4,3.75,0,0.5))
      plot.new()
      plot.window(xlim=c(1, 5), ylim=ylim)
      grid()
      axis(1, labels=Ntvar, at=1:5, las=1)
      axis(2, las=2)
      if(kld==3) title(xlab = expression(N[t]), line=2.4)
      if(p==1) title(ylab=ylab, line=2.75)
      abline(h=h_ab, col="grey")
      for(s in 1:4) lines(1:5, object[kld, p, , s],
                          col=cols[s], lwd=2, type="l")

      # Legend
      if(p==leg_box[1] & kld==leg_box[2]) {

        if(class(leg_pos) == "character") {
          legend(leg_pos, legend=paste0(c(15, 30, 60, 120), " Subjects"),
                 text.col = cols,
                 col = cols,
                 bty="n", cex=1)
        } else {
          legend(leg_pos[1], leg_pos[2], legend=paste0(c(15, 30, 60, 120), " Subjects"),
                 text.col = cols,
                 col = cols,
                 bty="n", cex=1)
        }



      }

    } #end for: p
  } # end for: kld

} # eoF




# ---------------------------------------------------------------------
# ---------- Make Plots for Transition Parameters ----------------------
# ---------------------------------------------------------------------

# object <- trans_bias_agg_iter

PlotSimTrans <- function(object, ylim,
                         ylab = NULL,
                         h_ab = NULL,
                         leg_pos = "topright",
                         leg_pos2="topleft") {

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
  plotLabel(expression("D"["KL"]*" = 3"))
  plotLabel(expression("D"["KL"]*" = 5"))
  plotLabel(expression("D"["KL"]*" = 7"))

  # Rows: True K
  plotLabel("p = 2", srt=90)
  plotLabel("p = 4", srt=90)
  plotLabel("p = 8", srt=90)

  # ----- Loop in Data -----
  for(p in 1:3) {
    for(kld in 1:3) {
      par(mar=c(4,3.75,0,0.5))
      plot.new()
      plot.window(xlim=c(1, 5), ylim=ylim)
      grid()
      axis(1, labels=Ntvar, at=1:5, las=1)
      axis(2, las=2)
      abline(h=h_ab, col="grey")

      if(kld==1) title(ylab=ylab, line=2.5)

      for(s in 1:4) {

        v_diag <- v_off_diag <- rep(NA, n_Ntvar)
        for(n in 1:n_Ntvar) {
          trans_fix <-  object[kld, p, n, s, , ]
          v_diag[n] <- mean(diag(trans_fix)) # diagonal
          v_off_diag[n] <- mean(trans_fix[row(trans_fix) != col(trans_fix)]) # off-diagonal
        } # end for: n

        # Plot over N:
        lines(1:5, v_diag, col=cols[s], lwd=2, type="l", lty=1)
        lines(1:5, v_off_diag, col=cols[s], lwd=2, type="l", lty=2)

      } # end for: s

      # Legend
      if(p==1 & kld==3) {
        legend(leg_pos, legend=paste0(c(15, 30, 60, 120), " Subjects"),
               text.col = cols,
               col = cols,
               bty="n", cex=1)
        legend(leg_pos2, c("Self-Transition", "Between-Transition"), lty=1:2, bty="n")
      } # end if: legend

    } #end for: p
  } # end for: kld

} # eoF



# ---------------------------------------------------------------------
# ---------- Make Plots for Transition Parameters ----------------------
# ---------------------------------------------------------------------

# object <- trans_bias_agg_iter

PlotSimTrans_KLD23 <- function(object,
                               ylim,
                               ylab = NULL,
                               ylab_line=2.5,
                               h_ab = NULL,
                               leg_pos = "topright",
                               leg_pos2="topleft", leg2_box=c(1,3)) {

  # ----- Define Layout -----
  lmat <- rbind(c(0, 1:3),
                c(4, 6:8),
                c(5, 9:11))

  lo <- layout(mat = lmat,
               widths = c(.15, 1, 1, 1),
               heights = c(.15, 1, 1, 1))
  # layout.show(lo)

  # ----- Plot Labels -----
  # Cols: p
  plotLabel("         p = 2")
  plotLabel("         p = 4")
  plotLabel("         p = 8")

  # Rows: KLD
  plotLabel(expression("            D"["KL"]*" = 5"), srt=90)
  plotLabel(expression("            D"["KL"]*" = 7"), srt=90)


  # ----- Loop in Data -----
  for(kld in 2:3) {
    for(p in 1:3) {
      par(mar=c(4,3.75,0,0.5))
      plot.new()
      plot.window(xlim=c(1, 5), ylim=ylim)
      axis(1, labels=Ntvar, at=1:5, las=1)
      axis(2, las=2)
      grid()
      abline(h=h_ab, col="grey")
      if(kld==3) title(xlab = expression(N[t]), line=2.4)

      if(p==1) title(ylab=ylab, line=ylab_line)

      for(s in 1:4) {

        v_diag <- v_off_diag <- rep(NA, n_Ntvar)
        for(n in 1:n_Ntvar) {
          trans_fix <-  object[kld, p, n, s, , ]
          v_diag[n] <- mean(diag(trans_fix)) # diagonal
          v_off_diag[n] <- mean(trans_fix[row(trans_fix) != col(trans_fix)]) # off-diagonal
        } # end for: n

        # Plot over N:
        lines(1:5, v_diag, col=cols[s], lwd=2, type="l", lty=1)
        lines(1:5, v_off_diag, col=cols[s], lwd=2, type="l", lty=2)

      } # end for: s

      # Legend
      if(p==1 & kld==3) {
        legend(leg_pos, legend=paste0(c(15, 30, 60, 120), " Subjects"),
               text.col = cols,
               col = cols,
               bty="n", cex=0.9)
      } # end if: legend
      if(p==leg2_box[1] & kld==leg2_box[2]) {
        legend(leg_pos2, c("Self-Transition", "Between-Transition"), lty=1:2, bty="n", cex=0.9)
      }

    } #end for: p
  } # end for: kld

} # eoF







# ---------------------------------------------------------------------
# ---------- XXXXXX --------------------
# ---------------------------------------------------------------------




