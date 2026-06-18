library(mHMMbayes)
library(ggplot2)
library(reshape2)

# Simulstion settings
settings_SF <- read.csv("Sensitivity/CoreF_settings.csv") # 25
settings_SA <- read.csv("Sensitivity/CoreA_settings.csv") # 28
settings_SB <- read.csv("Sensitivity/CoreB_settings.csv") # 28
settings_SH <- read.csv("Sensitivity/CoreH_settings.csv") # 31
true_m  <- 3
n_dep   <- 4

settings_full <- rbind(settings_SF, settings_SH, settings_SA)

# KL_div  <- settings[num_argv, 3]
# n_t     <- settings[num_argv, 2]
# n       <- settings[num_argv, 1]

# mu_sc        <- settings[num_argv, 4]
# mu_K0_sc     <- settings[num_argv, 5]
# alph_unif    <- settings[num_argv, 6]
# alph_K0_sc   <- settings[num_argv, 7]
# var_tau_sig  <- settings[num_argv, 8]

# n_sim   <- 100
# J       <- 2000
# burn_in <- 500

## Reading in output files ####
# 3 states
# runs on fast core
id_F <- c(1:25)
filesF <- paste0("out_n", settings_SF[id_F, 1], "_nt", settings_SF[id_F, 2],"_KLD", settings_SF[id_F, 3], 
                   "_mu_sc", settings_SF[id_F, 4], "_mu_K0", settings_SF[id_F, 5], "_tans", settings_SF[id_F, 6], 
                 "_trans_K0", settings_SF[id_F, 7], "_varTauSig", settings_SF[id_F, 8]
)
outF <- vector("list", length(filesF))
names(outF) <- filesF

# runs on heavy core
id_H <- c(1:31)
filesH <- paste0("out_n", settings_SH[id_H, 1], "_nt", settings_SH[id_H, 2],"_KLD", settings_SH[id_H, 3], 
                 "_mu_sc", settings_SH[id_H, 4], "_mu_K0", settings_SH[id_H, 5], "_tans", settings_SH[id_H, 6], 
                 "_trans_K0", settings_SH[id_H, 7], "_varTauSig", settings_SH[id_H, 8]
)
outH <- vector("list", length(filesH))
names(outH) <- filesH

# runs splitted over 2 cores A and B
id_AB <- c(1:28)
filesA <-  paste0("out_n", settings_SA[id_AB, 1], "_nt", settings_SA[id_AB, 2],"_KLD", settings_SA[id_AB, 3], 
                  "_mu_sc", settings_SA[id_AB, 4], "_mu_K0", settings_SA[id_AB, 5], "_tans", settings_SA[id_AB, 6], 
                  "_trans_K0", settings_SA[id_AB, 7], "_varTauSig", settings_SA[id_AB, 8], "_A")
filesB <- paste0("out_n", settings_SB[id_AB, 1], "_nt", settings_SB[id_AB, 2],"_KLD", settings_SB[id_AB, 3], 
                 "_mu_sc", settings_SB[id_AB, 4], "_mu_K0", settings_SB[id_AB, 5], "_tans", settings_SB[id_AB, 6], 
                 "_trans_K0", settings_SB[id_AB, 7], "_varTauSig", settings_SB[id_AB, 8], "_B")
outA <- vector("list", length(filesA))
outB <- vector("list", length(filesB))
names(outA) <- filesA
names(outB) <- filesB



group_out_3st_emiss_mean <- vector("list", length = length(c(filesF, filesH, filesA)))
group_out_3st_emiss_sd <- vector("list", length = length(c(filesF, filesH, filesA)))
group_out_3st_gamma_prob <- vector("list", length = length(c(filesF, filesH, filesA)))
group_out_3st_gamma_int <- vector("list", length = length(c(filesF, filesH, filesA)))
means_true <-vector("list", length = length(c(filesF, filesH, filesA)))
selection_out <- vector("list", length = length(c(filesF, filesH, filesA)))
inferred_states <- vector("list", length = length(c(filesF, filesH, filesA)))
true_state_probs <- vector("list", length = length(c(filesF, filesH, filesA)))
true_states <- vector("list", length = length(c(filesF, filesH, filesA)))

names(group_out_3st_emiss_mean) <- c(filesF, filesH, filesA)
names(group_out_3st_emiss_sd) <-c(filesF, filesH, filesA)
names(group_out_3st_gamma_prob) <- c(filesF, filesH, filesA)
names(group_out_3st_gamma_int) <- c(filesF, filesH, filesA)
names(means_true) <- c(filesF, filesH, filesA)
names(selection_out) <- c(filesF, filesH, filesA)
names(inferred_states) <- c(filesF, filesH, filesA)
names(true_state_probs) <- c(filesF, filesH, filesA)
names(true_states) <- c(filesF, filesH, filesA)

for(sc in 1:length(filesF)){
  outF[[sc]] <- readRDS(paste0("Sensitivity/Data/", filesF[sc], ".rds")) 
  selection_out[[sc]]             <- outF[[sc]]$selection_out
  group_out_3st_emiss_mean[[sc]]  <- outF[[sc]]$group_out_3st$emiss_mean
  group_out_3st_emiss_sd[[sc]]    <- outF[[sc]]$group_out_3st$emiss_sd
  group_out_3st_gamma_prob[[sc]]  <- outF[[sc]]$group_out_3st$trans_prob
  group_out_3st_gamma_int[[sc]]   <- outF[[sc]]$group_out_3st$trans_interc
  means_true[[sc]]                <- outF[[sc]]$means_true
  inferred_states[[sc]]           <- outF[[sc]]$inferred_states
  true_state_probs[[sc]]          <- outF[[sc]]$true_state_probs
  true_states[[sc]]               <- outF[[sc]]$true_states
  if(sc > 1){
    outF[[sc]] <- NULL
  }
}

n_sim <- outF[[1]]$input$n_sim
J <- outF[[1]]$input$J
burn_in <- outF[[1]]$input$burn_in

# rm(outF)

l_files <- length(filesF)

for(sc in 1:length(filesH)){
  outH[[sc]] <- readRDS(paste0("Sensitivity/Data/", filesH[sc], ".rds")) 
  selection_out[[sc + l_files]]             <- outH[[sc]]$selection_out
  group_out_3st_emiss_mean[[sc + l_files]]  <- outH[[sc]]$group_out_3st$emiss_mean
  group_out_3st_emiss_sd[[sc + l_files]]    <- outH[[sc]]$group_out_3st$emiss_sd
  group_out_3st_gamma_prob[[sc + l_files]]  <- outH[[sc]]$group_out_3st$trans_prob
  group_out_3st_gamma_int[[sc + l_files]]   <- outH[[sc]]$group_out_3st$trans_interc
  means_true[[sc + l_files]]                <- outH[[sc]]$means_true
  inferred_states[[sc + l_files]]           <- outH[[sc]]$inferred_states
  true_state_probs[[sc + l_files]]          <- outH[[sc]]$true_state_probs
  true_states[[sc + l_files]]               <- outH[[sc]]$true_states
  if(sc > 1){
    outH[[sc]] <- NULL
  }
}

# rm(outH)

l_files_AB <- length(c(filesF, filesH))

for(sc in 1:length(filesA)){
  outA[[sc]] <- readRDS(paste0("Sensitivity/Data/", filesA[sc], ".rds")) 
  outB[[sc]] <- readRDS(paste0("Sensitivity/Data/", filesB[sc], ".rds")) 
  selection_out[[sc + l_files_AB]]                       <- Map(rbind,outA[[sc]]$selection_out, outB[[sc]]$selection_out)
  group_out_3st_emiss_mean[[sc + l_files_AB]]$median     <-  Map(rbind,outA[[sc]]$group_out_3st$emiss_mean$median, outB[[sc]]$group_out_3st$emiss_mean$median)
  group_out_3st_emiss_mean[[sc + l_files_AB]]$quant_2.5  <-  Map(rbind,outA[[sc]]$group_out_3st$emiss_mean$quant_2.5, outB[[sc]]$group_out_3st$emiss_mean$quant_2.5)
  group_out_3st_emiss_mean[[sc + l_files_AB]]$quant_5    <-  Map(rbind,outA[[sc]]$group_out_3st$emiss_mean$quant_5, outB[[sc]]$group_out_3st$emiss_mean$quant_5)
  group_out_3st_emiss_mean[[sc + l_files_AB]]$quant_95   <-  Map(rbind,outA[[sc]]$group_out_3st$emiss_mean$quant_95, outB[[sc]]$group_out_3st$emiss_mean$quant_95)
  group_out_3st_emiss_mean[[sc + l_files_AB]]$quant_97.5 <-  Map(rbind,outA[[sc]]$group_out_3st$emiss_mean$quant_97.5, outB[[sc]]$group_out_3st$emiss_mean$quant_97.5)
  group_out_3st_emiss_sd[[sc + l_files_AB]]$median       <-  Map(rbind,outA[[sc]]$group_out_3st$emiss_sd$median, outB[[sc]]$group_out_3st$emiss_sd$median)
  group_out_3st_emiss_sd[[sc + l_files_AB]]$quant_2.5    <-  Map(rbind,outA[[sc]]$group_out_3st$emiss_sd$quant_2.5, outB[[sc]]$group_out_3st$emiss_sd$quant_2.5)
  group_out_3st_emiss_sd[[sc + l_files_AB]]$quant_5      <-  Map(rbind,outA[[sc]]$group_out_3st$emiss_sd$quant_5, outB[[sc]]$group_out_3st$emiss_sd$quant_5)
  group_out_3st_emiss_sd[[sc + l_files_AB]]$quant_95     <-  Map(rbind,outA[[sc]]$group_out_3st$emiss_sd$quant_95, outB[[sc]]$group_out_3st$emiss_sd$quant_95)
  group_out_3st_emiss_sd[[sc + l_files_AB]]$quant_97.5   <-  Map(rbind,outA[[sc]]$group_out_3st$emiss_sd$quant_97.5, outB[[sc]]$group_out_3st$emiss_sd$quant_97.5)
  group_out_3st_gamma_prob[[sc + l_files_AB]]            <-  Map(rbind,outA[[sc]]$group_out_3st$trans_prob, outB[[sc]]$group_out_3st$trans_prob)
  group_out_3st_gamma_int[[sc + l_files_AB]]             <-  Map(rbind,outA[[sc]]$group_out_3st$trans_interc, outB[[sc]]$group_out_3st$trans_interc)
  means_true[[sc + l_files_AB]]                          <- Map(rbind,outA[[sc]]$means_true, outB[[sc]]$means_true)
  inferred_states[[sc + l_files_AB]]                     <- cbind(outA[[sc]]$inferred_states, outB[[sc]]$inferred_states[,-1])
  true_state_probs[[sc + l_files_AB]]                    <- cbind(outA[[sc]]$true_state_probs, outB[[sc]]$true_state_probs[,-1])
  true_states[[sc + l_files_AB]]                         <- cbind(outA[[sc]]$true_states, outB[[sc]]$true_states[,-1])
  outA[[sc]] <- NULL
  outB[[sc]] <- NULL
}

# rm(outA)
# rm(outB)



## SAVE EXTRACTED FILES ####

extracted_results_Sensitivity <- list(group_out_3st_emiss_mean = group_out_3st_emiss_mean,
                              group_out_3st_emiss_sd = group_out_3st_emiss_sd,
                              group_out_3st_gamma_prob = group_out_3st_gamma_prob,
                              group_out_3st_gamma_int = group_out_3st_gamma_int,
                              means_true = means_true,
                              selection_out = selection_out,
                              inferred_states = inferred_states,
                              true_state_probs = true_state_probs,
                              true_states = true_states, 
                              n_sim = n_sim, 
                              J = J, 
                              burn_in = burn_in)
saveRDS(extracted_results_Sensitivity, file = "Sensitivity/Data/extracted_results_Sensitivity.RDS")

# extracted_results_Sensitivity <- readRDS("Extracted_results/extracted_results_Sensitivity.RDS")


# Model performance ####

### Explaining bias in emission ####
#### Quantifying likely label switching problem occurences ####
rmse <- function(x){
  row_mean <- mean(x)
  rmse <- sqrt(sum((mean(x) - x)^2 / true_m))
  return(rmse)
}

Label_switch_proxy_Sensitivity <- data.frame(sim_iteration = 1:n_sim,
                                     n = rep(settings_full[, 1], each = n_sim * n_dep), 
                                     n_t = rep(settings_full[, 2], each = n_sim * n_dep),
                                     KL_div = rep(settings_full[, 3], each = n_sim * n_dep), 
                                     mu_sc = rep(settings_full[, 4], each = n_sim * n_dep),
                                     mu_K0_sc = rep(settings_full[, 5], each = n_sim * n_dep),
                                     trans_unif = rep(settings_full[, 6], each = n_sim * n_dep), 
                                     trans_K0_sc = rep(settings_full[, 7], each = n_sim * n_dep), 
                                     var_tau_sig = rep(settings_full[, 8], each = n_sim * n_dep),
                                     dep = rep(c(1:n_dep), each = n_sim),
                                     RMSE = NA)


it <- n_sim  * n_dep
cumsum_it <-  c(0,cumsum(rep(n_sim  * n_dep, dim(settings_full)[1])))

for(sc in c(1:84)){ # as we only still have core fast
# for(sc in 1:dim(settings_full)[1]){
  Label_switch_proxy_Sensitivity[(1 + cumsum_it[sc]) : (cumsum_it[sc + 1]), 11] <- unlist(lapply(group_out_3st_emiss_mean[[sc]]$median, apply, 1, rmse))
}

View(Label_switch_proxy_Sensitivity)

saveRDS(Label_switch_proxy_Sensitivity, file = "Sensitivity/Data/Label_switch_proxy_Sensitivity.RDS")
# Label_switch_proxy_Sensitivity <- readRDS("Sensitivity/Data/Label_switch_proxy_Sensitivity.RDS")

aggr_Label_switch_proxy_Sensitivity <- aggregate(Label_switch_proxy_Sensitivity, by = list(Label_switch_proxy_Sensitivity$sim_iteration, Label_switch_proxy_Sensitivity$var_tau_sig,
                                                                                          Label_switch_proxy_Sensitivity$trans_K0_sc, Label_switch_proxy_Sensitivity$trans_unif,
                                                                                          Label_switch_proxy_Sensitivity$mu_K0_sc, Label_switch_proxy_Sensitivity$mu_sc,
                                                                                          Label_switch_proxy_Sensitivity$KL_div,
                                                                                           Label_switch_proxy_Sensitivity$n_t, Label_switch_proxy_Sensitivity$n), FUN = mean)


aggr_Label_switch_proxy_Sensitivity$present <- (aggr_Label_switch_proxy_Sensitivity$RMSE < 0.20 ) * 1

mean(aggr_Label_switch_proxy_Sensitivity$present[aggr_Label_switch_proxy_Sensitivity$KL_div == 7], na.rm = TRUE)

mean(aggr_Label_switch_proxy_Sensitivity$present[aggr_Label_switch_proxy_Sensitivity$KL_div == 5], na.rm = TRUE)


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

ggplot(data = aggr_label2, mapping = aes(x = n_t, y = present, group = interaction(as.factor(n), as.factor(mu_sc), 
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



## Emission distribution ####
Performance_emission_Sensitivity <- data.frame(sim_iteration = 1:n_sim,
                                       n = rep(settings_full[, 1], each = n_sim * true_m * n_dep), 
                                       n_t = rep(settings_full[, 2], each = n_sim * true_m * n_dep), 
                                       KL_div = rep(settings_full[, 3], each = n_sim * true_m * n_dep), 
                                       mu_sc = rep(settings_full[, 4], each = n_sim * true_m * n_dep),
                                       mu_K0_sc = rep(settings_full[, 5], each = n_sim * true_m * n_dep),
                                       trans_unif = rep(settings_full[, 6], each = n_sim * true_m * n_dep), 
                                       trans_K0_sc = rep(settings_full[, 7], each = n_sim * true_m * n_dep), 
                                       var_tau_sig = rep(settings_full[, 8], each = n_sim * true_m * n_dep),
                                       k = rep(1:true_m, each = n_sim),
                                       dep = rep(1:n_dep, each = true_m * n_sim),
                                       mean_hat = NA,
                                       mean_true = NA, 
                                       mean_95_lower = NA, 
                                       mean_95_upper = NA, 
                                       SD_hat = NA, 
                                       SD_true = NA,
                                       SD_95_lower = NA,  
                                       SD_95_upper = NA
)



# group_out_3st_emiss_mean <- extracted_results_Sensitivity$group_out_3st_emiss_mean
# group_out_3st_emiss_sd <- extracted_results_Sensitivity$group_out_3st_emiss_sd
# means_true <- extracted_results_Sensitivity$means_true


cumsum_it <-  c(0,cumsum(rep(n_sim  * n_dep * true_m, dim(settings_full)[1])))

for(sc in c(1:84)){
  Performance_emission_Sensitivity[(1 + cumsum_it[sc]) : (cumsum_it[sc + 1]),12] <- unlist(group_out_3st_emiss_mean[[sc]]$median)
  Performance_emission_Sensitivity[(1 + cumsum_it[sc]) : (cumsum_it[sc + 1]),13] <- unlist(means_true[[sc]])  
  Performance_emission_Sensitivity[(1 + cumsum_it[sc]) : (cumsum_it[sc + 1]),14] <- unlist(group_out_3st_emiss_mean[[sc]]$quant_2.5)
  Performance_emission_Sensitivity[(1 + cumsum_it[sc]) : (cumsum_it[sc + 1]),15] <- unlist(group_out_3st_emiss_mean[[sc]]$quant_97.5)
  Performance_emission_Sensitivity[(1 + cumsum_it[sc]) : (cumsum_it[sc + 1]),16] <- unlist(group_out_3st_emiss_sd[[sc]]$median)
  Performance_emission_Sensitivity[(1 + cumsum_it[sc]) : (cumsum_it[sc + 1]),17] <- 1
  Performance_emission_Sensitivity[(1 + cumsum_it[sc]) : (cumsum_it[sc + 1]),18] <- unlist(group_out_3st_emiss_sd[[sc]]$quant_2.5)
  Performance_emission_Sensitivity[(1 + cumsum_it[sc]) : (cumsum_it[sc + 1]),19] <- unlist(group_out_3st_emiss_sd[[sc]]$quant_97.5)
}




View(Performance_emission_Sensitivity)

saveRDS(Performance_emission_Sensitivity, file = "Sensitivity/Data/Performance_emission_Sensitivity.RDS")
# Performance_emission_Sensitivity <- readRDS(file = "Sensitivity/Data/Performance_emission_Sensitivity.RDS")


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
  
ggplot(data = aggr_Performance_emission_Sensitivity, mapping = aes(x = n_t, y = abs_rel_bias, group = interaction(as.factor(n), as.factor(mu_sc), 
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




# View(aggr_bias_noL_Performance_emission_Sensitivity)

mean(aggr_bias_noL_Performance_emission_Sensitivity$abs_rel_bias)
mean(aggr_bias_noL_Performance_emission_Sensitivity$abs_rel_bias[aggr_bias_noL_Performance_emission_Sensitivity$KL_div == 5])
mean(aggr_bias_noL_Performance_emission_Sensitivity$abs_rel_bias[aggr_bias_noL_Performance_emission_Sensitivity$KL_div == 7])

mean(aggr_bias_noL_Performance_emission_Sensitivity$abs_rel_bias[aggr_bias_noL_Performance_emission_Sensitivity$n == 30])
mean(aggr_bias_noL_Performance_emission_Sensitivity$abs_rel_bias[aggr_bias_noL_Performance_emission_Sensitivity$n == 120])

aggr_bias_noL_Performance_emission_Sensitivity$block[aggr_bias_noL_Performance_emission_Sensitivity$var_tau_sig > 0] <- "BL3_var"
aggr_bias_noL_Performance_emission_Sensitivity$block[aggr_bias_noL_Performance_emission_Sensitivity$trans_unif > 0] <- "BL2_transitions"
aggr_bias_noL_Performance_emission_Sensitivity$block[aggr_bias_noL_Performance_emission_Sensitivity$mu_sc > 0] <- "BL1_emissions"

aggr_bias_noL_Performance_emission_Sensitivity$block <- as.factor(aggr_bias_noL_Performance_emission_Sensitivity$block)  
aggr_bias_noL_Performance_emission_Sensitivity$K0_sc <-aggr_bias_noL_Performance_emission_Sensitivity$mu_K0_sc + aggr_bias_noL_Performance_emission_Sensitivity$trans_K0_sc


ggplot(data = aggr_bias_noL_Performance_emission_Sensitivity, mapping = aes(x = n_t, y = abs_rel_bias, group = interaction(as.factor(n), as.factor(mu_sc), 
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


#### inspecting relative bias for gaussian emission SDs ####
mean(aggr_bias_noL_Performance_emission_Sensitivity$SD_rel_bias)
mean(aggr_bias_noL_Performance_emission_Sensitivity$SD_rel_bias[aggr_bias_noL_Performance_emission_Sensitivity$KL_div == 5])
mean(aggr_bias_noL_Performance_emission_Sensitivity$SD_rel_bias[aggr_bias_noL_Performance_emission_Sensitivity$KL_div == 7])

mean(aggr_bias_noL_Performance_emission_Sensitivity$SD_rel_bias[aggr_bias_noL_Performance_emission_Sensitivity$n_dep == 2])
mean(aggr_bias_noL_Performance_emission_Sensitivity$SD_rel_bias[aggr_bias_noL_Performance_emission_Sensitivity$n_dep == 4])
mean(aggr_bias_noL_Performance_emission_Sensitivity$SD_rel_bias[aggr_bias_noL_Performance_emission_Sensitivity$n_dep == 8])

mean(aggr_bias_noL_Performance_emission_Sensitivity$SD_rel_bias[aggr_bias_noL_Performance_emission_Sensitivity$n == 15])
mean(aggr_bias_noL_Performance_emission_Sensitivity$SD_rel_bias[aggr_bias_noL_Performance_emission_Sensitivity$n == 30])
mean(aggr_bias_noL_Performance_emission_Sensitivity$SD_rel_bias[aggr_bias_noL_Performance_emission_Sensitivity$n == 60])
mean(aggr_bias_noL_Performance_emission_Sensitivity$SD_rel_bias[aggr_bias_noL_Performance_emission_Sensitivity$n == 120])

mean(aggr_bias_noL_Performance_emission_Sensitivity$SD_rel_bias[aggr_bias_noL_Performance_emission_Sensitivity$n_t == 50])
mean(aggr_bias_noL_Performance_emission_Sensitivity$SD_rel_bias[aggr_bias_noL_Performance_emission_Sensitivity$n_t == 100])
mean(aggr_bias_noL_Performance_emission_Sensitivity$SD_rel_bias[aggr_bias_noL_Performance_emission_Sensitivity$n_t == 200])
mean(aggr_bias_noL_Performance_emission_Sensitivity$SD_rel_bias[aggr_bias_noL_Performance_emission_Sensitivity$n_t == 400])
mean(aggr_bias_noL_Performance_emission_Sensitivity$SD_rel_bias[aggr_bias_noL_Performance_emission_Sensitivity$n_t == 800])

mean(aggr_bias_noL_Performance_emission_Sensitivity$SD_rel_bias[aggr_bias_noL_Performance_emission_Sensitivity$KL_div == 5 & 
                                                          aggr_bias_noL_Performance_emission_Sensitivity$n_dep == 2 &
                                                          aggr_bias_noL_Performance_emission_Sensitivity$n != 15])



#### inspecting precision for gaussian emission SDs ####
aggr_precision_noL_Performance_emission_Sensitivity <- aggregate(bias_noL_Performance_emission_Sensitivity, by = list(bias_noL_Performance_emission_Sensitivity$var_tau_sig,
                                                                                                                      bias_noL_Performance_emission_Sensitivity$trans_K0_sc,
                                                                                                                      bias_noL_Performance_emission_Sensitivity$trans_unif,
                                                                                                                      bias_noL_Performance_emission_Sensitivity$mu_K0_sc,
                                                                                                                      bias_noL_Performance_emission_Sensitivity$mu_sc,
                                                                                                                      bias_noL_Performance_emission_Sensitivity$KL_div, 
                                                                                                      bias_noL_Performance_emission_Sensitivity$n_t, 
                                                                                                      bias_noL_Performance_emission_Sensitivity$n), FUN = var)


# aggr_precision_noL_Performance_emission_Sensitivity <- aggr_precision_noL_Performance_emission_Sensitivity[-which(aggr_precision_noL_Performance_emission_Sensitivity$Group.1 == 5 &
#                                                                                                    aggr_precision_noL_Performance_emission_Sensitivity$Group.2 == 2 &
#                                                                                                    aggr_precision_noL_Performance_emission_Sensitivity$Group.4 == 120 &
#                                                                                                    ( aggr_precision_noL_Performance_emission_Sensitivity$Group.3 == 200 |
#                                                                                                        aggr_precision_noL_Performance_emission_Sensitivity$Group.3 == 400)),]

aggr_precision_noL_Performance_emission_Sensitivity$SD_precision <- sqrt(aggr_precision_noL_Performance_emission_Sensitivity$SD_hat)


mean(aggr_precision_noL_Performance_emission_Sensitivity$SD_precision)

max(aggr_precision_noL_Performance_emission_Sensitivity$SD_precision[aggr_precision_noL_Performance_emission_Sensitivity$Group.2 == 2 &
                                                               aggr_precision_noL_Performance_emission_Sensitivity$Group.1 == 5] )
max(aggr_precision_noL_Performance_emission_Sensitivity$SD_precision[aggr_precision_noL_Performance_emission_Sensitivity$Group.2 == 4&
                                                               aggr_precision_noL_Performance_emission_Sensitivity$Group.1 == 5])
max(aggr_precision_noL_Performance_emission_Sensitivity$SD_precision[aggr_precision_noL_Performance_emission_Sensitivity$Group.2 == 8&
                                                               aggr_precision_noL_Performance_emission_Sensitivity$Group.1 == 5])

max(aggr_precision_noL_Performance_emission_Sensitivity$SD_precision[aggr_precision_noL_Performance_emission_Sensitivity$Group.2 == 2 &
                                                               aggr_precision_noL_Performance_emission_Sensitivity$Group.1 == 7] )
max(aggr_precision_noL_Performance_emission_Sensitivity$SD_precision[aggr_precision_noL_Performance_emission_Sensitivity$Group.2 == 4&
                                                               aggr_precision_noL_Performance_emission_Sensitivity$Group.1 == 7])
max(aggr_precision_noL_Performance_emission_Sensitivity$SD_precision[aggr_precision_noL_Performance_emission_Sensitivity$Group.2 == 8&
                                                               aggr_precision_noL_Performance_emission_Sensitivity$Group.1 == 7])


## Transition probabilities ####
# group_out_3st_gamma_prob <- extracted_results_Sensitivity$group_out_3st_gamma_prob

true_gamma <-  matrix(c(0.7, 0.2, 0.1, 
                        0.1, 0.8, 0.1, 
                        0.1, 0.1, 0.8), ncol = true_m, byrow = TRUE)

Performance_gamma_Sensitivity <- data.frame(sim_iteration = 1:n_sim,
                                    n = rep(settings_full[, 1], each = n_sim * true_m * true_m), 
                                    n_t = rep(settings_full[, 2], each = n_sim * true_m * true_m),
                                    KL_div = rep(settings_full[, 3], each = n_sim * true_m * true_m), 
                                    mu_sc = rep(settings_full[, 4], each = n_sim * true_m * true_m),
                                    mu_K0_sc = rep(settings_full[, 5], each = n_sim * true_m * true_m),
                                    trans_unif = rep(settings_full[, 6], each = n_sim * true_m * true_m), 
                                    trans_K0_sc = rep(settings_full[, 7], each = n_sim * true_m * true_m), 
                                    var_tau_sig = rep(settings_full[, 8], each = n_sim * true_m * true_m),
                                    from_state_i = rep(1:true_m, each = true_m * n_sim),
                                    to_state_j = rep(1:true_m, each = n_sim),
                                    gamma_ij_hat = NA,
                                    gamma_ij_true = rep(as.vector(t(true_gamma)), each = n_sim), 
                                    gamma_ij_95_lower = NA, 
                                    gamma_ij_95_upper = NA) 

it <-  n_sim * true_m * true_m
for(sc in c(1:84)){
  Performance_gamma_Sensitivity[ (1 + (sc-1) * it) : (sc * it),12] <- as.vector(group_out_3st_gamma_prob[[sc]]$median)
  Performance_gamma_Sensitivity[ (1 + (sc-1) * it) : (sc * it),14] <- as.vector(group_out_3st_gamma_prob[[sc]]$quant_2.5)
  Performance_gamma_Sensitivity[ (1 + (sc-1) * it) : (sc * it),15] <- as.vector(group_out_3st_gamma_prob[[sc]]$quant_97.5)
}

View(Performance_gamma_Sensitivity)
saveRDS(Performance_gamma_Sensitivity, file = "Sensitivity/Data/Performance_gamma_Sensitivity.RDS")

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



View(bias_noL_summary)
mean(bias_noL_summary$mean_bias[bias_noL_summary$Diag == 1])
mean(bias_noL_summary$mean_bias[bias_noL_summary$Diag == 0])

mean(bias_noL_summary$mean_bias[bias_noL_summary$Diag == 1 & bias_noL_summary$KL_div == 5])
mean(bias_noL_summary$mean_bias[bias_noL_summary$Diag == 0 & bias_noL_summary$KL_div == 5])

mean(bias_noL_summary$mean_bias[bias_noL_summary$Diag == 1 & bias_noL_summary$KL_div == 7])
mean(bias_noL_summary$mean_bias[bias_noL_summary$Diag == 0 & bias_noL_summary$KL_div == 7])




summary <- aggregate(Performance_gamma_3st, by = list(Performance_gamma_3st$Diag, # without label switching correction 
                                                      Performance_gamma_3st$n_dep, 
                                                      Performance_gamma_3st$KL_div, 
                                                      Performance_gamma_3st$n_t, 
                                                      Performance_gamma_3st$n), FUN = mean)
View(summary)


bias_noL_summary$block[bias_noL_summary$var_tau_sig > 0] <- "BL3_var"
bias_noL_summary$block[bias_noL_summary$trans_unif > 0] <- "BL2_transitions"
bias_noL_summary$block[bias_noL_summary$mu_sc > 0] <- "BL1_emissions"

bias_noL_summary$block <- as.factor(bias_noL_summary$block)  
bias_noL_summary$K0_sc <-bias_noL_summary$mu_K0_sc + bias_noL_summary$trans_K0_sc


ggplot(data = bias_noL_summary, mapping = aes(x = n_t, y = mean_bias, group = interaction(as.factor(n), as.factor(mu_sc),  as.factor(Diag),
                                              as.factor(K0_sc), as.factor(var_tau_sig)), shape = as.factor(K0_sc), color = interaction(as.factor(n)), linetype = interaction(as.factor(mu_sc), as.factor(var_tau_sig)))) +
  geom_point() + 
  geom_line() +
  ylab("Mean bias transition probabilities") +
 # scale_color_discrete(name = "Number of\nsubjects") +
 # scale_linetype_manual(name = "Transition type", 
 #                       values = c(1, 3),
  #                      labels = c("Self-transition", "Between state transition")) +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(block)) +
  xlab("number of observations per subject") +
  geom_hline(yintercept=0, linetype="solid", color = "grey") +
  geom_hline(yintercept=c(-0.04, 0.04), linetype="dashed", color = "grey") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

#### inspecting precision in transtion probs ####


precision_noL_summary <- aggregate(bias_noL_Performance_gamma_Sensitivity, by = list(bias_noL_Performance_gamma_Sensitivity$Diag,
                                                                                     bias_noL_Performance_gamma_Sensitivity$var_tau_sig,
                                                                                     bias_noL_Performance_gamma_Sensitivity$trans_K0_sc,
                                                                                     bias_noL_Performance_gamma_Sensitivity$trans_unif,
                                                                                     bias_noL_Performance_gamma_Sensitivity$mu_K0_sc,
                                                                                     bias_noL_Performance_gamma_Sensitivity$mu_sc,
                                                                                     bias_noL_Performance_gamma_Sensitivity$KL_div, 
                                                                                     bias_noL_Performance_gamma_Sensitivity$n_t,
                                                                             bias_noL_Performance_gamma_Sensitivity$n,
                                                                             bias_noL_Performance_gamma_Sensitivity$from_state_i,
                                                                             bias_noL_Performance_gamma_Sensitivity$to_state_j), FUN = var)

precision_noL_summary <- precision_noL_summary[-which(precision_noL_summary$Group.3 == 5 &
                                                        precision_noL_summary$Group.2 == 2 &
                                                        precision_noL_summary$Group.5 == 120 &
                                                        ( precision_noL_summary$Group.4 == 200 |
                                                            precision_noL_summary$Group.4 == 400)),]


precision_noL_summary$precision <- sqrt(precision_noL_summary$gamma_ij_hat)


mean(precision_noL_summary$precision)
mean(precision_noL_summary$precision[precision_noL_summary$Group.1 == 1 & 
                                       precision_noL_summary$Group.4 == 50])
mean(precision_noL_summary$precision[precision_noL_summary$Group.1 == 1 & 
                                       precision_noL_summary$Group.4 == 200])


## Coverage ####

#### EMission means ####

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


ggplot(data = aggr_cov_noL_Performance_emission_Sensitivity, mapping = aes(x = n_t, y = cov_mean, group = interaction(as.factor(n), as.factor(mu_sc), 
                                                                                                                           as.factor(K0_sc), as.factor(var_tau_sig)), shape = as.factor(K0_sc), color = interaction(as.factor(n)), linetype = interaction(as.factor(mu_sc), as.factor(var_tau_sig)))) +
  geom_point() + 
  geom_line() +
  ylab("Coverage emission distribution") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(block)) +
  xlab("number of observations per subject") +
  geom_hline(yintercept=c(.90, .95), linetype="dashed", color = "grey") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))



#### Emission SDs ####

mean(aggr_cov_noL_Performance_emission_Sensitivity$cov_SD)
mean(aggr_cov_noL_Performance_emission_Sensitivity$cov_SD[aggr_cov_noL_Performance_emission_Sensitivity$KL_div == 5])
mean(aggr_cov_noL_Performance_emission_Sensitivity$cov_SD[aggr_cov_noL_Performance_emission_Sensitivity$KL_div == 7])

mean(aggr_cov_noL_Performance_emission_Sensitivity$cov_SD[aggr_cov_noL_Performance_emission_Sensitivity$n == 30])
mean(aggr_cov_noL_Performance_emission_Sensitivity$cov_SD[aggr_cov_noL_Performance_emission_Sensitivity$n == 120])

mean(aggr_cov_noL_Performance_emission_Sensitivity$cov_SD_width[aggr_cov_noL_Performance_emission_Sensitivity$n_t == 50])
mean(aggr_cov_noL_Performance_emission_Sensitivity$cov_SD_width[aggr_cov_noL_Performance_emission_Sensitivity$n_t == 200])


ggplot(data = aggr_cov_noL_Performance_emission_Sensitivity, mapping = aes(x = n_t, y = cov_SD, group = interaction(as.factor(n), as.factor(mu_sc), 
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


#### Transition probs ####
cov_noL_Performance_gamma_Sensitivity <- bias_noL_Performance_gamma_Sensitivity
cov_noL_Performance_gamma_Sensitivity$cov_trans <- (cov_noL_Performance_gamma_Sensitivity$gamma_ij_true > cov_noL_Performance_gamma_Sensitivity$gamma_ij_95_lower &
                                                      cov_noL_Performance_gamma_Sensitivity$gamma_ij_true < cov_noL_Performance_gamma_Sensitivity$gamma_ij_95_upper) * 1 


#### inspecting coverage in transition probs ####
cov_noL_summary <- aggregate(cov_noL_Performance_gamma_Sensitivity, by = list(cov_noL_Performance_gamma_Sensitivity$Diag, 
                                                                              cov_noL_Performance_gamma_Sensitivity$var_tau_sig,
                                                                              cov_noL_Performance_gamma_Sensitivity$trans_K0_sc,
                                                                              cov_noL_Performance_gamma_Sensitivity$trans_unif,
                                                                              cov_noL_Performance_gamma_Sensitivity$mu_K0_sc,
                                                                              cov_noL_Performance_gamma_Sensitivity$mu_sc,
                                                                              cov_noL_Performance_gamma_Sensitivity$KL_div, 
                                                                              cov_noL_Performance_gamma_Sensitivity$n_t,
                                                                              cov_noL_Performance_gamma_Sensitivity$n), FUN = mean)



mean(cov_noL_summary$cov_trans)
mean(cov_noL_summary$cov_trans[cov_noL_summary$Diag == TRUE])
mean(cov_noL_summary$cov_trans[cov_noL_summary$Diag == FALSE])

cov_noL_summary$block[cov_noL_summary$var_tau_sig > 0] <- "BL3_var"
cov_noL_summary$block[cov_noL_summary$trans_unif > 0] <- "BL2_transitions"
cov_noL_summary$block[cov_noL_summary$mu_sc > 0] <- "BL1_emissions"

cov_noL_summary$block <- as.factor(cov_noL_summary$block)  
cov_noL_summary$K0_sc <-cov_noL_summary$mu_K0_sc + cov_noL_summary$trans_K0_sc



ggplot(data = cov_noL_summary, mapping = aes(x = n_t, y = cov_trans, group = interaction(as.factor(n), as.factor(mu_sc),  as.factor(Diag),
                                                                                          as.factor(K0_sc), as.factor(var_tau_sig)), shape = as.factor(K0_sc), color = interaction(as.factor(n)), linetype = interaction(as.factor(mu_sc), as.factor(var_tau_sig)))) +
  geom_point() + 
  geom_line() +
  ylab("Coverage transition probabilities") +
  # scale_color_discrete(name = "Number of\nsubjects") +
  # scale_linetype_manual(name = "Transition type", 
  #                       values = c(1, 3),
  #                      labels = c("Self-transition", "Between state transition")) +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(block)) +
  xlab("number of observations per subject") +
  geom_hline(yintercept=0, linetype="solid", color = "grey") +
  geom_hline(yintercept=c(.90, .95), linetype="dashed", color = "grey") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

## State decoding  ####

true_states <- extracted_results_Sensitivity$true_states
inferred_states <- extracted_results_Sensitivity$inferred_states
true_state_probs <- extracted_results_Sensitivity$true_state_probs

State_decoding_Sensitivity <- data.frame(sim_iteration = 1:n_sim,
                                 n = rep(settings_full[, 1], each = n_sim), 
                                 n_t = rep(settings_full[, 2], each = n_sim),
                                 KL_div = rep(settings_full[, 3], each = n_sim), 
                                 mu_sc = rep(settings_full[, 4], each = n_sim),
                                 mu_K0_sc = rep(settings_full[, 5], each = n_sim),
                                 trans_unif = rep(settings_full[, 6], each = n_sim), 
                                 trans_K0_sc = rep(settings_full[, 7], each = n_sim), 
                                 var_tau_sig = rep(settings_full[, 8], each = n_sim),
                                 mean_prop_correct = NA,
                                 mean_kappa = NA, 
                                 mean_prop_over_20 = NA
)

for(sc in c(1:84)){
  State_decoding_Sensitivity[(1 + (sc-1) * 100) : (sc * 100),10] <- apply(inferred_states[[sc]][,-1] == true_states[[sc]][,-1], 2, sum) / (dim(true_states[[sc]][,-1])[1])
  inferred_totals <- t(apply(inferred_states[[sc]][,-1],2, table) / 100)
  true_totals <- t(apply(true_states[[sc]][,-1],2, table) / 100)
  matrix_true_inferred <- matrix(NA, nrow = n_sim, ncol = true_m * true_m)
  for(i in 1:n_sim){
    matrix_true_inferred[i,] <- as.vector(table(inferred_states[[sc]][,i + 1], true_states[[sc]][,i + 1])) / 100
  }
  agreements <- matrix(apply(matrix_true_inferred[, c(1, 5, 9)], 1, sum), ncol = 1)
  chance_agreements <- matrix(apply(inferred_totals * true_totals / sum(inferred_totals[1,]), 1, sum), ncol = 1)
  State_decoding_Sensitivity[(1 + (sc-1) * 100) : (sc * 100),11] <- (agreements - chance_agreements) / ( matrix(apply(true_totals, 1, sum), ncol = 1) - chance_agreements)
  State_decoding_Sensitivity[(1 + (sc-1) * 100) : (sc * 100),12] <- apply(true_state_probs[[sc]][,-1] > 0.20, 2, sum) / (dim(true_states[[sc]][,-1])[1])
}

View(State_decoding_Sensitivity)

saveRDS(State_decoding_Sensitivity, file = "Sensitivity/Data/State_decoding_Sensitivity.RDS")


decoding_noL_Sensitivity <- State_decoding_Sensitivity[order(
  State_decoding_Sensitivity$var_tau_sig,
  State_decoding_Sensitivity$trans_K0_sc,
  State_decoding_Sensitivity$trans_unif,
  State_decoding_Sensitivity$mu_K0_sc,
  State_decoding_Sensitivity$mu_sc,
  State_decoding_Sensitivity$KL_div, 
  State_decoding_Sensitivity$n_t, 
  State_decoding_Sensitivity$n,
  State_decoding_Sensitivity$sim_iteration
),]


decoding_noL_Sensitivity$RMSE <- order_aggr_Label_switch_proxy_Sensitivity$RMSE

decoding_noL_Sensitivity <- decoding_noL_Sensitivity[decoding_noL_Sensitivity$RMSE > 0.20,]

#### inspecting state decoding ####
aggr_decoding_noL_Sensitivity <- aggregate(decoding_noL_Sensitivity, by = list( decoding_noL_Sensitivity$var_tau_sig,
                                                                                decoding_noL_Sensitivity$trans_K0_sc,
                                                                                decoding_noL_Sensitivity$trans_unif,
                                                                                decoding_noL_Sensitivity$mu_K0_sc,
                                                                                decoding_noL_Sensitivity$mu_sc,
                                                                                decoding_noL_Sensitivity$KL_div, 
                                                               decoding_noL_Sensitivity$n_t, 
                                                               decoding_noL_Sensitivity$n), FUN = mean)




mean(aggr_decoding_noL_Sensitivity$mean_prop_correct)
min(aggr_decoding_noL_Sensitivity$mean_prop_correct)
max(aggr_decoding_noL_Sensitivity$mean_prop_correct)

mean(aggr_decoding_noL_Sensitivity$mean_kappa)
min(aggr_decoding_noL_Sensitivity$mean_kappa)
max(aggr_decoding_noL_Sensitivity$mean_kappa)

sum(aggr_decoding_noL_Sensitivity$mean_kappa >= 0.80) / dim(aggr_decoding_noL_Sensitivity)[1]

aggr_decoding_noL_Sensitivity$block[aggr_decoding_noL_Sensitivity$var_tau_sig > 0] <- "BL3_var"
aggr_decoding_noL_Sensitivity$block[aggr_decoding_noL_Sensitivity$trans_unif > 0] <- "BL2_transitions"
aggr_decoding_noL_Sensitivity$block[aggr_decoding_noL_Sensitivity$mu_sc > 0] <- "BL1_emissions"

aggr_decoding_noL_Sensitivity$block <- as.factor(aggr_decoding_noL_Sensitivity$block)  
aggr_decoding_noL_Sensitivity$K0_sc <-aggr_decoding_noL_Sensitivity$mu_K0_sc + aggr_decoding_noL_Sensitivity$trans_K0_sc


ggplot(data = aggr_decoding_noL_Sensitivity, mapping = aes(x = n_t, y = mean_prop_correct, group = interaction(as.factor(n), as.factor(mu_sc), 
                                                                                                                      as.factor(K0_sc), as.factor(var_tau_sig)), shape = as.factor(K0_sc), color = interaction(as.factor(n)), linetype = interaction(as.factor(mu_sc), as.factor(var_tau_sig)))) +
  geom_point() + 
  geom_line() +
  ylab("Accuracy") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(block)) +
  xlab("number of observations per subject") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

ggplot(data = aggr_decoding_noL_Sensitivity, mapping = aes(x = n_t, y = mean_kappa, group = interaction(as.factor(n), as.factor(mu_sc), 
                                                                                                               as.factor(K0_sc), as.factor(var_tau_sig)), shape = as.factor(K0_sc), color = interaction(as.factor(n)), linetype = interaction(as.factor(mu_sc), as.factor(var_tau_sig)))) +
  geom_point() + 
  geom_line() +
  ylab("Kappa") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(block)) +
  xlab("number of observations per subject") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))



