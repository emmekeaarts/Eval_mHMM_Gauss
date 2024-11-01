library(mHMMbayes)
library(ggplot2)
library(reshape2)

# Simulstion settings
settings4 <- read.csv("sim_scenarios_V4.csv")
true_m  <- 2

# scenario <- 1
# KL_div  <- settings[scenario, 4]
# n_t     <- settings[scenario, 2]
# n       <- settings[scenario, 1]
# n_dep   <- settings[scenario, 3]
# out_file <- paste0("out_n", n, "_nt", n_t, "_ndep", n_dep, "_", true_m, "st_KLD", KL_div)

## Reading in output files ####
# 2 states
# runs that could be done in 1 go
selected <- sort(c(1:90))
files <- paste0("out_n", settings4[selected, 1], "_nt", settings4[selected, 2], "_ndep", settings4[selected, 3], "_", true_m, "st_KLD", settings4[selected, 4])
out <- vector("list", length(files))
names(out) <- files

# runs splitted over 2 cores 
selected_p2 <- sort(c(121:144))
files_p2_1 <-  paste0("out_n", settings4[selected_p2, 1], "_nt", settings4[selected_p2, 2], "_ndep", settings4[selected_p2, 3], "_", true_m, "st_KLD", settings4[selected_p2, 4], "_1")
files_p2_2 <- paste0("out_n", settings4[selected_p2, 1], "_nt", settings4[selected_p2, 2], "_ndep", settings4[selected_p2, 3], "_", true_m, "st_KLD", settings4[selected_p2, 4], "_2")
out_p2_1 <- vector("list", length(files_p2_1))
out_p2_2 <- vector("list", length(files_p2_2))
names(out_p2_1) <- files_p2_1
names(out_p2_2) <- files_p2_2


selected_p4 <- sort(c(91:120, 145:171))
files_p4_1 <-  paste0("out_n", settings4[selected_p4, 1], "_nt", settings4[selected_p4, 2], "_ndep", settings4[selected_p4, 3], "_", true_m, "st_KLD", settings4[selected_p4, 4], "_1")
files_p4_2 <- paste0("out_n", settings4[selected_p4, 1], "_nt", settings4[selected_p4, 2], "_ndep", settings4[selected_p4, 3], "_", true_m, "st_KLD", settings4[selected_p4, 4], "_2")
files_p4_3 <- paste0("out_n", settings4[selected_p4, 1], "_nt", settings4[selected_p4, 2], "_ndep", settings4[selected_p4, 3], "_", true_m, "st_KLD", settings4[selected_p4, 4], "_3")
files_p4_4 <- paste0("out_n", settings4[selected_p4, 1], "_nt", settings4[selected_p4, 2], "_ndep", settings4[selected_p4, 3], "_", true_m, "st_KLD", settings4[selected_p4, 4], "_4")
out_p4_1 <- vector("list", length(files_p4_1))
out_p4_2 <- vector("list", length(files_p4_2))
out_p4_3 <- vector("list", length(files_p4_3))
out_p4_4 <- vector("list", length(files_p4_4))
names(out_p4_1) <- files_p4_1
names(out_p4_2) <- files_p4_2
names(out_p4_3) <- files_p4_3
names(out_p4_4) <- files_p4_4

group_out_2st_emiss_mean <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
group_out_2st_emiss_sd <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
group_out_2st_gamma_prob <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
group_out_2st_gamma_int <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
means_true <-vector("list", length = length(c(files, files_p2_1, files_p4_1)))
selection_out <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
inferred_states <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
true_state_probs <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
true_states <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))

names(group_out_2st_emiss_mean) <- c(files, files_p2_1, files_p4_1)
names(group_out_2st_emiss_sd) <-c(files, files_p2_1, files_p4_1)
names(group_out_2st_gamma_prob) <- c(files, files_p2_1, files_p4_1)
names(group_out_2st_gamma_int) <- c(files, files_p2_1, files_p4_1)
names(means_true) <- c(files, files_p2_1, files_p4_1)
names(selection_out) <- c(files, files_p2_1, files_p4_1)
names(inferred_states) <- c(files, files_p2_1, files_p4_1)
names(true_state_probs) <- c(files, files_p2_1, files_p4_1)
names(true_states) <- c(files, files_p2_1, files_p4_1)

for(sc in 1:length(files)){
  out[[sc]] <- readRDS(paste0("out_V4/", files[sc], ".rds")) 
  selection_out[[sc]]             <- out[[sc]]$selection_out
  group_out_2st_emiss_mean[[sc]]  <- out[[sc]]$group_out_2st$emiss_mean
  group_out_2st_emiss_sd[[sc]]    <- out[[sc]]$group_out_2st$emiss_sd
  group_out_2st_gamma_prob[[sc]]  <- out[[sc]]$group_out_2st$trans_prob
  group_out_2st_gamma_int[[sc]]   <- out[[sc]]$group_out_2st$trans_interc
  means_true[[sc]]                <- out[[sc]]$means_true
  inferred_states[[sc]]           <- out[[sc]]$inferred_states
  true_state_probs[[sc]]          <- out[[sc]]$true_state_probs
  true_states[[sc]]               <- out[[sc]]$true_states
  if(sc > 1){
    out[[sc]] <- NULL
  }
}

n_sim <- out[[1]]$input$n_sim
J <- out[[1]]$input$J
burn_in <- out[[1]]$input$burn_in

# rm(out)

l_files <- length(files)

for(sc in 1:length(files_p2_1)){
  out_p2_1[[sc]] <- readRDS(paste0("out_V4/", files_p2_1[sc], ".rds")) 
  out_p2_2[[sc]] <- readRDS(paste0("out_V4/", files_p2_2[sc], ".rds")) 
  selection_out[[sc + l_files]]                       <- Map(rbind, out_p2_1[[sc]]$selection_out, out_p2_2[[sc]]$selection_out)
  group_out_2st_emiss_mean[[sc + l_files]]$median     <-  Map(rbind,out_p2_1[[sc]]$group_out_2st$emiss_mean$median, out_p2_2[[sc]]$group_out_2st$emiss_mean$median)
  group_out_2st_emiss_mean[[sc + l_files]]$quant_2.5  <-  Map(rbind,out_p2_1[[sc]]$group_out_2st$emiss_mean$quant_2.5, out_p2_2[[sc]]$group_out_2st$emiss_mean$quant_2.5)
  group_out_2st_emiss_mean[[sc + l_files]]$quant_5    <-  Map(rbind,out_p2_1[[sc]]$group_out_2st$emiss_mean$quant_5, out_p2_2[[sc]]$group_out_2st$emiss_mean$quant_5)
  group_out_2st_emiss_mean[[sc + l_files]]$quant_95   <-  Map(rbind,out_p2_1[[sc]]$group_out_2st$emiss_mean$quant_95, out_p2_2[[sc]]$group_out_2st$emiss_mean$quant_95)
  group_out_2st_emiss_mean[[sc + l_files]]$quant_97.5 <-  Map(rbind,out_p2_1[[sc]]$group_out_2st$emiss_mean$quant_97.5, out_p2_2[[sc]]$group_out_2st$emiss_mean$quant_97.5)
  group_out_2st_emiss_sd[[sc + l_files]]$median       <-  Map(rbind,out_p2_1[[sc]]$group_out_2st$emiss_sd$median, out_p2_2[[sc]]$group_out_2st$emiss_sd$median)
  group_out_2st_emiss_sd[[sc + l_files]]$quant_2.5    <-  Map(rbind,out_p2_1[[sc]]$group_out_2st$emiss_sd$quant_2.5, out_p2_2[[sc]]$group_out_2st$emiss_sd$quant_2.5)
  group_out_2st_emiss_sd[[sc + l_files]]$quant_5      <-  Map(rbind,out_p2_1[[sc]]$group_out_2st$emiss_sd$quant_5, out_p2_2[[sc]]$group_out_2st$emiss_sd$quant_5)
  group_out_2st_emiss_sd[[sc + l_files]]$quant_95     <-  Map(rbind,out_p2_1[[sc]]$group_out_2st$emiss_sd$quant_95, out_p2_2[[sc]]$group_out_2st$emiss_sd$quant_95)
  group_out_2st_emiss_sd[[sc + l_files]]$quant_97.5   <-  Map(rbind,out_p2_1[[sc]]$group_out_2st$emiss_sd$quant_97.5, out_p2_2[[sc]]$group_out_2st$emiss_sd$quant_97.5)
  group_out_2st_gamma_prob[[sc + l_files]]            <-  Map(rbind, out_p2_1[[sc]]$group_out_2st$trans_prob, out_p2_2[[sc]]$group_out_2st$trans_prob)
  group_out_2st_gamma_int[[sc + l_files]]             <-  Map(rbind, out_p2_1[[sc]]$group_out_2st$trans_interc, out_p2_2[[sc]]$group_out_2st$trans_interc)
  means_true[[sc + l_files]]                          <- Map(rbind, out_p2_1[[sc]]$means_true, out_p2_2[[sc]]$means_true)
  inferred_states[[sc + l_files]]                     <- cbind(out_p2_1[[sc]]$inferred_states, out_p2_2[[sc]]$inferred_states[,-1])
  true_state_probs[[sc + l_files]]                    <- cbind(out_p2_1[[sc]]$true_state_probs, out_p2_2[[sc]]$true_state_probs[,-1])
  true_states[[sc + l_files]]                         <- cbind(out_p2_1[[sc]]$true_states, out_p2_2[[sc]]$true_states[,-1])
  out_p2_1[[sc]] <- NULL
  out_p2_2[[sc]] <- NULL
}


# rm(out_p2_1)
# rm(out_p2_2)

l_files_p2 <- length(c(files, files_p2_1))

for(sc in 1:length(files_p4_1)){
  out_p4_1[[sc]] <- readRDS(paste0("out_V4/", files_p4_1[sc], ".rds")) 
  out_p4_2[[sc]] <- readRDS(paste0("out_V4/", files_p4_2[sc], ".rds")) 
  out_p4_3[[sc]] <- readRDS(paste0("out_V4/", files_p4_3[sc], ".rds")) 
  out_p4_4[[sc]] <- readRDS(paste0("out_V4/", files_p4_4[sc], ".rds")) 
  selection_out[[sc + l_files_p2]]          <- Map(rbind, out_p4_1[[sc]]$selection_out, 
                                                   out_p4_2[[sc]]$selection_out, 
                                                   out_p4_3[[sc]]$selection_out, 
                                                   out_p4_4[[sc]]$selection_out)
  group_out_2st_emiss_mean[[sc + l_files_p2]]$median     <-  Map(rbind,out_p4_1[[sc]]$group_out_2st$emiss_mean$median, 
                                                                 out_p4_2[[sc]]$group_out_2st$emiss_mean$median,
                                                                 out_p4_3[[sc]]$group_out_2st$emiss_mean$median,
                                                                 out_p4_4[[sc]]$group_out_2st$emiss_mean$median)
  group_out_2st_emiss_mean[[sc + l_files_p2]]$quant_2.5  <-  Map(rbind,out_p4_1[[sc]]$group_out_2st$emiss_mean$quant_2.5, 
                                                                 out_p4_2[[sc]]$group_out_2st$emiss_mean$quant_2.5,
                                                                 out_p4_3[[sc]]$group_out_2st$emiss_mean$quant_2.5,
                                                                 out_p4_4[[sc]]$group_out_2st$emiss_mean$quant_2.5)
  group_out_2st_emiss_mean[[sc + l_files_p2]]$quant_5    <-  Map(rbind,out_p4_1[[sc]]$group_out_2st$emiss_mean$quant_5, 
                                                                 out_p4_2[[sc]]$group_out_2st$emiss_mean$quant_5,
                                                                 out_p4_3[[sc]]$group_out_2st$emiss_mean$quant_5, 
                                                                 out_p4_4[[sc]]$group_out_2st$emiss_mean$quant_5)
  group_out_2st_emiss_mean[[sc + l_files_p2]]$quant_95   <-  Map(rbind,out_p4_1[[sc]]$group_out_2st$emiss_mean$quant_95, 
                                                                 out_p4_2[[sc]]$group_out_2st$emiss_mean$quant_95,
                                                                 out_p4_3[[sc]]$group_out_2st$emiss_mean$quant_95, 
                                                                 out_p4_4[[sc]]$group_out_2st$emiss_mean$quant_95)
  group_out_2st_emiss_mean[[sc + l_files_p2]]$quant_97.5 <-  Map(rbind,out_p4_1[[sc]]$group_out_2st$emiss_mean$quant_97.5, 
                                                                 out_p4_2[[sc]]$group_out_2st$emiss_mean$quant_97.5,
                                                                 out_p4_3[[sc]]$group_out_2st$emiss_mean$quant_97.5, 
                                                                 out_p4_4[[sc]]$group_out_2st$emiss_mean$quant_97.5)
  group_out_2st_emiss_sd[[sc + l_files_p2]]$median       <-  Map(rbind,out_p4_1[[sc]]$group_out_2st$emiss_sd$median, 
                                                                 out_p4_2[[sc]]$group_out_2st$emiss_sd$median,
                                                                 out_p4_3[[sc]]$group_out_2st$emiss_sd$median, 
                                                                 out_p4_4[[sc]]$group_out_2st$emiss_sd$median)
  group_out_2st_emiss_sd[[sc + l_files_p2]]$quant_2.5    <-  Map(rbind,out_p4_1[[sc]]$group_out_2st$emiss_sd$quant_2.5, 
                                                                 out_p4_2[[sc]]$group_out_2st$emiss_sd$quant_2.5,
                                                                 out_p4_3[[sc]]$group_out_2st$emiss_sd$quant_2.5, 
                                                                 out_p4_4[[sc]]$group_out_2st$emiss_sd$quant_2.5)
  group_out_2st_emiss_sd[[sc + l_files_p2]]$quant_5      <-  Map(rbind,out_p4_1[[sc]]$group_out_2st$emiss_sd$quant_5, 
                                                                 out_p4_2[[sc]]$group_out_2st$emiss_sd$quant_5,
                                                                 out_p4_3[[sc]]$group_out_2st$emiss_sd$quant_5, 
                                                                 out_p4_4[[sc]]$group_out_2st$emiss_sd$quant_5)
  group_out_2st_emiss_sd[[sc + l_files_p2]]$quant_95     <-  Map(rbind,out_p4_1[[sc]]$group_out_2st$emiss_sd$quant_95, 
                                                                 out_p4_2[[sc]]$group_out_2st$emiss_sd$quant_95,
                                                                 out_p4_3[[sc]]$group_out_2st$emiss_sd$quant_95, 
                                                                 out_p4_4[[sc]]$group_out_2st$emiss_sd$quant_95)
  group_out_2st_emiss_sd[[sc + l_files_p2]]$quant_97.5   <-  Map(rbind,out_p4_1[[sc]]$group_out_2st$emiss_sd$quant_97.5, 
                                                                 out_p4_2[[sc]]$group_out_2st$emiss_sd$quant_97.5,
                                                                 out_p4_3[[sc]]$group_out_2st$emiss_sd$quant_97.5, 
                                                                 out_p4_4[[sc]]$group_out_2st$emiss_sd$quant_97.5)
  group_out_2st_gamma_prob[[sc + l_files_p2]]            <-  Map(rbind, out_p4_1[[sc]]$group_out_2st$trans_prob, 
                                                                 out_p4_2[[sc]]$group_out_2st$trans_prob,
                                                                 out_p4_3[[sc]]$group_out_2st$trans_prob, 
                                                                 out_p4_4[[sc]]$group_out_2st$trans_prob)
  group_out_2st_gamma_int[[sc + l_files_p2]]             <-  Map(rbind, out_p4_1[[sc]]$group_out_2st$trans_interc, 
                                                                 out_p4_2[[sc]]$group_out_2st$trans_interc,
                                                                 out_p4_3[[sc]]$group_out_2st$trans_interc, 
                                                                 out_p4_4[[sc]]$group_out_2st$trans_interc)
  means_true[[sc + l_files_p2]]                          <- Map(rbind, out_p4_1[[sc]]$means_true, 
                                                                out_p4_2[[sc]]$means_true,
                                                                out_p4_3[[sc]]$means_true, 
                                                                out_p4_4[[sc]]$means_true)
  inferred_states[[sc + l_files_p2]]                     <- cbind(out_p4_1[[sc]]$inferred_states, 
                                                                  out_p4_2[[sc]]$inferred_states[,-1],
                                                                  out_p4_3[[sc]]$inferred_states[,-1],
                                                                  out_p4_4[[sc]]$inferred_states[,-1])
  true_state_probs[[sc + l_files_p2]]                    <- cbind(out_p4_1[[sc]]$true_state_probs, 
                                                                  out_p4_2[[sc]]$true_state_probs[,-1],
                                                                  out_p4_3[[sc]]$true_state_probs[,-1],
                                                                  out_p4_4[[sc]]$true_state_probs[,-1])
  true_states[[sc + l_files_p2]]                         <- cbind(out_p4_1[[sc]]$true_states, 
                                                                  out_p4_2[[sc]]$true_states[,-1],
                                                                  out_p4_3[[sc]]$true_states[,-1],
                                                                  out_p4_4[[sc]]$true_states[,-1])
  out_p4_1[[sc]] <- NULL
  out_p4_2[[sc]] <- NULL
  out_p4_3[[sc]] <- NULL
  out_p4_4[[sc]] <- NULL
}

# rm(out_p4_1)
# rm(out_p4_2)
# rm(out_p4_3)
# rm(out_p4_4)


## SAVE EXTRACTED FILES ####

extracted_results_2st <- list(group_out_2st_emiss_mean = group_out_2st_emiss_mean,
                              group_out_2st_emiss_sd = group_out_2st_emiss_sd,
                              group_out_2st_gamma_prob = group_out_2st_gamma_prob,
                              group_out_2st_gamma_int = group_out_2st_gamma_int,
                              means_true = means_true,
                              selection_out = selection_out,
                              inferred_states = inferred_states,
                              true_state_probs = true_state_probs,
                              true_states = true_states, 
                              n_sim = n_sim, 
                              J = J, 
                              burn_in = burn_in)
saveRDS(extracted_results_2st, file = "extracted_results_2st.RDS")

# extracted_results_2st <- readRDS("extracted_results_2st.RDS")

# Model selection ####

Model_selection_2st <- data.frame(sim_iteration = 1:n_sim,
                                  true_m = true_m,
                                  n = rep(settings4[c(selected, selected_p2, selected_p4), 1], each = n_sim * 4), 
                                  n_t = rep(settings4[c(selected, selected_p2, selected_p4), 2], each = n_sim * 4),
                                  n_dep = rep(settings4[c(selected, selected_p2, selected_p4), 3], each = n_sim * 4),
                                  KL_div = rep(settings4[c(selected, selected_p2, selected_p4), 4], each = n_sim * 4), 
                                  K_candidate = rep(1:4, each = n_sim),
                                  AIC = NA,
                                  AICc = NA
)

for(sc in 1:171){
  for(kc in 1:4){
    Model_selection_2st[(1 + (sc - 1) * n_sim * 4 + n_sim * (kc - 1)) : 
                          (100 + (sc - 1) * n_sim * 4 + n_sim * (kc - 1)) ,8] <- selection_out[[sc]]$AIC_mean[,kc]
    Model_selection_2st[(1 + (sc - 1) * n_sim * 4 + n_sim * (kc - 1)) : 
                          (100 + (sc - 1) * n_sim * 4 + n_sim * (kc - 1)) ,9] <- selection_out[[sc]]$AICc_mean[,kc]
  }
}

View(Model_selection_2st)
saveRDS(Model_selection_2st, file = "result_tables/Model_selection_2st.RDS")



# Model performance 

## Emission distribution ####
Performance_emission_2st <- data.frame(sim_iteration = 1:n_sim,
                                       true_m = true_m,
                                       n = rep(settings4[c(selected, selected_p2, selected_p4), 1], times = n_sim * true_m * settings4[c(selected, selected_p2, selected_p4), 3]), 
                                       n_t = rep(settings4[c(selected, selected_p2, selected_p4), 2], times = n_sim * true_m * settings4[c(selected, selected_p2, selected_p4), 3]),
                                       n_dep = rep(settings4[c(selected, selected_p2, selected_p4), 3], times = n_sim * true_m * settings4[c(selected, selected_p2, selected_p4), 3]),
                                       KL_div = rep(settings4[c(selected, selected_p2, selected_p4), 4], times = n_sim * true_m * settings4[c(selected, selected_p2, selected_p4), 3]), 
                                       k = rep(1:true_m, each = n_sim),
                                       dep = rep(sequence(settings4[c(selected, selected_p2, selected_p4), 3], 1, 1), each = true_m * n_sim),
                                       mean_hat = NA,
                                       mean_true = NA, 
                                       mean_95_lower = NA, 
                                       mean_95_upper = NA, 
                                       SD_hat = NA, 
                                       SD_true = NA,
                                       SD_95_lower = NA,  
                                       SD_95_upper = NA
)


# group_out_2st_emiss_mean <- extracted_results_2st$group_out_2st_emiss_mean
# group_out_2st_emiss_sd <- extracted_results_2st$group_out_2st_emiss_sd
# means_true <- extracted_results_2st$means_true

it <- n_sim * true_m * settings4[c(selected, selected_p2, selected_p4), 3]
cumsum_it <- c(0,cumsum(it))

for(sc in 1:171){
  Performance_emission_2st[(1 + cumsum_it[sc]) : (it[sc] + cumsum_it[sc]),9] <- unlist(group_out_2st_emiss_mean[[sc]]$median)
  Performance_emission_2st[(1 + cumsum_it[sc]) : (it[sc] + cumsum_it[sc]),10] <- unlist(means_true[[sc]])  
  Performance_emission_2st[(1 + cumsum_it[sc]) : (it[sc] + cumsum_it[sc]),11] <- unlist(group_out_2st_emiss_mean[[sc]]$quant_2.5)
  Performance_emission_2st[(1 + cumsum_it[sc]) : (it[sc] + cumsum_it[sc]),12] <- unlist(group_out_2st_emiss_mean[[sc]]$quant_97.5)
  Performance_emission_2st[(1 + cumsum_it[sc]) : (it[sc] + cumsum_it[sc]),13] <- unlist(group_out_2st_emiss_sd[[sc]]$median)
  Performance_emission_2st[(1 + cumsum_it[sc]) : (it[sc] + cumsum_it[sc]),14] <- 1
  Performance_emission_2st[(1 + cumsum_it[sc]) : (it[sc] + cumsum_it[sc]),15] <- unlist(group_out_2st_emiss_sd[[sc]]$quant_2.5)
  Performance_emission_2st[(1 + cumsum_it[sc]) : (it[sc] + cumsum_it[sc]),16] <- unlist(group_out_2st_emiss_sd[[sc]]$quant_97.5)
}
View(Performance_emission_2st)

saveRDS(Performance_emission_2st, file = "result_tables/Performance_emission_2st.RDS")

### Explaining bias in emission ####
#### Quantifying likely label switching problem occurences ####
rmse <- function(x){
  row_mean <- mean(x)
  rmse <- sqrt(sum((mean(x) - x)^2 / true_m))
  return(rmse)
}

Label_switch_proxy_2st <- data.frame(sim_iteration = 1:n_sim,
                                     true_m = true_m,
                                     n = rep(settings4[c(selected, selected_p2, selected_p4), 1], times = n_sim  * settings4[c(selected, selected_p2, selected_p4), 3]), 
                                     n_t = rep(settings4[c(selected, selected_p2, selected_p4), 2], times = n_sim  * settings4[c(selected, selected_p2, selected_p4), 3]),
                                     n_dep = rep(settings4[c(selected, selected_p2, selected_p4), 3], times = n_sim  * settings4[c(selected, selected_p2, selected_p4), 3]),
                                     KL_div = rep(settings4[c(selected, selected_p2, selected_p4), 4], times = n_sim  * settings4[c(selected, selected_p2, selected_p4), 3]), 
                                     dep = rep(sequence(settings4[c(selected, selected_p2, selected_p4), 3], 1, 1), each = n_sim),
                                     RMSE = NA)


it <- n_sim  * settings4[c(selected, selected_p2, selected_p4), 3]
cumsum_it <-  c(0,cumsum(n_sim  * settings4[c(selected, selected_p2, selected_p4), 3]))

for(sc in 1:length(c(selected, selected_p2, selected_p4))){
  Label_switch_proxy_2st[(1 + cumsum_it[sc]) : (it[sc] + cumsum_it[sc]),8] <- unlist(lapply(group_out_2st_emiss_mean[[sc]]$median, apply, 1, rmse))
}

aggregate(Label_switch_proxy_2st, by = list(Label_switch_proxy_2st$n, 
                                            Label_switch_proxy_2st$n_t, 
                                            Label_switch_proxy_2st$n_dep, 
                                            Label_switch_proxy_2st$KL_div), FUN = mean)


saveRDS(Label_switch_proxy_2st, file = "result_tables/Label_switch_proxy_2st.RDS")



## Transition probabilities ####
# group_out_2st_gamma_prob <- extracted_results_2st$group_out_2st_gamma_prob

true_gamma <- matrix(c(0.7, 0.3, 
                       0.2, 0.8), ncol = true_m, byrow = TRUE)

Performance_gamma_2st <- data.frame(sim_iteration = 1:n_sim,
                                    true_m = true_m,
                                    n = rep(settings4[c(selected, selected_p2, selected_p4), 1], each = n_sim * true_m * true_m), 
                                    n_t = rep(settings4[c(selected, selected_p2, selected_p4), 2], each = n_sim * true_m * true_m),
                                    n_dep = rep(settings4[c(selected, selected_p2, selected_p4), 3], each = n_sim * true_m * true_m),
                                    KL_div = rep(settings4[c(selected, selected_p2, selected_p4), 4], each = n_sim * true_m * true_m), 
                                    from_state_i = rep(1:true_m, each = true_m * n_sim),
                                    to_state_j = rep(1:true_m, each = n_sim),
                                    gamma_ij_hat = NA,
                                    gamma_ij_true = rep(as.vector(t(true_gamma)), each = n_sim), 
                                    gamma_ij_95_lower = NA, 
                                    gamma_ij_95_upper = NA) 

it <-  n_sim * true_m * true_m
for(sc in 1:171){
  Performance_gamma_2st[ (1 + (sc-1) * it) : (sc * it),9] <- as.vector(group_out_2st_gamma_prob[[sc]]$median)
  Performance_gamma_2st[ (1 + (sc-1) * it) : (sc * it),11] <- as.vector(group_out_2st_gamma_prob[[sc]]$quant_2.5)
  Performance_gamma_2st[ (1 + (sc-1) * it) : (sc * it),12] <- as.vector(group_out_2st_gamma_prob[[sc]]$quant_97.5)
}

View(Performance_gamma_2st)

saveRDS(Performance_gamma_2st, file = "result_tables/Performance_gamma_2st.RDS")


## State decoding  ####

true_states <- extracted_results_2st$true_states
inferred_states <- extracted_results_2st$inferred_states
true_state_probs <- extracted_results_2st$true_state_probs

State_decoding_2st <- data.frame(sim_iteration = 1:n_sim,
                                 true_m = true_m,
                                 n = rep(settings4[c(selected, selected_p2, selected_p4), 1], each = n_sim), 
                                 n_t = rep(settings4[c(selected, selected_p2, selected_p4), 2], each = n_sim),
                                 n_dep = rep(settings4[c(selected, selected_p2, selected_p4), 3], each = n_sim),
                                 KL_div = rep(settings4[c(selected, selected_p2, selected_p4), 4], each = n_sim), 
                                 mean_prop_correct = NA,
                                 mean_kappa = NA, 
                                 mean_prop_over_20 = NA
)

for(sc in 1:171){
  State_decoding_2st[(1 + (sc-1) * 100) : (sc * 100),7] <- apply(inferred_states[[sc]][,-1] == true_states[[sc]][,-1], 2, sum) / (dim(true_states[[sc]][,-1])[1])
  inferred_totals <- t(apply(inferred_states[[sc]][,-1],2, table) / 100)
  true_totals <- t(apply(true_states[[sc]][,-1],2, table) / 100)
  matrix_true_inferred <- matrix(NA, nrow = n_sim, ncol = true_m * true_m)
  for(i in 1:n_sim){
    matrix_true_inferred[i,] <- as.vector(table(inferred_states[[sc]][,i + 1], true_states[[sc]][,i + 1])) / 100
  }
  agreements <- matrix(apply(matrix_true_inferred[, c(1, 4)], 1, sum), ncol = 1)
  chance_agreements <- matrix(apply(inferred_totals * true_totals / sum(inferred_totals[1,]), 1, sum), ncol = 1)
  State_decoding_2st[(1 + (sc-1) * 100) : (sc * 100),8] <- (agreements - chance_agreements) / ( matrix(apply(true_totals, 1, sum), ncol = 1) - chance_agreements)
  State_decoding_2st[(1 + (sc-1) * 100) : (sc * 100),9] <- apply(true_state_probs[[sc]][,-1] > 0.20, 2, sum) / (dim(true_states[[sc]][,-1])[1])
}

View(State_decoding_2st)

saveRDS(State_decoding_2st, file = "result_tables/State_decoding_2st.RDS")




