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
# 3 states
# runs that could be done in 1 go
selected <- sort(c(1:90))
files <- paste0("out_n", settings4[selected, 1], "_nt", settings4[selected, 2], "_ndep", settings4[selected, 3], "_", true_m, "st_KLD", settings4[selected, 4])
out0 <- out1 <- out2 <- vector("list", length(files))
names(out0) <- files
names(out1) <- files
names(out2) <- files

# runs splitted over 2 cores 
selected_p2 <- sort(c(121:144))
files_p2_1 <-  paste0("out_n", settings4[selected_p2, 1], "_nt", settings4[selected_p2, 2], "_ndep", settings4[selected_p2, 3], "_", true_m, "st_KLD", settings4[selected_p2, 4])
files_p2_2 <- paste0("out_n", settings4[selected_p2, 1], "_nt", settings4[selected_p2, 2], "_ndep", settings4[selected_p2, 3], "_", true_m, "st_KLD", settings4[selected_p2, 4])
out0_p2_1 <- vector("list", length(files_p2_1))
out1_p2_1 <- vector("list", length(files_p2_1))
out2_p2_1 <- vector("list", length(files_p2_1))
out0_p2_2 <- vector("list", length(files_p2_2))
out1_p2_2 <- vector("list", length(files_p2_2))
out2_p2_2 <- vector("list", length(files_p2_2))
names(out0_p2_1) <- files_p2_1
names(out1_p2_1) <- files_p2_1
names(out2_p2_1) <- files_p2_1
names(out0_p2_2) <- files_p2_2
names(out1_p2_2) <- files_p2_2
names(out2_p2_2) <- files_p2_2


selected_p4 <- sort(c(91:120, 145:171))
files_p4_1 <-  paste0("out_n", settings4[selected_p4, 1], "_nt", settings4[selected_p4, 2], "_ndep", settings4[selected_p4, 3], "_", true_m, "st_KLD", settings4[selected_p4, 4])
files_p4_2 <- paste0("out_n", settings4[selected_p4, 1], "_nt", settings4[selected_p4, 2], "_ndep", settings4[selected_p4, 3], "_", true_m, "st_KLD", settings4[selected_p4, 4])
files_p4_3 <- paste0("out_n", settings4[selected_p4, 1], "_nt", settings4[selected_p4, 2], "_ndep", settings4[selected_p4, 3], "_", true_m, "st_KLD", settings4[selected_p4, 4])
files_p4_4 <- paste0("out_n", settings4[selected_p4, 1], "_nt", settings4[selected_p4, 2], "_ndep", settings4[selected_p4, 3], "_", true_m, "st_KLD", settings4[selected_p4, 4])
out0_p4_1 <- vector("list", length(files_p4_1))
out0_p4_2 <- vector("list", length(files_p4_2))
out0_p4_3 <- vector("list", length(files_p4_3))
out0_p4_4 <- vector("list", length(files_p4_4))
out1_p4_1 <- vector("list", length(files_p4_1))
out1_p4_2 <- vector("list", length(files_p4_2))
out1_p4_3 <- vector("list", length(files_p4_3))
out1_p4_4 <- vector("list", length(files_p4_4))
out2_p4_1 <- vector("list", length(files_p4_1))
out2_p4_2 <- vector("list", length(files_p4_2))
out2_p4_3 <- vector("list", length(files_p4_3))
out2_p4_4 <- vector("list", length(files_p4_4))
names(out0_p4_1) <- files_p4_1
names(out0_p4_2) <- files_p4_2
names(out0_p4_3) <- files_p4_3
names(out0_p4_4) <- files_p4_4
names(out1_p4_1) <- files_p4_1
names(out1_p4_2) <- files_p4_2
names(out1_p4_3) <- files_p4_3
names(out1_p4_4) <- files_p4_4
names(out2_p4_1) <- files_p4_1
names(out2_p4_2) <- files_p4_2
names(out2_p4_3) <- files_p4_3
names(out2_p4_4) <- files_p4_4

group_out_2st_emiss_mean0 <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
group_out_2st_emiss_mean1 <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
group_out_2st_emiss_mean2 <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
means_true0 <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
means_true1 <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
means_true2 <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
gamma_full_conv0 <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
gamma_full_conv1 <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
gamma_full_conv2 <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
emiss_mean_full_conv0 <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
emiss_mean_full_conv1 <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
emiss_mean_full_conv2 <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))

names(means_true0) <- c(files)
names(means_true1) <- c(files)
names(means_true2) <- c(files)
names(gamma_full_conv0) <- c(files)
names(gamma_full_conv1) <- c(files)
names(gamma_full_conv2) <- c(files)
names(emiss_mean_full_conv0) <- c(files)
names(emiss_mean_full_conv1) <- c(files)
names(emiss_mean_full_conv2) <- c(files)


for(sc in 1:length(files)){
  out0[[sc]] <- readRDS(paste0("out_V4/", files[sc], ".rds")) 
  out1[[sc]] <- readRDS(paste0("out_V4_conv/", files[sc], "_conv1.rds")) 
  out2[[sc]] <- readRDS(paste0("out_V4_conv/", files[sc], "_conv2.rds")) 
  means_true0[[sc]]                <- out0[[sc]]$means_true
  means_true1[[sc]]                <- out1[[sc]]$means_true_conv1
  means_true2[[sc]]                <- out2[[sc]]$means_true_conv2
  group_out_2st_emiss_mean0[[sc]]$median  <- out0[[sc]]$group_out_2st$emiss_mean$median
  group_out_2st_emiss_mean1[[sc]]$median  <- out1[[sc]]$group_out_2st_conv1$emiss_mean$median
  group_out_2st_emiss_mean2[[sc]]$median  <- out2[[sc]]$group_out_2st_conv2$emiss_mean$median
  gamma_full_conv0[[sc]]           <- list(out0[[sc]]$gamma_full[[1]][[true_m]],
                                           out0[[sc]]$gamma_full[[2]][[true_m]],
                                           out0[[sc]]$gamma_full[[3]][[true_m]],
                                           out0[[sc]]$gamma_full[[4]][[true_m]])
  gamma_full_conv1[[sc]]           <- out1[[sc]]$gamma_full_conv1
  gamma_full_conv2[[sc]]           <- out2[[sc]]$gamma_full_conv2
  emiss_mean_full_conv0[[sc]]     <- list(out0[[sc]]$emiss_mean_full[[1]][[true_m]],
                                          out0[[sc]]$emiss_mean_full[[2]][[true_m]],
                                          out0[[sc]]$emiss_mean_full[[3]][[true_m]],
                                          out0[[sc]]$emiss_mean_full[[4]][[true_m]])
  emiss_mean_full_conv1[[sc]]     <- out1[[sc]]$emiss_mean_full_conv1
  emiss_mean_full_conv2[[sc]]     <- out2[[sc]]$emiss_mean_full_conv2
  if(sc > 1){
    out0[[sc]] <- NULL
    out1[[sc]] <- NULL
    out2[[sc]] <- NULL
  }
}


n_sim <- out0[[1]]$input$n_sim
J <- out0[[1]]$input$J
burn_in <- out0[[1]]$input$burn_in

# rm(out0)
# rm(out1)
# rm(out2)

l_files <- length(files)

for(sc in 1:length(files_p2_1)){
  out0_p2_1[[sc]] <- readRDS(paste0("out_V4/", files_p2_1[sc], "_1.rds")) 
  out0_p2_2[[sc]] <- readRDS(paste0("out_V4/", files_p2_2[sc], "_2.rds")) 
  out1_p2_1[[sc]] <- readRDS(paste0("out_V4_conv/", files_p2_1[sc], "_conv1_1.rds")) 
  out1_p2_2[[sc]] <- readRDS(paste0("out_V4_conv/", files_p2_2[sc], "_conv1_2.rds")) 
  out2_p2_1[[sc]] <- readRDS(paste0("out_V4_conv/", files_p2_1[sc], "_conv2_1.rds")) 
  out2_p2_2[[sc]] <- readRDS(paste0("out_V4_conv/", files_p2_2[sc], "_conv2_2.rds")) 
  means_true0[[sc]]                <- Map(rbind, out0_p2_1[[sc]]$means_true, out0_p2_2[[sc]]$means_true)
  means_true1[[sc]]                <- Map(rbind, out1_p2_1[[sc]]$means_true, out1_p2_2[[sc]]$means_true)
  means_true2[[sc]]                <- Map(rbind, out2_p2_1[[sc]]$means_true, out2_p2_2[[sc]]$means_true)
  group_out_2st_emiss_mean0[[sc + l_files]]$median     <-  Map(rbind,out0_p2_1[[sc]]$group_out_2st$emiss_mean$median, 
                                                               out0_p2_2[[sc]]$group_out_2st$emiss_mean$median)
  group_out_2st_emiss_mean1[[sc + l_files]]$median  <- Map(rbind, out1_p2_1[[sc]]$group_out_2st_conv1$emiss_mean$median,
                                                           out1_p2_2[[sc]]$group_out_2st_conv1$emiss_mean$median)
  group_out_2st_emiss_mean2[[sc + l_files]]$median  <- Map(rbind, out2_p2_1[[sc]]$group_out_2st_conv2$emiss_mean$median,
                                                           out2_p2_2[[sc]]$group_out_2st_conv2$emiss_mean$median)
  gamma_full_conv0[[sc + l_files]]           <- list(out0_p2_1[[sc]]$gamma_full[[1]][[true_m]], out0_p2_1[[sc]]$gamma_full[[2]][[true_m]],
                                                     out0_p2_2[[sc]]$gamma_full[[1]][[true_m]], out0_p2_2[[sc]]$gamma_full[[2]][[true_m]])
  gamma_full_conv1[[sc + l_files]]           <- list(out1_p2_1[[sc]]$gamma_full[[1]], out1_p2_1[[sc]]$gamma_full[[2]],
                                                     out1_p2_2[[sc]]$gamma_full[[1]], out1_p2_2[[sc]]$gamma_full[[2]])
  gamma_full_conv2[[sc + l_files]]           <- list(out2_p2_1[[sc]]$gamma_full[[1]], out2_p2_1[[sc]]$gamma_full[[2]],
                                                     out2_p2_2[[sc]]$gamma_full[[1]], out2_p2_2[[sc]]$gamma_full[[2]])
  emiss_mean_full_conv0[[sc + l_files]]     <- list(out0_p2_1[[sc]]$emiss_mean_full[[1]][[true_m]], out0_p2_1[[sc]]$emiss_mean_full[[2]][[true_m]], 
                                                    out0_p2_2[[sc]]$emiss_mean_full[[1]][[true_m]], out0_p2_2[[sc]]$emiss_mean_full[[2]][[true_m]])
  emiss_mean_full_conv1[[sc + l_files]]     <- list(out1_p2_1[[sc]]$emiss_mean_full_conv1[[1]], out1_p2_1[[sc]]$emiss_mean_full_conv1[[2]],
                                                    out1_p2_2[[sc]]$emiss_mean_full_conv1[[1]], out1_p2_2[[sc]]$emiss_mean_full_conv1[[2]])
  emiss_mean_full_conv2[[sc + l_files]]     <- list(out2_p2_1[[sc]]$emiss_mean_full_conv2[[1]], out2_p2_1[[sc]]$emiss_mean_full_conv2[[2]],
                                                    out2_p2_2[[sc]]$emiss_mean_full_conv2[[1]], out2_p2_2[[sc]]$emiss_mean_full_conv2[[2]])
  out0_p2_1[[sc]] <- NULL
  out0_p2_1[[sc]] <- NULL
  out1_p2_1[[sc]] <- NULL
  out1_p2_2[[sc]] <- NULL
  out2_p2_1[[sc]] <- NULL
  out2_p2_2[[sc]] <- NULL
}

# rm(out0_p2_1)
# rm(out0_p2_2)
# rm(out1_p2_1)
# rm(out1_p2_2)
# rm(out2_p2_1)
# rm(out2_p2_2)

l_files_p2 <- length(c(files, files_p2_1))

for(sc in 1:length(files_p4_1)){
  out0_p4_1[[sc]] <- readRDS(paste0("out_V4/", files_p4_1[sc], "_1.rds")) 
  out0_p4_2[[sc]] <- readRDS(paste0("out_V4/", files_p4_2[sc], "_2.rds")) 
  out0_p4_3[[sc]] <- readRDS(paste0("out_V4/", files_p4_3[sc], "_3.rds")) 
  out0_p4_4[[sc]] <- readRDS(paste0("out_V4/", files_p4_4[sc], "_4.rds")) 
  
  out1_p4_1[[sc]] <- readRDS(paste0("out_V4_conv/", files_p4_1[sc], "_conv1_1.rds")) 
  out1_p4_2[[sc]] <- readRDS(paste0("out_V4_conv/", files_p4_2[sc], "_conv1_2.rds")) 
  out1_p4_3[[sc]] <- readRDS(paste0("out_V4_conv/", files_p4_3[sc], "_conv1_3.rds")) 
  out1_p4_4[[sc]] <- readRDS(paste0("out_V4_conv/", files_p4_4[sc], "_conv1_4.rds")) 
  
  out2_p4_1[[sc]] <- readRDS(paste0("out_V4_conv/", files_p4_1[sc], "_conv2_1.rds")) 
  out2_p4_2[[sc]] <- readRDS(paste0("out_V4_conv/", files_p4_2[sc], "_conv2_2.rds")) 
  out2_p4_3[[sc]] <- readRDS(paste0("out_V4_conv/", files_p4_3[sc], "_conv2_3.rds")) 
  out2_p4_4[[sc]] <- readRDS(paste0("out_V4_conv/", files_p4_4[sc], "_conv2_4.rds")) 
  
  means_true0[[sc  + l_files_p2]]                <- Map(rbind, out0_p4_1[[sc]]$means_true, out0_p4_2[[sc]]$means_true, 
                                                        out0_p4_3[[sc]]$means_true, out0_p4_4[[sc]]$means_true)
  means_true1[[sc  + l_files_p2]]                <- Map(rbind, out1_p4_1[[sc]]$means_true, out1_p4_2[[sc]]$means_true,
                                                        out1_p4_3[[sc]]$means_true, out1_p4_4[[sc]]$means_true)
  means_true2[[sc  + l_files_p2]]                <- Map(rbind, out2_p4_1[[sc]]$means_true, out2_p4_2[[sc]]$means_true,
                                                        out2_p4_3[[sc]]$means_true, out2_p4_4[[sc]]$means_true)
  group_out_2st_emiss_mean0[[sc + l_files_p2]]$median     <-  Map(rbind,out0_p4_1[[sc]]$group_out_2st$emiss_mean$median, 
                                                                  out0_p4_2[[sc]]$group_out_2st$emiss_mean$median,
                                                                  out0_p4_3[[sc]]$group_out_2st$emiss_mean$median,
                                                                  out0_p4_4[[sc]]$group_out_2st$emiss_mean$median)
  group_out_2st_emiss_mean1[[sc + l_files_p2]]$median  <- Map(rbind, out1_p4_1[[sc]]$group_out_2st_conv1$emiss_mean$median,
                                                              out1_p4_2[[sc]]$group_out_2st_conv1$emiss_mean$median,
                                                              out1_p4_3[[sc]]$group_out_2st_conv1$emiss_mean$median,
                                                              out1_p4_4[[sc]]$group_out_2st_conv1$emiss_mean$median)
  group_out_2st_emiss_mean2[[sc + l_files_p2]]$median  <- Map(rbind, out2_p4_1[[sc]]$group_out_2st_conv2$emiss_mean$median,
                                                              out2_p4_2[[sc]]$group_out_2st_conv2$emiss_mean$median,
                                                              out2_p4_3[[sc]]$group_out_2st_conv2$emiss_mean$median,
                                                              out2_p4_4[[sc]]$group_out_2st_conv2$emiss_mean$median)
  gamma_full_conv0[[sc  + l_files_p2]]           <- list(out0_p4_1[[sc]]$gamma_full[[1]][[true_m]], out0_p4_2[[sc]]$gamma_full[[1]][[true_m]],
                                                         out0_p4_3[[sc]]$gamma_full[[1]][[true_m]], out0_p4_4[[sc]]$gamma_full[[1]][[true_m]])
  gamma_full_conv1[[sc  + l_files_p2]]           <- list(out1_p4_1[[sc]]$gamma_full[[1]], out1_p4_2[[sc]]$gamma_full[[1]],
                                                         out1_p4_3[[sc]]$gamma_full[[1]], out1_p4_4[[sc]]$gamma_full[[1]])
  gamma_full_conv2[[sc  + l_files_p2]]           <- list(out2_p4_1[[sc]]$gamma_full[[1]], out2_p4_2[[sc]]$gamma_full[[1]],
                                                         out2_p4_3[[sc]]$gamma_full[[1]], out2_p4_4[[sc]]$gamma_full[[1]])
  emiss_mean_full_conv0[[sc  + l_files_p2]]     <- list(out0_p4_1[[sc]]$emiss_mean_full[[1]][[true_m]], out0_p4_2[[sc]]$emiss_mean_full[[1]][[true_m]], 
                                                        out0_p4_3[[sc]]$emiss_mean_full[[1]][[true_m]], out0_p4_4[[sc]]$emiss_mean_full[[1]][[true_m]])
  emiss_mean_full_conv1[[sc  + l_files_p2]]     <- list(out1_p4_1[[sc]]$emiss_mean_full_conv1[[1]], out1_p4_2[[sc]]$emiss_mean_full_conv1[[1]],
                                                        out1_p4_3[[sc]]$emiss_mean_full_conv1[[1]], out1_p4_4[[sc]]$emiss_mean_full_conv1[[1]])
  emiss_mean_full_conv2[[sc  + l_files_p2]]     <- list(out2_p4_1[[sc]]$emiss_mean_full_conv2[[1]], out2_p4_2[[sc]]$emiss_mean_full_conv2[[1]],
                                                        out2_p4_3[[sc]]$emiss_mean_full_conv2[[1]], out2_p4_4[[sc]]$emiss_mean_full_conv2[[1]])
  
  out0_p4_1[[sc]] <- NULL
  out0_p4_2[[sc]] <- NULL
  out0_p4_3[[sc]] <- NULL
  out0_p4_4[[sc]] <- NULL
  
  out1_p4_1[[sc]] <- NULL
  out1_p4_2[[sc]] <- NULL
  out1_p4_3[[sc]] <- NULL
  out1_p4_4[[sc]] <- NULL
  
  out2_p4_1[[sc]] <- NULL
  out2_p4_2[[sc]] <- NULL
  out2_p4_3[[sc]] <- NULL
  out2_p4_4[[sc]] <- NULL
}

# rm(out0_p4_1)
# rm(out0_p4_2)
# rm(out0_p4_3)
# rm(out0_p4_4)
# rm(out1_p4_1)
# rm(out1_p4_2)
# rm(out1_p4_3)
# rm(out1_p4_4)
# rm(out2_p4_1)
# rm(out2_p4_2)
# rm(out2_p4_3)
# rm(out2_p4_4)

## SAVE EXTRACTED FILES ####

extracted_convergence_2st <- list(means_true0 = means_true0,
                                  means_true1 = means_true1,
                                  means_true2 = means_true2,
                                  group_out_2st_emiss_mean0 = group_out_2st_emiss_mean0,
                                  group_out_2st_emiss_mean1 = group_out_2st_emiss_mean1,
                                  group_out_2st_emiss_mean2 = group_out_2st_emiss_mean2,
                                  gamma_full_conv0 = gamma_full_conv0,
                                  gamma_full_conv1 = gamma_full_conv1,
                                  gamma_full_conv2 = gamma_full_conv2,
                                  emiss_mean_full_conv0 = emiss_mean_full_conv0,
                                  emiss_mean_full_conv1 = emiss_mean_full_conv1,
                                  emiss_mean_full_conv2 = emiss_mean_full_conv2, 
                                  n_sim = n_sim, 
                                  J = J, 
                                  burn_in = burn_in)
saveRDS(extracted_convergence_2st, file = "extracted_convergence_2st.RDS")

## First ASSESS label switching problem ####
#### Quantifying likely label switching problem occurences ####
rmse <- function(x){
  row_mean <- mean(x)
  rmse <- sqrt(sum((mean(x) - x)^2 / true_m))
  return(rmse)
}

# Do per convergence chain and partition 
# RUN 1, FULL RUN
Label_switch_proxy_2st0p1 <- data.frame(sim_iteration = 1:n_sim,
                                        true_m = true_m,
                                        n = rep(settings4[c(selected), 1], times = n_sim  * settings4[c(selected), 3]), 
                                        n_t = rep(settings4[c(selected), 2], times = n_sim  * settings4[c(selected), 3]),
                                        n_dep = rep(settings4[c(selected), 3], times = n_sim  * settings4[c(selected), 3]),
                                        KL_div = rep(settings4[c(selected), 4], times = n_sim  * settings4[c(selected), 3]), 
                                        dep = rep(sequence(settings4[c(selected), 3], 1, 1), each = n_sim),
                                        RMSE = NA)


it <- n_sim  * settings4[c(selected), 3]
cumsum_it <-  c(0,cumsum(n_sim  * settings4[c(selected), 3]))

for(sc in selected){
  Label_switch_proxy_2st0p1[(1 + cumsum_it[sc]) : (it[sc] + cumsum_it[sc]),8] <- unlist(lapply(group_out_2st_emiss_mean0[[sc]]$median, apply, 1, rmse))
}

View(Label_switch_proxy_2st0p1)



Label_switch_proxy_2st0p2 <- data.frame(sim_iteration = 1:n_sim,
                                        true_m = true_m,
                                        n = rep(settings4[c(selected_p2), 1], times = n_sim  * settings4[c(selected_p2), 3]), 
                                        n_t = rep(settings4[c(selected_p2), 2], times = n_sim  * settings4[c(selected_p2), 3]),
                                        n_dep = rep(settings4[c(selected_p2), 3], times = n_sim  * settings4[c(selected_p2), 3]),
                                        KL_div = rep(settings4[c(selected_p2), 4], times = n_sim  * settings4[c(selected_p2), 3]), 
                                        dep = rep(sequence(settings4[c(selected_p2), 3], 1, 1), each = n_sim),
                                        RMSE = NA)


it <- n_sim  * settings4[c(selected_p2), 3]
cumsum_it <-  c(0,cumsum(n_sim  * settings4[c(selected_p2), 3]))

for(sc in 1:length(selected_p2)){
  Label_switch_proxy_2st0p2[(1 + cumsum_it[sc]) : (it[sc] + cumsum_it[sc]),8] <- unlist(lapply(group_out_2st_emiss_mean0[[sc + l_files]]$median, apply, 1, rmse))
}

View(Label_switch_proxy_2st0p2)



Label_switch_proxy_2st0p4 <- data.frame(sim_iteration = 1:n_sim,
                                        true_m = true_m,
                                        n = rep(settings4[c(selected_p4), 1], times = n_sim  * settings4[c(selected_p4), 3]), 
                                        n_t = rep(settings4[c(selected_p4), 2], times = n_sim  * settings4[c(selected_p4), 3]),
                                        n_dep = rep(settings4[c(selected_p4), 3], times = n_sim  * settings4[c(selected_p4), 3]),
                                        KL_div = rep(settings4[c(selected_p4), 4], times = n_sim  * settings4[c(selected_p4), 3]), 
                                        dep = rep(sequence(settings4[c(selected_p4), 3], 1, 1), each = n_sim),
                                        RMSE = NA)


it <- n_sim  * settings4[c(selected_p4), 3]
cumsum_it <-  c(0,cumsum(n_sim  * settings4[c(selected_p4), 3]))

for(sc in 1:length(c(selected_p4))){
  Label_switch_proxy_2st0p4[(1 + cumsum_it[sc]) : (it[sc] + cumsum_it[sc]),8] <- unlist(lapply(group_out_2st_emiss_mean0[[sc + l_files_p2]]$median, apply, 1, rmse))
}

View(Label_switch_proxy_2st0p4)



# RUN 2, only 4 runs (for p1 and p2, by accident 6 runs for p4)
Label_switch_proxy_2st1p1 <- data.frame(sim_iteration = 1:4,
                                        true_m = true_m,
                                        n = rep(settings4[c(selected), 1], times = 4  * settings4[c(selected), 3]), 
                                        n_t = rep(settings4[c(selected), 2], times = 4  * settings4[c(selected), 3]),
                                        n_dep = rep(settings4[c(selected), 3], times = 4  * settings4[c(selected), 3]),
                                        KL_div = rep(settings4[c(selected), 4], times = 4  * settings4[c(selected), 3]), 
                                        dep = rep(sequence(settings4[c(selected), 3], 1, 1), each = 4),
                                        RMSE = NA)


it <- 4  * settings4[c(selected), 3]
cumsum_it <-  c(0,cumsum(4  * settings4[c(selected), 3]))

for(sc in selected){
  Label_switch_proxy_2st1p1[(1 + cumsum_it[sc]) : (it[sc] + cumsum_it[sc]),8] <- unlist(lapply(group_out_2st_emiss_mean1[[sc]]$median, apply, 1, rmse))
}

View(Label_switch_proxy_2st1p1)



Label_switch_proxy_2st1p2 <- data.frame(sim_iteration = 1:4,
                                        true_m = true_m,
                                        n = rep(settings4[c(selected_p2), 1], times = 4  * settings4[c(selected_p2), 3]), 
                                        n_t = rep(settings4[c(selected_p2), 2], times = 4  * settings4[c(selected_p2), 3]),
                                        n_dep = rep(settings4[c(selected_p2), 3], times = 4  * settings4[c(selected_p2), 3]),
                                        KL_div = rep(settings4[c(selected_p2), 4], times = 4  * settings4[c(selected_p2), 3]), 
                                        dep = rep(sequence(settings4[c(selected_p2), 3], 1, 1), each = 4),
                                        RMSE = NA)


it <- 4  * settings4[c(selected_p2), 3]
cumsum_it <-  c(0,cumsum(4  * settings4[c(selected_p2), 3]))

for(sc in 1:length(selected_p2)){
  Label_switch_proxy_2st1p2[(1 + cumsum_it[sc]) : (it[sc] + cumsum_it[sc]),8] <- unlist(lapply(group_out_2st_emiss_mean1[[sc + l_files]]$median, apply, 1, rmse))
}

View(Label_switch_proxy_2st1p2)



Label_switch_proxy_2st1p4 <- data.frame(sim_iteration = 1:4,
                                        true_m = true_m,
                                        n = rep(settings4[c(selected_p4), 1], times = 4  * settings4[c(selected_p4), 3]), 
                                        n_t = rep(settings4[c(selected_p4), 2], times = 4  * settings4[c(selected_p4), 3]),
                                        n_dep = rep(settings4[c(selected_p4), 3], times = 4  * settings4[c(selected_p4), 3]),
                                        KL_div = rep(settings4[c(selected_p4), 4], times = 4  * settings4[c(selected_p4), 3]), 
                                        dep = rep(sequence(settings4[c(selected_p4), 3], 1, 1), each = 4),
                                        RMSE = NA)


it <- 4  * settings4[c(selected_p4), 3]
cumsum_it <-  c(0,cumsum(4  * settings4[c(selected_p4), 3]))

for(sc in 1:length(c(selected_p4))){
  Label_switch_proxy_2st1p4[(1 + cumsum_it[sc]) : (it[sc] + cumsum_it[sc]),8] <- unlist(lapply(group_out_2st_emiss_mean1[[sc + l_files_p2]]$median, apply, 1, rmse))
}

View(Label_switch_proxy_2st1p4)

# RUN 3, only 4 runs (for p1 and p2, by accident 6 runs for p4)
Label_switch_proxy_2st2p1 <- data.frame(sim_iteration = 1:4,
                                        true_m = true_m,
                                        n = rep(settings4[c(selected), 1], times = 4  * settings4[c(selected), 3]), 
                                        n_t = rep(settings4[c(selected), 2], times = 4  * settings4[c(selected), 3]),
                                        n_dep = rep(settings4[c(selected), 3], times = 4  * settings4[c(selected), 3]),
                                        KL_div = rep(settings4[c(selected), 4], times = 4  * settings4[c(selected), 3]), 
                                        dep = rep(sequence(settings4[c(selected), 3], 1, 1), each = 4),
                                        RMSE = NA)


it <- 4  * settings4[c(selected), 3]
cumsum_it <-  c(0,cumsum(4  * settings4[c(selected), 3]))

for(sc in selected){
  Label_switch_proxy_2st2p1[(1 + cumsum_it[sc]) : (it[sc] + cumsum_it[sc]),8] <- unlist(lapply(group_out_2st_emiss_mean2[[sc]]$median, apply, 1, rmse))
}

View(Label_switch_proxy_2st2p1)



Label_switch_proxy_2st2p2 <- data.frame(sim_iteration = 1:4,
                                        true_m = true_m,
                                        n = rep(settings4[c(selected_p2), 1], times = 4  * settings4[c(selected_p2), 3]), 
                                        n_t = rep(settings4[c(selected_p2), 2], times = 4  * settings4[c(selected_p2), 3]),
                                        n_dep = rep(settings4[c(selected_p2), 3], times = 4  * settings4[c(selected_p2), 3]),
                                        KL_div = rep(settings4[c(selected_p2), 4], times = 4  * settings4[c(selected_p2), 3]), 
                                        dep = rep(sequence(settings4[c(selected_p2), 3], 1, 1), each = 4),
                                        RMSE = NA)


it <- 4  * settings4[c(selected_p2), 3]
cumsum_it <-  c(0,cumsum(4  * settings4[c(selected_p2), 3]))

for(sc in 1:length(selected_p2)){
  Label_switch_proxy_2st2p2[(1 + cumsum_it[sc]) : (it[sc] + cumsum_it[sc]),8] <- unlist(lapply(group_out_2st_emiss_mean2[[sc + l_files]]$median, apply, 1, rmse))
}

View(Label_switch_proxy_2st2p2)



Label_switch_proxy_2st2p4 <- data.frame(sim_iteration = 1:4,
                                        true_m = true_m,
                                        n = rep(settings4[c(selected_p4), 1], times = 4  * settings4[c(selected_p4), 3]), 
                                        n_t = rep(settings4[c(selected_p4), 2], times = 4  * settings4[c(selected_p4), 3]),
                                        n_dep = rep(settings4[c(selected_p4), 3], times = 4  * settings4[c(selected_p4), 3]),
                                        KL_div = rep(settings4[c(selected_p4), 4], times = 4  * settings4[c(selected_p4), 3]), 
                                        dep = rep(sequence(settings4[c(selected_p4), 3], 1, 1), each = 4),
                                        RMSE = NA)


it <- 4  * settings4[c(selected_p4), 3]
cumsum_it <-  c(0,cumsum(4  * settings4[c(selected_p4), 3]))

for(sc in 1:length(c(selected_p4))){
  Label_switch_proxy_2st2p4[(1 + cumsum_it[sc]) : (it[sc] + cumsum_it[sc]),8] <- unlist(lapply(group_out_2st_emiss_mean2[[sc + l_files_p2]]$median, apply, 1, rmse))
}

View(Label_switch_proxy_2st2p4)



### aggregating over dependent variables within a simrun, convergence run 1 ####

aggr_Label_switch_proxy_2st0p1 <- aggregate(Label_switch_proxy_2st0p1, by = list(Label_switch_proxy_2st0p1$sim_iteration, 
                                                                                 Label_switch_proxy_2st0p1$KL_div,
                                                                                 Label_switch_proxy_2st0p1$n_dep, 
                                                                                 Label_switch_proxy_2st0p1$n_t, 
                                                                                 Label_switch_proxy_2st0p1$n), FUN = mean)

aggr_Label_switch_proxy_2st0p2 <- aggregate(Label_switch_proxy_2st0p2, by = list(Label_switch_proxy_2st0p2$sim_iteration, 
                                                                                 Label_switch_proxy_2st0p2$KL_div,
                                                                                 Label_switch_proxy_2st0p2$n_dep, 
                                                                                 Label_switch_proxy_2st0p2$n_t, 
                                                                                 Label_switch_proxy_2st0p2$n), FUN = mean)

aggr_Label_switch_proxy_2st0p4 <- aggregate(Label_switch_proxy_2st0p4, by = list(Label_switch_proxy_2st0p4$sim_iteration, 
                                                                                 Label_switch_proxy_2st0p4$KL_div,
                                                                                 Label_switch_proxy_2st0p4$n_dep, 
                                                                                 Label_switch_proxy_2st0p4$n_t, 
                                                                                 Label_switch_proxy_2st0p4$n), FUN = mean)

#### Making correct sub selection 
aggr_Label_switch_proxy_2st0p1_select <- aggr_Label_switch_proxy_2st0p1[aggr_Label_switch_proxy_2st0p1$sim_iteration %in% c(20, 40, 60, 80),]
aggr_Label_switch_proxy_2st0p2_select <- aggr_Label_switch_proxy_2st0p2[aggr_Label_switch_proxy_2st0p2$sim_iteration %in% c(20, 40, (50 + 20), (50 + 40)),]
aggr_Label_switch_proxy_2st0p4_select <- aggr_Label_switch_proxy_2st0p4[aggr_Label_switch_proxy_2st0p4$sim_iteration %in% c(20, (25+20), (50 + 20), (75 + 20)),]


### aggregating over dependent variables within a simrun, convergence run 2 ####

aggr_Label_switch_proxy_2st1p1 <- aggregate(Label_switch_proxy_2st1p1, by = list(Label_switch_proxy_2st1p1$sim_iteration, 
                                                                                 Label_switch_proxy_2st1p1$KL_div,
                                                                                 Label_switch_proxy_2st1p1$n_dep, 
                                                                                 Label_switch_proxy_2st1p1$n_t, 
                                                                                 Label_switch_proxy_2st1p1$n), FUN = mean)

aggr_Label_switch_proxy_2st1p2 <- aggregate(Label_switch_proxy_2st1p2, by = list(Label_switch_proxy_2st1p2$sim_iteration, 
                                                                                 Label_switch_proxy_2st1p2$KL_div,
                                                                                 Label_switch_proxy_2st1p2$n_dep, 
                                                                                 Label_switch_proxy_2st1p2$n_t, 
                                                                                 Label_switch_proxy_2st1p2$n), FUN = mean)

aggr_Label_switch_proxy_2st1p4 <- aggregate(Label_switch_proxy_2st1p4, by = list(Label_switch_proxy_2st1p4$sim_iteration, 
                                                                                 Label_switch_proxy_2st1p4$KL_div,
                                                                                 Label_switch_proxy_2st1p4$n_dep, 
                                                                                 Label_switch_proxy_2st1p4$n_t, 
                                                                                 Label_switch_proxy_2st1p4$n), FUN = mean)

### aggregating over dependent variables within a simrun, convergence run 3 ####
aggr_Label_switch_proxy_2st2p1 <- aggregate(Label_switch_proxy_2st2p1, by = list(Label_switch_proxy_2st2p1$sim_iteration, 
                                                                                 Label_switch_proxy_2st2p1$KL_div,
                                                                                 Label_switch_proxy_2st2p1$n_dep, 
                                                                                 Label_switch_proxy_2st2p1$n_t, 
                                                                                 Label_switch_proxy_2st2p1$n), FUN = mean)

aggr_Label_switch_proxy_2st2p2 <- aggregate(Label_switch_proxy_2st2p2, by = list(Label_switch_proxy_2st2p2$sim_iteration, 
                                                                                 Label_switch_proxy_2st2p2$KL_div,
                                                                                 Label_switch_proxy_2st2p2$n_dep, 
                                                                                 Label_switch_proxy_2st2p2$n_t, 
                                                                                 Label_switch_proxy_2st2p2$n), FUN = mean)

aggr_Label_switch_proxy_2st2p4 <- aggregate(Label_switch_proxy_2st2p4, by = list(Label_switch_proxy_2st2p4$sim_iteration, 
                                                                                 Label_switch_proxy_2st2p4$KL_div,
                                                                                 Label_switch_proxy_2st2p4$n_dep, 
                                                                                 Label_switch_proxy_2st2p4$n_t, 
                                                                                 Label_switch_proxy_2st2p4$n), FUN = mean)

# Creating label switching indicator

aggr_Label_switch_proxy_2st0p1_select$present <- (aggr_Label_switch_proxy_2st0p1_select$RMSE < 0.20 ) * 1
aggr_Label_switch_proxy_2st0p2_select$present <- (aggr_Label_switch_proxy_2st0p2_select$RMSE < 0.20 ) * 1
aggr_Label_switch_proxy_2st0p4_select$present <- (aggr_Label_switch_proxy_2st0p4_select$RMSE < 0.20 ) * 1

aggr_Label_switch_proxy_2st1p1$present1 <- (aggr_Label_switch_proxy_2st1p1$RMSE < 0.20 ) * 1
aggr_Label_switch_proxy_2st1p2$present1 <- (aggr_Label_switch_proxy_2st1p2$RMSE < 0.20 ) * 1
aggr_Label_switch_proxy_2st1p4$present1 <- (aggr_Label_switch_proxy_2st1p4$RMSE < 0.20 ) * 1

aggr_Label_switch_proxy_2st2p1$present2 <- (aggr_Label_switch_proxy_2st2p1$RMSE < 0.20 ) * 1
aggr_Label_switch_proxy_2st2p2$present2 <- (aggr_Label_switch_proxy_2st2p2$RMSE < 0.20 ) * 1
aggr_Label_switch_proxy_2st2p4$present2 <- (aggr_Label_switch_proxy_2st2p4$RMSE < 0.20 ) * 1


aggr_Label_switch_proxy_2st1p1$present0 <- aggr_Label_switch_proxy_2st0p1_select$present
aggr_Label_switch_proxy_2st1p2$present0 <- aggr_Label_switch_proxy_2st0p2_select$present 
aggr_Label_switch_proxy_2st1p4$present0 <- aggr_Label_switch_proxy_2st0p4_select$present 


aggr_Label_switch_proxy_2st1p1$present2 <- aggr_Label_switch_proxy_2st2p1$present2
aggr_Label_switch_proxy_2st1p2$present2 <- aggr_Label_switch_proxy_2st2p2$present2 
aggr_Label_switch_proxy_2st1p4$present2 <- aggr_Label_switch_proxy_2st2p4$present2


aggr_Label_switch_proxy_2st1p1$presentSUM <- aggr_Label_switch_proxy_2st1p1$present0 + 
  aggr_Label_switch_proxy_2st1p1$present1 +
  aggr_Label_switch_proxy_2st1p1$present2

aggr_Label_switch_proxy_2st1p1$no_chain <- (aggr_Label_switch_proxy_2st1p1$presentSUM != 0) * 1
aggr2_Label_switch_proxy_2st1p1 <- aggregate(aggr_Label_switch_proxy_2st1p1, by = list(aggr_Label_switch_proxy_2st1p1$KL_div,
                                                                                       aggr_Label_switch_proxy_2st1p1$n_dep, 
                                                                                       aggr_Label_switch_proxy_2st1p1$n_t, 
                                                                                       aggr_Label_switch_proxy_2st1p1$n), FUN = sum)

table(aggr_Label_switch_proxy_2st1p1$presentSUM)

aggr_Label_switch_proxy_2st1p2$presentSUM <- aggr_Label_switch_proxy_2st1p2$present0 + 
  aggr_Label_switch_proxy_2st1p2$present1 +
  aggr_Label_switch_proxy_2st1p2$present2

aggr_Label_switch_proxy_2st1p2$no_chain <- (aggr_Label_switch_proxy_2st1p2$presentSUM != 0) * 1
aggr2_Label_switch_proxy_2st1p2 <- aggregate(aggr_Label_switch_proxy_2st1p2, by = list(aggr_Label_switch_proxy_2st1p2$KL_div,
                                                                                       aggr_Label_switch_proxy_2st1p2$n_dep, 
                                                                                       aggr_Label_switch_proxy_2st1p2$n_t, 
                                                                                       aggr_Label_switch_proxy_2st1p2$n), FUN = sum)


table(aggr_Label_switch_proxy_2st1p2$presentSUM)

aggr_Label_switch_proxy_2st1p4$presentSUM <- aggr_Label_switch_proxy_2st1p4$present0 + 
  aggr_Label_switch_proxy_2st1p4$present1 +
  aggr_Label_switch_proxy_2st1p4$present2

table(aggr_Label_switch_proxy_2st1p4$presentSUM)

aggr_Label_switch_proxy_2st1p4$no_chain <- (aggr_Label_switch_proxy_2st1p4$presentSUM != 0) * 1
aggr2_Label_switch_proxy_2st1p4 <- aggregate(aggr_Label_switch_proxy_2st1p4, by = list(aggr_Label_switch_proxy_2st1p4$KL_div,
                                                                                       aggr_Label_switch_proxy_2st1p4$n_dep, 
                                                                                       aggr_Label_switch_proxy_2st1p4$n_t, 
                                                                                       aggr_Label_switch_proxy_2st1p4$n), FUN = sum)


# extracted_convergence_2st <- readRDS("extracted_convergence_2st.RDS")
## GR statistics transition probabilities ####
gamma_full_conv0 <- extracted_convergence_2st$gamma_full_conv0
gamma_full_conv1 <- extracted_convergence_2st$gamma_full_conv1
gamma_full_conv2 <- extracted_convergence_2st$gamma_full_conv2

Gelman_rubin_trans <- function(data1, data2, data3, m, J, burn_in, n_chain){
  L <- J - burn_in
  gamma_mean_data1 <- if(!is.null(data1)){apply(data1[((burn_in + 1) : J),], 2, mean)} else{NA}
  gamma_mean_data2 <- if(!is.null(data2)){apply(data2[((burn_in + 1) : J),], 2, mean)} else{NA}
  gamma_mean_data3 <- if(!is.null(data3)){apply(data3[((burn_in + 1) : J),], 2, mean)} else{NA}
  gamma_grand_mean <- apply(cbind(gamma_mean_data1,gamma_mean_data2,gamma_mean_data3),1, sum, na.rm = TRUE) / n_chain
  
  
  
  gamma_between_chain_var <- L / (n_chain-1) * apply(cbind((gamma_mean_data1 - gamma_grand_mean)^2, 
                                                           (gamma_mean_data2 - gamma_grand_mean)^2,
                                                           (gamma_mean_data3 - gamma_grand_mean)^2), 1, sum, na.rm = TRUE)
  gamma_within_chain_var_d1 <- if(!is.null(data1)){(1 / (n_chain-1)) * apply(data1[((burn_in + 1) : J),], 2, var)} else{NA}
  gamma_within_chain_var_d2 <- if(!is.null(data2)){(1 / (n_chain-1)) * apply(data2[((burn_in + 1) : J),], 2, var)} else{NA}
  gamma_within_chain_var_d3 <- if(!is.null(data3)){(1 / (n_chain-1)) * apply(data3[((burn_in + 1) : J),], 2, var)} else{NA}
  gamma_W <- (1/n_chain) * apply(cbind(gamma_within_chain_var_d1, gamma_within_chain_var_d2, gamma_within_chain_var_d3), 1, sum, na.rm = TRUE)
  
  gamma_GelRub <- sqrt(((L - 1)/ L * gamma_W + 1/L *  gamma_between_chain_var ) / gamma_W)
  return(gamma_GelRub)
}

sc <- 1
sim_it <- 1
n_sim <- 4
Convergence_gamma_2st_p1 <- data.frame(sim_rep = rep(1:4, each = true_m * true_m),
                                       true_m = true_m,
                                       n = rep(settings4[c(selected), 1], each = true_m * true_m * n_sim), 
                                       n_t = rep(settings4[c(selected), 2], each = true_m * true_m * n_sim),
                                       n_dep = rep(settings4[c(selected), 3], each = true_m * true_m * n_sim),
                                       KL_div = rep(settings4[c(selected), 4], each = true_m * true_m * n_sim), 
                                       from_state_i = rep(1:true_m, each = true_m),
                                       to_state_j = seq(1:true_m),
                                       GR_stat = NA,
                                       valid_runs = NA) 
Convergence_gamma_2st_p2 <- data.frame(sim_rep = rep(1:4, each = true_m * true_m),
                                       true_m = true_m,
                                       n = rep(settings4[c(selected_p2), 1], each = true_m * true_m * n_sim), 
                                       n_t = rep(settings4[c(selected_p2), 2], each = true_m * true_m * n_sim),
                                       n_dep = rep(settings4[c(selected_p2), 3], each = true_m * true_m * n_sim),
                                       KL_div = rep(settings4[c(selected_p2), 4], each = true_m * true_m * n_sim), 
                                       from_state_i = rep(1:true_m, each = true_m),
                                       to_state_j = seq(1:true_m),
                                       GR_stat = NA,
                                       valid_runs = NA) 
Convergence_gamma_2st_p4 <- data.frame(sim_rep = rep(1:4, each = true_m * true_m),
                                       true_m = true_m,
                                       n = rep(settings4[c(selected_p4), 1], each = true_m * true_m * n_sim), 
                                       n_t = rep(settings4[c(selected_p4), 2], each = true_m * true_m * n_sim),
                                       n_dep = rep(settings4[c(selected_p4), 3], each = true_m * true_m * n_sim),
                                       KL_div = rep(settings4[c(selected_p4), 4], each = true_m * true_m * n_sim), 
                                       from_state_i = rep(1:true_m, each = true_m),
                                       to_state_j = seq(1:true_m),
                                       GR_stat = NA,
                                       valid_runs = NA) 


for(sc in 1:length(selected)){
  for(sim_it in 1:4){
    spacer <- (sc-1) * true_m * true_m * n_sim
    ind <- (1 + (sim_it-1) * true_m * true_m + spacer)
    
    valid <- 4- aggr2_Label_switch_proxy_2st1p1$no_chain[which(aggr2_Label_switch_proxy_2st1p1$Group.4 == Convergence_gamma_2st_p1$n[ind] &
                                                                 aggr2_Label_switch_proxy_2st1p1$Group.3 == Convergence_gamma_2st_p1$n_t[ind] &
                                                                 aggr2_Label_switch_proxy_2st1p1$Group.2 == Convergence_gamma_2st_p1$n_dep[ind] &
                                                                 aggr2_Label_switch_proxy_2st1p1$Group.1 == Convergence_gamma_2st_p1$KL_div[ind])]
    Convergence_gamma_2st_p1[ind: 
                               ((sim_it) * true_m * true_m + spacer), 10] <- valid
    
    sum_present <- aggr_Label_switch_proxy_2st1p1$presentSUM[which(aggr_Label_switch_proxy_2st1p1$sim_iteration == sim_it &
                                                                     aggr_Label_switch_proxy_2st1p1$n == Convergence_gamma_2st_p1$n[ind] &
                                                                     aggr_Label_switch_proxy_2st1p1$n_t == Convergence_gamma_2st_p1$n_t[ind] &
                                                                     aggr_Label_switch_proxy_2st1p1$n_dep == Convergence_gamma_2st_p1$n_dep[ind] &
                                                                     aggr_Label_switch_proxy_2st1p1$KL_div == Convergence_gamma_2st_p1$KL_div[ind])]
    if(sum_present == 0){
      Convergence_gamma_2st_p1[ind:
                                 ((sim_it) * true_m * true_m + spacer) ,9] <- Gelman_rubin_trans(data1 = gamma_full_conv0[[sc]][[sim_it]],  
                                                                                                 data2 = gamma_full_conv1[[sc]][[sim_it]],  
                                                                                                 data3 = gamma_full_conv2[[sc]][[sim_it]], 
                                                                                                 m = true_m, 
                                                                                                 J = J, 
                                                                                                 burn_in = burn_in, 
                                                                                                 n_chain = 3)
    } 
  }
}

View(Convergence_gamma_2st_p1)


for(sc in 1:length(selected_p2)){
  for(sim_it in 1:4){
    spacer <- (sc-1) * true_m * true_m * n_sim
    ind <- (1 + (sim_it-1) * true_m * true_m + spacer)
    
    valid <- 4- aggr2_Label_switch_proxy_2st1p2$no_chain[which(aggr2_Label_switch_proxy_2st1p2$Group.4 == Convergence_gamma_2st_p2$n[ind] &
                                                                 aggr2_Label_switch_proxy_2st1p2$Group.3 == Convergence_gamma_2st_p2$n_t[ind] &
                                                                 aggr2_Label_switch_proxy_2st1p2$Group.2 == Convergence_gamma_2st_p2$n_dep[ind] &
                                                                 aggr2_Label_switch_proxy_2st1p2$Group.1 == Convergence_gamma_2st_p2$KL_div[ind])]
    Convergence_gamma_2st_p2[ind: 
                               ((sim_it) * true_m * true_m + spacer), 10] <- valid
    
    sum_present <- aggr_Label_switch_proxy_2st1p2$presentSUM[which(aggr_Label_switch_proxy_2st1p2$sim_iteration == sim_it &
                                                                     aggr_Label_switch_proxy_2st1p2$n == Convergence_gamma_2st_p2$n[ind] &
                                                                     aggr_Label_switch_proxy_2st1p2$n_t == Convergence_gamma_2st_p2$n_t[ind] &
                                                                     aggr_Label_switch_proxy_2st1p2$n_dep == Convergence_gamma_2st_p2$n_dep[ind] &
                                                                     aggr_Label_switch_proxy_2st1p2$KL_div == Convergence_gamma_2st_p2$KL_div[ind])]
    if(sum_present == 0){
      Convergence_gamma_2st_p2[ind:
                                 ((sim_it) * true_m * true_m + spacer) ,9] <- Gelman_rubin_trans(data1 = gamma_full_conv0[[sc + l_files]][[sim_it]],  
                                                                                                 data2 = gamma_full_conv1[[sc + l_files]][[sim_it]],  
                                                                                                 data3 = gamma_full_conv2[[sc + l_files]][[sim_it]], 
                                                                                                 m = true_m, 
                                                                                                 J = J, 
                                                                                                 burn_in = burn_in, 
                                                                                                 n_chain = 3)
    } 
  }
}

View(Convergence_gamma_2st_p2)


for(sc in 1:length(selected_p4)){
  for(sim_it in 1:4){
    spacer <- (sc-1) * true_m * true_m * n_sim
    ind <- (1 + (sim_it-1) * true_m * true_m + spacer)
    
    valid <- 4- aggr2_Label_switch_proxy_2st1p4$no_chain[which(aggr2_Label_switch_proxy_2st1p4$Group.4 == Convergence_gamma_2st_p4$n[ind] &
                                                                 aggr2_Label_switch_proxy_2st1p4$Group.3 == Convergence_gamma_2st_p4$n_t[ind] &
                                                                 aggr2_Label_switch_proxy_2st1p4$Group.2 == Convergence_gamma_2st_p4$n_dep[ind] &
                                                                 aggr2_Label_switch_proxy_2st1p4$Group.1 == Convergence_gamma_2st_p4$KL_div[ind])]
    Convergence_gamma_2st_p4[ind: 
                               ((sim_it) * true_m * true_m + spacer), 10] <- valid
    
    sum_present <- aggr_Label_switch_proxy_2st1p4$presentSUM[which(aggr_Label_switch_proxy_2st1p4$sim_iteration == sim_it &
                                                                     aggr_Label_switch_proxy_2st1p4$n == Convergence_gamma_2st_p4$n[ind] &
                                                                     aggr_Label_switch_proxy_2st1p4$n_t == Convergence_gamma_2st_p4$n_t[ind] &
                                                                     aggr_Label_switch_proxy_2st1p4$n_dep == Convergence_gamma_2st_p4$n_dep[ind] &
                                                                     aggr_Label_switch_proxy_2st1p4$KL_div == Convergence_gamma_2st_p4$KL_div[ind])]
    if(sum_present == 0){
      Convergence_gamma_2st_p4[ind:
                                 ((sim_it) * true_m * true_m + spacer) ,9] <- Gelman_rubin_trans(data1 = gamma_full_conv0[[sc + l_files]][[sim_it]],  
                                                                                                 data2 = gamma_full_conv1[[sc + l_files]][[sim_it]],  
                                                                                                 data3 = gamma_full_conv2[[sc + l_files]][[sim_it]], 
                                                                                                 m = true_m, 
                                                                                                 J = J, 
                                                                                                 burn_in = burn_in, 
                                                                                                 n_chain = 3)
    } 
  }
}

View(Convergence_gamma_2st_p4)


### merging transition convergence results over partitions ####
Convergence_gamma_2st_corrected <- rbind(Convergence_gamma_2st_p1,
                                         Convergence_gamma_2st_p2,
                                         Convergence_gamma_2st_p4)

saveRDS(Convergence_gamma_2st_corrected, file = "result_tables/Convergence_gamma_2st_corrected_V2.RDS")


Convergence_gamma_2st_corrected$prop <- (Convergence_gamma_2st_corrected$GR_stat > 1.1) * 1
Convergence_gamma_2st_corrected$prop_conv <- (Convergence_gamma_2st_corrected$GR_stat < 1.1) * 1

sum(Convergence_gamma_2st_corrected$prop_conv, na.rm = TRUE) / sum(!is.na(Convergence_gamma_2st_corrected$prop_conv))
sum(Convergence_gamma_2st_corrected$prop_conv[Convergence_gamma_2st_corrected$KL_div==3], na.rm = TRUE) / sum(!is.na(Convergence_gamma_2st_corrected$prop_conv[Convergence_gamma_2st_corrected$KL_div==3]))
sum(Convergence_gamma_2st_corrected$prop_conv[Convergence_gamma_2st_corrected$KL_div==5], na.rm = TRUE) / sum(!is.na(Convergence_gamma_2st_corrected$prop_conv[Convergence_gamma_2st_corrected$KL_div==5]))
sum(Convergence_gamma_2st_corrected$prop_conv[Convergence_gamma_2st_corrected$KL_div==7], na.rm = TRUE) / sum(!is.na(Convergence_gamma_2st_corrected$prop_conv[Convergence_gamma_2st_corrected$KL_div==7]))

Subset_Convergence_gamma_2st_corrected <- Convergence_gamma_2st_corrected[Convergence_gamma_2st_corrected$KL_div > 3 &
                                                                            Convergence_gamma_2st_corrected$valid_runs > 2,]

sum(Subset_Convergence_gamma_2st_corrected$prop_conv, na.rm = TRUE) /sum(!is.na(Subset_Convergence_gamma_2st_corrected$prop_conv))
sum(Subset_Convergence_gamma_2st_corrected$prop_conv[Subset_Convergence_gamma_2st_corrected$KL_div==5], na.rm = TRUE) / sum(!is.na(Subset_Convergence_gamma_2st_corrected$prop_conv[Subset_Convergence_gamma_2st_corrected$KL_div==5]))
sum(Subset_Convergence_gamma_2st_corrected$prop_conv[Subset_Convergence_gamma_2st_corrected$KL_div==7], na.rm = TRUE) / sum(!is.na(Subset_Convergence_gamma_2st_corrected$prop_conv[Subset_Convergence_gamma_2st_corrected$KL_div==7]))



summary <- aggregate(Convergence_gamma_2st_corrected, by = list(Convergence_gamma_2st$n_dep, Convergence_gamma_2st$KL_div, Convergence_gamma_2st$n_t, Convergence_gamma_2st$n), FUN = mean, na.rm = TRUE)
View(summary)

sum_is.na <- function(x){
  sum(is.na(x))
}

summary_labelsw <- aggregate(Convergence_gamma_2st_corrected, by = list(Convergence_gamma_2st$n_dep, Convergence_gamma_2st$KL_div, Convergence_gamma_2st$n_t, Convergence_gamma_2st$n), FUN = sum_is.na)

View(summary_labelsw)

# 9 panels 
ggplot(data = summary[summary$valid_runs > 2,], mapping = aes(x = n_t, y = prop_conv, group = interaction(as.factor(n), as.factor(n_dep)), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  geom_line() +
  ylab("Prop converged transitions") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(n_dep)), cols = vars(as.factor(KL_div))) +
  # geom_hline(yintercept=1.1, linetype="dashed", color = "grey") + 
  xlab("number of observations per subject") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))


# 9 panels 
ggplot(data = summary, mapping = aes(x = n_t, y = GR_stat, group = interaction(as.factor(n), as.factor(n_dep)), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  geom_line() +
  ylab("GR_stat transitions") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(n_dep)), cols = vars(as.factor(KL_div))) +
  # geom_hline(yintercept=1.1, linetype="dashed", color = "grey") + 
  xlab("number of observations per subject") +
  theme_minimal() +
  geom_hline(yintercept = 1.1) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))




### OLD VERSION ####



Gelman_rubin_trans <- function(data1, data2, data3, m, J, burn_in, n_chain){
  L <- J - burn_in
  gamma_mean_data1 <- apply(data1[((burn_in + 1) : J),], 2, mean)
  gamma_mean_data2 <- apply(data2[((burn_in + 1) : J),], 2, mean)
  gamma_mean_data3 <- apply(data3[((burn_in + 1) : J),], 2, mean)
  gamma_grand_mean <- (gamma_mean_data1 + gamma_mean_data2 + gamma_mean_data3) / n_chain
  gamma_between_chain_var <- L / (n_chain-1) * ((gamma_mean_data1 - gamma_grand_mean)^2 + 
                                                  (gamma_mean_data2 - gamma_grand_mean)^2 +
                                                  (gamma_mean_data3 - gamma_grand_mean)^2)
  gamma_within_chain_var_d1 <- (1 / (n_chain-1)) * apply(data1[((burn_in + 1) : J),], 2, var)
  gamma_within_chain_var_d2 <- (1 / (n_chain-1)) * apply(data2[((burn_in + 1) : J),], 2, var)
  gamma_within_chain_var_d3 <- (1 / (n_chain-1)) * apply(data3[((burn_in + 1) : J),], 2, var)
  gamma_W <- (1/n_chain) * (gamma_within_chain_var_d1 + gamma_within_chain_var_d2 + gamma_within_chain_var_d3)
  
  gamma_GelRub <- sqrt(((L - 1)/ L * gamma_W + 1/L *  gamma_between_chain_var ) / gamma_W)
  return(gamma_GelRub)
}

sc <- 1
sim_it <- 1
n_sim <- 4

Convergence_gamma_2st <- data.frame(sim_rep = rep(1:4, each = true_m * true_m),
                                    true_m = true_m,
                                    n = rep(settings4[c(selected, selected_p2, selected_p4), 1], each = true_m * true_m * n_sim), 
                                    n_t = rep(settings4[c(selected, selected_p2, selected_p4), 2], each = true_m * true_m * n_sim),
                                    n_dep = rep(settings4[c(selected, selected_p2, selected_p4), 3], each = true_m * true_m * n_sim),
                                    KL_div = rep(settings4[c(selected, selected_p2, selected_p4), 4], each = true_m * true_m * n_sim), 
                                    from_state_i = rep(1:true_m, each = true_m),
                                    to_state_j = seq(1:true_m),
                                    GR_stat = NA) 

View(Convergence_gamma_2st)

for(sc in 1:171){
  for(sim_it in 1:4){
    Convergence_gamma_2st[(1 + (sim_it-1) * true_m * true_m + (sc-1) * true_m * true_m * n_sim):
                            ((sim_it) * true_m * true_m + (sc-1) * true_m * true_m * n_sim) ,9] <- Gelman_rubin_trans(data1 = gamma_full_conv0[[sc]][[sim_it]],  
                                                                                                                      data2 = gamma_full_conv1[[sc]][[sim_it]],  
                                                                                                                      data3 = gamma_full_conv2[[sc]][[sim_it]], 
                                                                                                                      m = true_m, 
                                                                                                                      J = J, 
                                                                                                                      burn_in = burn_in, 
                                                                                                                      n_chain = 3)
  }
}


View(Convergence_gamma_2st)

saveRDS(Convergence_gamma_2st, file = "result_tables/Convergence_gamma_2st.RDS")



Convergence_gamma_2st$prop <- (Convergence_gamma_2st$GR_stat > 1.1) * 1
summary <- aggregate(Convergence_gamma_2st, by = list(Convergence_gamma_2st$n_dep, Convergence_gamma_2st$KL_div, Convergence_gamma_2st$n_t, Convergence_gamma_2st$n), FUN = mean)
View(summary)

# 9 panels 
ggplot(data = summary, mapping = aes(x = n_t, y = prop, group = interaction(as.factor(n), as.factor(n_dep)), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  geom_line() +
  ylab("GR_stat") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep))) +
  # geom_hline(yintercept=1.1, linetype="dashed", color = "grey") + 
  xlab("number of observations per subject") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))



###################
## GR statistics emission distribution ####
##################

Gelman_rubin_emiss <- function(data1, data2, data3, m, J, burn_in, n_chain, n_dep){
  L <- J - burn_in
  emiss_mean_data1 <- emiss_mean_data2 <- emiss_mean_data3 <- emiss_grand_mean <- 
    emiss_between_chain_var <- emiss_within_chain_var_d1 <- emiss_within_chain_var_d2 <-
    emiss_within_chain_var_d3 <-emiss_W <- emiss_GelRub <- vector(mode = "list", n_dep)
  for(j in 1:n_dep){
    emiss_mean_data1[[j]] <- if(!is.null(data1)){apply(data1[[j]][((burn_in + 1) : J),], 2, mean)} else{NA}
    emiss_mean_data2[[j]] <- if(!is.null(data2)){apply(data2[[j]][((burn_in + 1) : J),], 2, mean)} else{NA}
    emiss_mean_data3[[j]] <- if(!is.null(data3)){apply(data3[[j]][((burn_in + 1) : J),], 2, mean)} else{NA}
    emiss_grand_mean[[j]] <- apply(cbind(emiss_mean_data1[[j]],emiss_mean_data2[[j]], emiss_mean_data3[[j]]),1, sum, na.rm = TRUE) / n_chain
    emiss_between_chain_var[[j]] <- L / (n_chain-1) * apply(cbind((emiss_mean_data1[[j]] - emiss_grand_mean[[j]])^2, 
                                                                  (emiss_mean_data2[[j]] - emiss_grand_mean[[j]])^2,
                                                                  (emiss_mean_data3[[j]] - emiss_grand_mean[[j]])^2), 1, sum, na.rm = TRUE)
    emiss_within_chain_var_d1[[j]] <-  if(!is.null(data1)){(1 / (n_chain-1)) * apply(data1[[j]][((burn_in + 1) : J),], 2, var)} else{NA}
    emiss_within_chain_var_d2[[j]] <-  if(!is.null(data2)){(1 / (n_chain-1)) * apply(data2[[j]][((burn_in + 1) : J),], 2, var)} else{NA}
    emiss_within_chain_var_d3[[j]] <-  if(!is.null(data3)){(1 / (n_chain-1)) * apply(data3[[j]][((burn_in + 1) : J),], 2, var)} else{NA}
    emiss_W[[j]] <- (1/n_chain) * apply(cbind(emiss_within_chain_var_d1[[j]], emiss_within_chain_var_d2[[j]], emiss_within_chain_var_d3[[j]]), 1, sum, na.rm = TRUE)
    
    emiss_GelRub[[j]] <- sqrt(((L - 1)/ L * emiss_W[[j]] + 1/L *  emiss_between_chain_var[[j]] ) / emiss_W[[j]])
  }
  return(emiss_GelRub)
}


Convergence_emission_2st_p1 <- data.frame(sim_rep = rep(rep(1:n_sim, times = length(c(settings4[c(selected), 3]))), times = rep(c(settings4[c(selected), 3] * true_m), each = 4)),
                                          true_m = true_m,
                                          n = rep(settings4[c(selected), 1], times = n_sim * true_m * settings4[c(selected), 3]), 
                                          n_t = rep(settings4[c(selected), 2], times = n_sim * true_m * settings4[c(selected), 3]),
                                          n_dep = rep(settings4[c(selected), 3], times = n_sim * true_m * settings4[c(selected), 3]),
                                          KL_div = rep(settings4[c(selected), 4], times = n_sim * true_m * settings4[c(selected), 3]), 
                                          dep = rep(sequence(rep(settings4[c(selected), 3], each = n_sim), 1, 1), each = true_m),
                                          k = 1:true_m,
                                          GR_stat = NA, 
                                          valid_runs = NA)

Convergence_emission_2st_p2 <- data.frame(sim_rep = rep(rep(1:n_sim, times = length(c(settings4[c(selected_p2), 3]))), times = rep(c(settings4[c(selected_p2), 3] * true_m), each = 4)),
                                          true_m = true_m,
                                          n = rep(settings4[c(selected_p2), 1], times = n_sim * true_m * settings4[c(selected_p2), 3]), 
                                          n_t = rep(settings4[c(selected_p2), 2], times = n_sim * true_m * settings4[c(selected_p2), 3]),
                                          n_dep = rep(settings4[c(selected_p2), 3], times = n_sim * true_m * settings4[c(selected_p2), 3]),
                                          KL_div = rep(settings4[c(selected_p2), 4], times = n_sim * true_m * settings4[c(selected_p2), 3]), 
                                          dep = rep(sequence(rep(settings4[c(selected_p2), 3], each = n_sim), 1, 1), each = true_m),
                                          k = 1:true_m,
                                          GR_stat = NA, 
                                          valid_runs = NA
)
Convergence_emission_2st_p4 <- data.frame(sim_rep = rep(rep(1:n_sim, times = length(c(settings4[c(selected_p4), 3]))), times = rep(c(settings4[c(selected_p4), 3] * true_m), each = 4)),
                                          true_m = true_m,
                                          n = rep(settings4[c(selected_p4), 1], times = n_sim * true_m * settings4[c(selected_p4), 3]), 
                                          n_t = rep(settings4[c(selected_p4), 2], times = n_sim * true_m * settings4[c(selected_p4), 3]),
                                          n_dep = rep(settings4[c(selected_p4), 3], times = n_sim * true_m * settings4[c(selected_p4), 3]),
                                          KL_div = rep(settings4[c(selected_p4), 4], times = n_sim * true_m * settings4[c(selected_p4), 3]), 
                                          dep = rep(sequence(rep(settings4[c(selected_p4), 3], each = n_sim), 1, 1), each = true_m),
                                          k = 1:true_m,
                                          GR_stat = NA, 
                                          valid_runs = NA
)

emiss_mean_full_conv0 <- extracted_convergence_2st$emiss_mean_full_conv0
emiss_mean_full_conv1 <- extracted_convergence_2st$emiss_mean_full_conv1
emiss_mean_full_conv2 <- extracted_convergence_2st$emiss_mean_full_conv2


it <- n_sim * true_m * settings4[c(selected), 3]
cumsum_it <- c(0,cumsum(it))

for(sc in 1:length(selected)){
  for(sim_it in 1:4){
    spacer <- true_m * settings4[c(selected), 3][sc]
    ind <- (1 + cumsum_it[sc] + (sim_it - 1) * spacer)
    sum_present <- aggr_Label_switch_proxy_2st1p1$presentSUM[which(aggr_Label_switch_proxy_2st1p1$sim_iteration == sim_it &
                                                                     aggr_Label_switch_proxy_2st1p1$n == Convergence_emission_2st_p1$n[ind] &
                                                                     aggr_Label_switch_proxy_2st1p1$n_t == Convergence_emission_2st_p1$n_t[ind] &
                                                                     aggr_Label_switch_proxy_2st1p1$n_dep == Convergence_emission_2st_p1$n_dep[ind] &
                                                                     aggr_Label_switch_proxy_2st1p1$KL_div == Convergence_emission_2st_p1$KL_div[ind])]
    
    valid <- 4- aggr2_Label_switch_proxy_2st1p1$no_chain[which(aggr2_Label_switch_proxy_2st1p1$Group.4 == Convergence_emission_2st_p1$n[ind] &
                                                                 aggr2_Label_switch_proxy_2st1p1$Group.3 == Convergence_emission_2st_p1$n_t[ind] &
                                                                 aggr2_Label_switch_proxy_2st1p1$Group.2 == Convergence_emission_2st_p1$n_dep[ind] &
                                                                 aggr2_Label_switch_proxy_2st1p1$Group.1 == Convergence_emission_2st_p1$KL_div[ind])]
    Convergence_emission_2st_p1[ind: 
                                  (cumsum_it[sc] + sim_it * spacer), 10] <- valid
    if(sum_present == 0){
      Convergence_emission_2st_p1[ind: 
                                    (cumsum_it[sc] + sim_it * spacer), 9] <- unlist(Gelman_rubin_emiss(data1 = emiss_mean_full_conv0[[sc]][[sim_it]],  
                                                                                                       data2 = emiss_mean_full_conv1[[sc]][[sim_it]],  
                                                                                                       data3 = emiss_mean_full_conv2[[sc]][[sim_it]], 
                                                                                                       m = true_m, 
                                                                                                       J = J, 
                                                                                                       burn_in = burn_in, 
                                                                                                       n_chain = 3, 
                                                                                                       n_dep = settings4[c(selected), 3][sc]))
      
    } 
  }
}

# View(Convergence_emission_2st_p1)


it <- n_sim * true_m * settings4[c(selected_p2), 3]
cumsum_it <- c(0,cumsum(it))

for(sc in 1:length(selected_p2)){
  for(sim_it in 1:4){
    spacer <-true_m * settings4[c(selected_p2), 3][sc]
    ind <- (1 + cumsum_it[sc] + (sim_it - 1) * spacer)
    sum_present <- aggr_Label_switch_proxy_2st1p2$presentSUM[which(aggr_Label_switch_proxy_2st1p2$sim_iteration == sim_it &
                                                                     aggr_Label_switch_proxy_2st1p2$n == Convergence_emission_2st_p2$n[ind] &
                                                                     aggr_Label_switch_proxy_2st1p2$n_t == Convergence_emission_2st_p2$n_t[ind] &
                                                                     aggr_Label_switch_proxy_2st1p2$n_dep == Convergence_emission_2st_p2$n_dep[ind] &
                                                                     aggr_Label_switch_proxy_2st1p2$KL_div == Convergence_emission_2st_p2$KL_div[ind])]
    valid <- 4- aggr2_Label_switch_proxy_2st1p2$no_chain[which(aggr2_Label_switch_proxy_2st1p2$Group.4 == Convergence_emission_2st_p2$n[ind] &
                                                                 aggr2_Label_switch_proxy_2st1p2$Group.3 == Convergence_emission_2st_p2$n_t[ind] &
                                                                 aggr2_Label_switch_proxy_2st1p2$Group.2 == Convergence_emission_2st_p2$n_dep[ind] &
                                                                 aggr2_Label_switch_proxy_2st1p2$Group.1 == Convergence_emission_2st_p2$KL_div[ind])]
    Convergence_emission_2st_p2[ind: 
                                  (cumsum_it[sc] + sim_it * spacer), 10] <- valid
    
    if(sum_present == 0){
      Convergence_emission_2st_p2[ind: 
                                    (cumsum_it[sc] + sim_it * spacer), 9] <- unlist(Gelman_rubin_emiss(data1 = emiss_mean_full_conv0[[sc + l_files]][[sim_it]],  
                                                                                                       data2 = emiss_mean_full_conv1[[sc + l_files]][[sim_it]],  
                                                                                                       data3 = emiss_mean_full_conv2[[sc + l_files]][[sim_it]], 
                                                                                                       m = true_m, 
                                                                                                       J = J, 
                                                                                                       burn_in = burn_in, 
                                                                                                       n_chain = 3, 
                                                                                                       n_dep = settings4[c(selected_p2), 3][sc]))
      
    }
  }
}

# View(Convergence_emission_2st_p2)

it <- n_sim * true_m * settings4[c(selected_p4), 3]
cumsum_it <- c(0,cumsum(it))

for(sc in 1:length(selected_p4)){
  for(sim_it in 1:4){
    spacer <-true_m * settings4[c(selected_p4), 3][sc]
    ind <- (1 + cumsum_it[sc] + (sim_it - 1) * spacer)
    sum_present <- aggr_Label_switch_proxy_2st1p4$presentSUM[which(aggr_Label_switch_proxy_2st1p4$sim_iteration == sim_it &
                                                                     aggr_Label_switch_proxy_2st1p4$n == Convergence_emission_2st_p4$n[ind] &
                                                                     aggr_Label_switch_proxy_2st1p4$n_t == Convergence_emission_2st_p4$n_t[ind] &
                                                                     aggr_Label_switch_proxy_2st1p4$n_dep == Convergence_emission_2st_p4$n_dep[ind] &
                                                                     aggr_Label_switch_proxy_2st1p4$KL_div == Convergence_emission_2st_p4$KL_div[ind])]
    valid <- 4- aggr2_Label_switch_proxy_2st1p4$no_chain[which(aggr2_Label_switch_proxy_2st1p4$Group.4 == Convergence_emission_2st_p4$n[ind] &
                                                                 aggr2_Label_switch_proxy_2st1p4$Group.3 == Convergence_emission_2st_p4$n_t[ind] &
                                                                 aggr2_Label_switch_proxy_2st1p4$Group.2 == Convergence_emission_2st_p4$n_dep[ind] &
                                                                 aggr2_Label_switch_proxy_2st1p4$Group.1 == Convergence_emission_2st_p4$KL_div[ind])]
    Convergence_emission_2st_p4[ind: 
                                  (cumsum_it[sc] + sim_it * spacer), 10] <- valid
    if(sum_present == 0){
      Convergence_emission_2st_p4[ind: 
                                    (cumsum_it[sc] + sim_it * spacer), 9] <- unlist(Gelman_rubin_emiss(data1 = emiss_mean_full_conv0[[sc + l_files_p2]][[sim_it]],  
                                                                                                       data2 = emiss_mean_full_conv1[[sc + l_files_p2]][[sim_it]],  
                                                                                                       data3 = emiss_mean_full_conv2[[sc + l_files_p2]][[sim_it]], 
                                                                                                       m = true_m, 
                                                                                                       J = J, 
                                                                                                       burn_in = burn_in, 
                                                                                                       n_chain = 3, 
                                                                                                       n_dep = settings4[c(selected_p4), 3][sc]))
      
    }
  }
}

# View(Convergence_emission_2st_p4)

### merging transition convergence results over partitions ####
Convergence_emission_2st_corrected <- rbind(Convergence_emission_2st_p1,
                                            Convergence_emission_2st_p2,
                                            Convergence_emission_2st_p4)

saveRDS(Convergence_emission_2st_corrected, file = "result_tables/Convergence_emission_2st_corrected_V2.RDS")


Convergence_emission_2st_corrected$prop <- (Convergence_emission_2st_corrected$GR_stat > 1.1) * 1
Convergence_emission_2st_corrected$prop_conv <- (Convergence_emission_2st_corrected$GR_stat < 1.1) * 1

sum(Convergence_emission_2st_corrected$prop_conv, na.rm = TRUE) / sum(!is.na(Convergence_emission_2st_corrected$prop_conv))

sum(Convergence_emission_2st_corrected$prop_conv[Convergence_emission_2st_corrected$KL_div == 3], na.rm = TRUE) / sum(!is.na(Convergence_emission_2st_corrected$prop_conv[Convergence_emission_2st_corrected$KL_div == 3]))
sum(Convergence_emission_2st_corrected$prop_conv[Convergence_emission_2st_corrected$KL_div == 5], na.rm = TRUE) / sum(!is.na(Convergence_emission_2st_corrected$prop_conv[Convergence_emission_2st_corrected$KL_div == 5]))
sum(Convergence_emission_2st_corrected$prop_conv[Convergence_emission_2st_corrected$KL_div == 7], na.rm = TRUE) / sum(!is.na(Convergence_emission_2st_corrected$prop_conv[Convergence_emission_2st_corrected$KL_div == 7]))


Subset_Convergence_emission_2st_corrected <- Convergence_emission_2st_corrected[Convergence_emission_2st_corrected$KL_div > 3 &
                                                                                  Convergence_emission_2st_corrected$valid_runs > 2,]

sum(Subset_Convergence_emission_2st_corrected$prop_conv, na.rm = TRUE) /sum(!is.na(Subset_Convergence_emission_2st_corrected$prop_conv))
sum(Subset_Convergence_emission_2st_corrected$prop_conv[Subset_Convergence_emission_2st_corrected$KL_div==5], na.rm = TRUE) / sum(!is.na(Subset_Convergence_emission_2st_corrected$prop_conv[Subset_Convergence_emission_2st_corrected$KL_div==5]))
sum(Subset_Convergence_emission_2st_corrected$prop_conv[Subset_Convergence_emission_2st_corrected$KL_div==7], na.rm = TRUE) / sum(!is.na(Subset_Convergence_emission_2st_corrected$prop_conv[Subset_Convergence_emission_2st_corrected$KL_div==7]))



summary_emiss <- aggregate(Convergence_emission_2st_corrected, by = list(Convergence_emission_2st_corrected$n_dep, 
                                                                         Convergence_emission_2st_corrected$KL_div, 
                                                                         Convergence_emission_2st_corrected$n_t, 
                                                                         Convergence_emission_2st_corrected$n), FUN = mean, na.rm = TRUE)

summary_emiss_prop_conv <- aggregate(Convergence_emission_2st_corrected, by = list(Convergence_emission_2st_corrected$n_dep, 
                                                                                   Convergence_emission_2st_corrected$KL_div, 
                                                                                   Convergence_emission_2st_corrected$n_t, 
                                                                                   Convergence_emission_2st_corrected$n), FUN = sum_is.na)
View(summary_emiss_prop_conv)


# 9 panels 
ggplot(data = summary_emiss[summary_emiss$valid_runs>1,], mapping = aes(x = n_t, y = GR_stat, group = interaction(as.factor(n), as.factor(n_dep)), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  geom_line() +
  ylab("GR_stat") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(n_dep)), cols = vars(as.factor(KL_div))) +
  # geom_hline(yintercept=1.1, linetype="dashed", color = "grey") + 
  xlab("number of observations per subject") +
  theme_minimal() +
  geom_hline(yintercept = 1.1) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

# 9 panels 
ggplot(data = summary_emiss[summary_emiss$valid_runs>1,], mapping = aes(x = n_t, y = prop_conv, group = interaction(as.factor(n), as.factor(n_dep)), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  geom_line() +
  ylab("Proportion converged") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(n_dep)), cols = vars(as.factor(KL_div))) +
  # geom_hline(yintercept=1.1, linetype="dashed", color = "grey") + 
  xlab("number of observations per subject") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))



### OLD VERSION ####
Gelman_rubin_emiss <- function(data1, data2, data3, m, J, burn_in, n_chain, n_dep){
  L <- J - burn_in
  emiss_mean_data1 <- emiss_mean_data2 <- emiss_mean_data3 <- emiss_grand_mean <- 
    emiss_between_chain_var <- emiss_within_chain_var_d1 <- emiss_within_chain_var_d2 <-
    emiss_within_chain_var_d3 <-emiss_W <- emiss_GelRub <- vector(mode = "list", n_dep)
  for(j in 1:n_dep){
    emiss_mean_data1[[j]] <- apply(data1[[j]][((burn_in + 1) : J),], 2, mean)
    emiss_mean_data2[[j]] <- apply(data2[[j]][((burn_in + 1) : J),], 2, mean)
    emiss_mean_data3[[j]] <- apply(data3[[j]][((burn_in + 1) : J),], 2, mean)
    emiss_grand_mean[[j]] <- (emiss_mean_data1[[j]] + emiss_mean_data2[[j]] + emiss_mean_data3[[j]]) / n_chain
    emiss_between_chain_var[[j]] <- L / (n_chain-1) * ((emiss_mean_data1[[j]] - emiss_grand_mean[[j]])^2 + 
                                                         (emiss_mean_data2[[j]] - emiss_grand_mean[[j]])^2 +
                                                         (emiss_mean_data3[[j]] - emiss_grand_mean[[j]])^2)
    emiss_within_chain_var_d1[[j]] <- (1 / (n_chain-1)) * apply(data1[[j]][((burn_in + 1) : J),], 2, var)
    emiss_within_chain_var_d2[[j]] <- (1 / (n_chain-1)) * apply(data2[[j]][((burn_in + 1) : J),], 2, var)
    emiss_within_chain_var_d3[[j]] <- (1 / (n_chain-1)) * apply(data3[[j]][((burn_in + 1) : J),], 2, var)
    emiss_W[[j]] <- (1/n_chain) * (emiss_within_chain_var_d1[[j]] + emiss_within_chain_var_d2[[j]] + emiss_within_chain_var_d3[[j]])
    
    emiss_GelRub[[j]] <- sqrt(((L - 1)/ L * emiss_W[[j]] + 1/L *  emiss_between_chain_var[[j]] ) / emiss_W[[j]])
  }
  return(emiss_GelRub)
}


Convergence_emission_2st <- data.frame(sim_rep = rep(rep(1:n_sim, times = length(c(settings4[c(selected, selected_p2, selected_p4), 3]))), times = rep(c(settings4[c(selected, selected_p2, selected_p4), 3] * true_m), each = 4)),
                                       true_m = true_m,
                                       n = rep(settings4[c(selected, selected_p2, selected_p4), 1], times = n_sim * true_m * settings4[c(selected, selected_p2, selected_p4), 3]), 
                                       n_t = rep(settings4[c(selected, selected_p2, selected_p4), 2], times = n_sim * true_m * settings4[c(selected, selected_p2, selected_p4), 3]),
                                       n_dep = rep(settings4[c(selected, selected_p2, selected_p4), 3], times = n_sim * true_m * settings4[c(selected, selected_p2, selected_p4), 3]),
                                       KL_div = rep(settings4[c(selected, selected_p2, selected_p4), 4], times = n_sim * true_m * settings4[c(selected, selected_p2, selected_p4), 3]), 
                                       dep = rep(sequence(rep(settings4[c(selected, selected_p2, selected_p4), 3], each = n_sim), 1, 1), each = true_m),
                                       k = 1:true_m,
                                       GR_stat = NA
)

View(Convergence_emission_2st)

it <- n_sim * true_m * settings4[c(selected, selected_p2, selected_p4), 3]
cumsum_it <- c(0,cumsum(it))

for(sc in 1:171){
  for(sim_it in 1:4){
    Convergence_emission_2st[(1 + cumsum_it[sc] + (sim_it - 1) * true_m * settings4[c(selected, selected_p2, selected_p4), 3][sc]): 
                               (cumsum_it[sc] + sim_it * true_m * settings4[c(selected, selected_p2, selected_p4), 3][sc]),9] <- unlist(Gelman_rubin_emiss(data1 = emiss_mean_full_conv0[[sc]][[sim_it]],  
                                                                                                                                                           data2 = emiss_mean_full_conv1[[sc]][[sim_it]],  
                                                                                                                                                           data3 = emiss_mean_full_conv2[[sc]][[sim_it]], 
                                                                                                                                                           m = true_m, 
                                                                                                                                                           J = J, 
                                                                                                                                                           burn_in = burn_in, 
                                                                                                                                                           n_chain = 3, 
                                                                                                                                                           n_dep = settings4[c(selected, selected_p2, selected_p4), 3][sc]))
  }
}


View(Convergence_emission_2st)

saveRDS(Convergence_emission_2st, file = "result_tables/Convergence_emission_2st.RDS")




Convergence_emission_2st$prop <- (Convergence_emission_2st$GR_stat > 1.1) * 1
summary_emiss <- aggregate(Convergence_emission_2st, by = list(Convergence_emission_2st$n_dep, Convergence_emission_2st$KL_div, Convergence_emission_2st$n_t, Convergence_emission_2st$n), FUN = mean)
View(summary_emiss)




# 9 panels 
ggplot(data = summary_emiss, mapping = aes(x = n_t, y = prop, group = interaction(as.factor(n), as.factor(n_dep)), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  geom_line() +
  ylab("GR_stat") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep))) +
  # geom_hline(yintercept=1.1, linetype="dashed", color = "grey") + 
  xlab("number of observations per subject") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))
