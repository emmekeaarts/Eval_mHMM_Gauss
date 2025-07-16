library(mHMMbayes)
library(ggplot2)
library(reshape2)

# Simulstion settings
settings4 <- read.csv("sim_scenarios_V4.csv")
true_m  <- 1

# scenario <- 1
# KL_div  <- settings[scenario, 4]
# n_t     <- settings[scenario, 2]
# n       <- settings[scenario, 1]
# n_dep   <- settings[scenario, 3]
# out_file <- paste0("out_n", n, "_nt", n_t, "_ndep", n_dep, "_", true_m, "st_KLD", KL_div)

## Reading in output files ####
# 1 states
# runs that could be done in 1 go
selected <- sort(c(1:120))
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


selected_p4 <- sort(c(145:171))
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

group_out_1st_emiss_mean <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
group_out_1st_emiss_sd <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
means_true <-vector("list", length = length(c(files, files_p2_1, files_p4_1)))
selection_out <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))

names(group_out_1st_emiss_mean) <- c(files, files_p2_1, files_p4_1)
names(group_out_1st_emiss_sd) <-c(files, files_p2_1, files_p4_1)
names(means_true) <- c(files, files_p2_1, files_p4_1)
names(selection_out) <- c(files, files_p2_1, files_p4_1)

for(sc in 1:length(files)){
  out[[sc]] <- readRDS(paste0("out_V4/", files[sc], ".rds")) 
  selection_out[[sc]]             <- out[[sc]]$selection_out
  group_out_1st_emiss_mean[[sc]]  <- out[[sc]]$group_out_1st$emiss_mean
  group_out_1st_emiss_sd[[sc]]    <- out[[sc]]$group_out_1st$emiss_sd
  means_true[[sc]]                <- out[[sc]]$means_true
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
  group_out_1st_emiss_mean[[sc + l_files]]$median     <-  Map(rbind,out_p2_1[[sc]]$group_out_1st$emiss_mean$median, out_p2_2[[sc]]$group_out_1st$emiss_mean$median)
  group_out_1st_emiss_mean[[sc + l_files]]$quant_2.5  <-  Map(rbind,out_p2_1[[sc]]$group_out_1st$emiss_mean$quant_2.5, out_p2_2[[sc]]$group_out_1st$emiss_mean$quant_2.5)
  group_out_1st_emiss_mean[[sc + l_files]]$quant_5    <-  Map(rbind,out_p2_1[[sc]]$group_out_1st$emiss_mean$quant_5, out_p2_2[[sc]]$group_out_1st$emiss_mean$quant_5)
  group_out_1st_emiss_mean[[sc + l_files]]$quant_95   <-  Map(rbind,out_p2_1[[sc]]$group_out_1st$emiss_mean$quant_95, out_p2_2[[sc]]$group_out_1st$emiss_mean$quant_95)
  group_out_1st_emiss_mean[[sc + l_files]]$quant_97.5 <-  Map(rbind,out_p2_1[[sc]]$group_out_1st$emiss_mean$quant_97.5, out_p2_2[[sc]]$group_out_1st$emiss_mean$quant_97.5)
  group_out_1st_emiss_sd[[sc + l_files]]$median       <-  Map(rbind,out_p2_1[[sc]]$group_out_1st$emiss_sd$median, out_p2_2[[sc]]$group_out_1st$emiss_sd$median)
  group_out_1st_emiss_sd[[sc + l_files]]$quant_2.5    <-  Map(rbind,out_p2_1[[sc]]$group_out_1st$emiss_sd$quant_2.5, out_p2_2[[sc]]$group_out_1st$emiss_sd$quant_2.5)
  group_out_1st_emiss_sd[[sc + l_files]]$quant_5      <-  Map(rbind,out_p2_1[[sc]]$group_out_1st$emiss_sd$quant_5, out_p2_2[[sc]]$group_out_1st$emiss_sd$quant_5)
  group_out_1st_emiss_sd[[sc + l_files]]$quant_95     <-  Map(rbind,out_p2_1[[sc]]$group_out_1st$emiss_sd$quant_95, out_p2_2[[sc]]$group_out_1st$emiss_sd$quant_95)
  group_out_1st_emiss_sd[[sc + l_files]]$quant_97.5   <-  Map(rbind,out_p2_1[[sc]]$group_out_1st$emiss_sd$quant_97.5, out_p2_2[[sc]]$group_out_1st$emiss_sd$quant_97.5)
  means_true[[sc + l_files]]                          <- Map(rbind, out_p2_1[[sc]]$means_true, out_p2_2[[sc]]$means_true)
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
  group_out_1st_emiss_mean[[sc + l_files_p2]]$median     <-  Map(rbind,out_p4_1[[sc]]$group_out_1st$emiss_mean$median, 
                                                                 out_p4_2[[sc]]$group_out_1st$emiss_mean$median,
                                                                 out_p4_3[[sc]]$group_out_1st$emiss_mean$median,
                                                                 out_p4_4[[sc]]$group_out_1st$emiss_mean$median)
  group_out_1st_emiss_mean[[sc + l_files_p2]]$quant_2.5  <-  Map(rbind,out_p4_1[[sc]]$group_out_1st$emiss_mean$quant_2.5, 
                                                                 out_p4_2[[sc]]$group_out_1st$emiss_mean$quant_2.5,
                                                                 out_p4_3[[sc]]$group_out_1st$emiss_mean$quant_2.5,
                                                                 out_p4_4[[sc]]$group_out_1st$emiss_mean$quant_2.5)
  group_out_1st_emiss_mean[[sc + l_files_p2]]$quant_5    <-  Map(rbind,out_p4_1[[sc]]$group_out_1st$emiss_mean$quant_5, 
                                                                 out_p4_2[[sc]]$group_out_1st$emiss_mean$quant_5,
                                                                 out_p4_3[[sc]]$group_out_1st$emiss_mean$quant_5, 
                                                                 out_p4_4[[sc]]$group_out_1st$emiss_mean$quant_5)
  group_out_1st_emiss_mean[[sc + l_files_p2]]$quant_95   <-  Map(rbind,out_p4_1[[sc]]$group_out_1st$emiss_mean$quant_95, 
                                                                 out_p4_2[[sc]]$group_out_1st$emiss_mean$quant_95,
                                                                 out_p4_3[[sc]]$group_out_1st$emiss_mean$quant_95, 
                                                                 out_p4_4[[sc]]$group_out_1st$emiss_mean$quant_95)
  group_out_1st_emiss_mean[[sc + l_files_p2]]$quant_97.5 <-  Map(rbind,out_p4_1[[sc]]$group_out_1st$emiss_mean$quant_97.5, 
                                                                 out_p4_2[[sc]]$group_out_1st$emiss_mean$quant_97.5,
                                                                 out_p4_3[[sc]]$group_out_1st$emiss_mean$quant_97.5, 
                                                                 out_p4_4[[sc]]$group_out_1st$emiss_mean$quant_97.5)
  group_out_1st_emiss_sd[[sc + l_files_p2]]$median       <-  Map(rbind,out_p4_1[[sc]]$group_out_1st$emiss_sd$median, 
                                                                 out_p4_2[[sc]]$group_out_1st$emiss_sd$median,
                                                                 out_p4_3[[sc]]$group_out_1st$emiss_sd$median, 
                                                                 out_p4_4[[sc]]$group_out_1st$emiss_sd$median)
  group_out_1st_emiss_sd[[sc + l_files_p2]]$quant_2.5    <-  Map(rbind,out_p4_1[[sc]]$group_out_1st$emiss_sd$quant_2.5, 
                                                                 out_p4_2[[sc]]$group_out_1st$emiss_sd$quant_2.5,
                                                                 out_p4_3[[sc]]$group_out_1st$emiss_sd$quant_2.5, 
                                                                 out_p4_4[[sc]]$group_out_1st$emiss_sd$quant_2.5)
  group_out_1st_emiss_sd[[sc + l_files_p2]]$quant_5      <-  Map(rbind,out_p4_1[[sc]]$group_out_1st$emiss_sd$quant_5, 
                                                                 out_p4_2[[sc]]$group_out_1st$emiss_sd$quant_5,
                                                                 out_p4_3[[sc]]$group_out_1st$emiss_sd$quant_5, 
                                                                 out_p4_4[[sc]]$group_out_1st$emiss_sd$quant_5)
  group_out_1st_emiss_sd[[sc + l_files_p2]]$quant_95     <-  Map(rbind,out_p4_1[[sc]]$group_out_1st$emiss_sd$quant_95, 
                                                                 out_p4_2[[sc]]$group_out_1st$emiss_sd$quant_95,
                                                                 out_p4_3[[sc]]$group_out_1st$emiss_sd$quant_95, 
                                                                 out_p4_4[[sc]]$group_out_1st$emiss_sd$quant_95)
  group_out_1st_emiss_sd[[sc + l_files_p2]]$quant_97.5   <-  Map(rbind,out_p4_1[[sc]]$group_out_1st$emiss_sd$quant_97.5, 
                                                                 out_p4_2[[sc]]$group_out_1st$emiss_sd$quant_97.5,
                                                                 out_p4_3[[sc]]$group_out_1st$emiss_sd$quant_97.5, 
                                                                 out_p4_4[[sc]]$group_out_1st$emiss_sd$quant_97.5)
  means_true[[sc + l_files_p2]]                          <- Map(rbind, out_p4_1[[sc]]$means_true, 
                                                                out_p4_2[[sc]]$means_true,
                                                                out_p4_3[[sc]]$means_true, 
                                                                out_p4_4[[sc]]$means_true)
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

extracted_results_1st <- list(group_out_1st_emiss_mean = group_out_1st_emiss_mean,
                              group_out_1st_emiss_sd = group_out_1st_emiss_sd,
                              means_true = means_true,
                              selection_out = selection_out)
saveRDS(extracted_results_1st, file = "extracted_results_1st.RDS")


# Model selection ####

Model_selection_1st <- data.frame(sim_iteration = 1:n_sim,
                                  true_m = 1,
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
    Model_selection_1st[(1 + (sc - 1) * n_sim * 4 + n_sim * (kc - 1)) : 
                          (100 + (sc - 1) * n_sim * 4 + n_sim * (kc - 1)) ,8] <- selection_out[[sc]]$AIC_mean[,kc]
    Model_selection_1st[(1 + (sc - 1) * n_sim * 4 + n_sim * (kc - 1)) : 
                          (100 + (sc - 1) * n_sim * 4 + n_sim * (kc - 1)) ,9] <- selection_out[[sc]]$AICc_mean[,kc]
  }
}
View(Model_selection_1st)
saveRDS(Model_selection_1st, file = "result_tables/Model_selection_1st.RDS")
