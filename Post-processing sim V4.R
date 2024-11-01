library(mHMMbayes)
library(ggplot2)
library(reshape2)

# Simulstion settings
settings4 <- read.csv("sim_scenarios_V4.csv")
true_m  <- 3

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

group_out_3st_emiss_mean <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
group_out_3st_emiss_sd <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
group_out_3st_gamma_prob <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
group_out_3st_gamma_int <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
means_true <-vector("list", length = length(c(files, files_p2_1, files_p4_1)))
selection_out <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
inferred_states <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
true_state_probs <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))
true_states <- vector("list", length = length(c(files, files_p2_1, files_p4_1)))

names(group_out_3st_emiss_mean) <- c(files, files_p2_1, files_p4_1)
names(group_out_3st_emiss_sd) <-c(files, files_p2_1, files_p4_1)
names(group_out_3st_gamma_prob) <- c(files, files_p2_1, files_p4_1)
names(group_out_3st_gamma_int) <- c(files, files_p2_1, files_p4_1)
names(means_true) <- c(files, files_p2_1, files_p4_1)
names(selection_out) <- c(files, files_p2_1, files_p4_1)
names(inferred_states) <- c(files, files_p2_1, files_p4_1)
names(true_state_probs) <- c(files, files_p2_1, files_p4_1)
names(true_states) <- c(files, files_p2_1, files_p4_1)

for(sc in 1:length(files)){
  out[[sc]] <- readRDS(paste0("out_V4/", files[sc], ".rds")) 
  selection_out[[sc]]             <- out[[sc]]$selection_out
  group_out_3st_emiss_mean[[sc]]  <- out[[sc]]$group_out_3st$emiss_mean
  group_out_3st_emiss_sd[[sc]]    <- out[[sc]]$group_out_3st$emiss_sd
  group_out_3st_gamma_prob[[sc]]  <- out[[sc]]$group_out_3st$trans_prob
  group_out_3st_gamma_int[[sc]]   <- out[[sc]]$group_out_3st$trans_interc
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
  group_out_3st_emiss_mean[[sc + l_files]]$median     <-  Map(rbind,out_p2_1[[sc]]$group_out_3st$emiss_mean$median, out_p2_2[[sc]]$group_out_3st$emiss_mean$median)
  group_out_3st_emiss_mean[[sc + l_files]]$quant_2.5  <-  Map(rbind,out_p2_1[[sc]]$group_out_3st$emiss_mean$quant_2.5, out_p2_2[[sc]]$group_out_3st$emiss_mean$quant_2.5)
  group_out_3st_emiss_mean[[sc + l_files]]$quant_5    <-  Map(rbind,out_p2_1[[sc]]$group_out_3st$emiss_mean$quant_5, out_p2_2[[sc]]$group_out_3st$emiss_mean$quant_5)
  group_out_3st_emiss_mean[[sc + l_files]]$quant_95   <-  Map(rbind,out_p2_1[[sc]]$group_out_3st$emiss_mean$quant_95, out_p2_2[[sc]]$group_out_3st$emiss_mean$quant_95)
  group_out_3st_emiss_mean[[sc + l_files]]$quant_97.5 <-  Map(rbind,out_p2_1[[sc]]$group_out_3st$emiss_mean$quant_97.5, out_p2_2[[sc]]$group_out_3st$emiss_mean$quant_97.5)
  group_out_3st_emiss_sd[[sc + l_files]]$median       <-  Map(rbind,out_p2_1[[sc]]$group_out_3st$emiss_sd$median, out_p2_2[[sc]]$group_out_3st$emiss_sd$median)
  group_out_3st_emiss_sd[[sc + l_files]]$quant_2.5    <-  Map(rbind,out_p2_1[[sc]]$group_out_3st$emiss_sd$quant_2.5, out_p2_2[[sc]]$group_out_3st$emiss_sd$quant_2.5)
  group_out_3st_emiss_sd[[sc + l_files]]$quant_5      <-  Map(rbind,out_p2_1[[sc]]$group_out_3st$emiss_sd$quant_5, out_p2_2[[sc]]$group_out_3st$emiss_sd$quant_5)
  group_out_3st_emiss_sd[[sc + l_files]]$quant_95     <-  Map(rbind,out_p2_1[[sc]]$group_out_3st$emiss_sd$quant_95, out_p2_2[[sc]]$group_out_3st$emiss_sd$quant_95)
  group_out_3st_emiss_sd[[sc + l_files]]$quant_97.5   <-  Map(rbind,out_p2_1[[sc]]$group_out_3st$emiss_sd$quant_97.5, out_p2_2[[sc]]$group_out_3st$emiss_sd$quant_97.5)
  group_out_3st_gamma_prob[[sc + l_files]]            <-  Map(rbind, out_p2_1[[sc]]$group_out_3st$trans_prob, out_p2_2[[sc]]$group_out_3st$trans_prob)
  group_out_3st_gamma_int[[sc + l_files]]             <-  Map(rbind, out_p2_1[[sc]]$group_out_3st$trans_interc, out_p2_2[[sc]]$group_out_3st$trans_interc)
  means_true[[sc + l_files]]                          <- Map(rbind, out_p2_1[[sc]]$means_true, out_p2_2[[sc]]$means_true)
  inferred_states[[sc + l_files]]                     <- cbind(out_p2_1[[sc]]$inferred_states, out_p2_2[[sc]]$inferred_states[,-1])
  true_state_probs[[sc + l_files]]                    <- cbind(out_p2_1[[sc]]$true_state_probs, out_p2_2[[sc]]$true_state_probs[,-1])
  true_states[[sc + l_files]]                         <- cbind(out_p2_1[[sc]]$true_states, out_p2_2[[sc]]$true_states[,-1])
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
  group_out_3st_emiss_mean[[sc + l_files_p2]]$median     <-  Map(rbind,out_p4_1[[sc]]$group_out_3st$emiss_mean$median, 
                                                                    out_p4_2[[sc]]$group_out_3st$emiss_mean$median,
                                                                    out_p4_3[[sc]]$group_out_3st$emiss_mean$median,
                                                                    out_p4_4[[sc]]$group_out_3st$emiss_mean$median)
  group_out_3st_emiss_mean[[sc + l_files_p2]]$quant_2.5  <-  Map(rbind,out_p4_1[[sc]]$group_out_3st$emiss_mean$quant_2.5, 
                                                                    out_p4_2[[sc]]$group_out_3st$emiss_mean$quant_2.5,
                                                                    out_p4_3[[sc]]$group_out_3st$emiss_mean$quant_2.5,
                                                                    out_p4_4[[sc]]$group_out_3st$emiss_mean$quant_2.5)
  group_out_3st_emiss_mean[[sc + l_files_p2]]$quant_5    <-  Map(rbind,out_p4_1[[sc]]$group_out_3st$emiss_mean$quant_5, 
                                                                    out_p4_2[[sc]]$group_out_3st$emiss_mean$quant_5,
                                                                    out_p4_3[[sc]]$group_out_3st$emiss_mean$quant_5, 
                                                                    out_p4_4[[sc]]$group_out_3st$emiss_mean$quant_5)
  group_out_3st_emiss_mean[[sc + l_files_p2]]$quant_95   <-  Map(rbind,out_p4_1[[sc]]$group_out_3st$emiss_mean$quant_95, 
                                                                    out_p4_2[[sc]]$group_out_3st$emiss_mean$quant_95,
                                                                    out_p4_3[[sc]]$group_out_3st$emiss_mean$quant_95, 
                                                                    out_p4_4[[sc]]$group_out_3st$emiss_mean$quant_95)
  group_out_3st_emiss_mean[[sc + l_files_p2]]$quant_97.5 <-  Map(rbind,out_p4_1[[sc]]$group_out_3st$emiss_mean$quant_97.5, 
                                                                    out_p4_2[[sc]]$group_out_3st$emiss_mean$quant_97.5,
                                                                    out_p4_3[[sc]]$group_out_3st$emiss_mean$quant_97.5, 
                                                                    out_p4_4[[sc]]$group_out_3st$emiss_mean$quant_97.5)
  group_out_3st_emiss_sd[[sc + l_files_p2]]$median       <-  Map(rbind,out_p4_1[[sc]]$group_out_3st$emiss_sd$median, 
                                                                    out_p4_2[[sc]]$group_out_3st$emiss_sd$median,
                                                                    out_p4_3[[sc]]$group_out_3st$emiss_sd$median, 
                                                                    out_p4_4[[sc]]$group_out_3st$emiss_sd$median)
  group_out_3st_emiss_sd[[sc + l_files_p2]]$quant_2.5    <-  Map(rbind,out_p4_1[[sc]]$group_out_3st$emiss_sd$quant_2.5, 
                                                                    out_p4_2[[sc]]$group_out_3st$emiss_sd$quant_2.5,
                                                                    out_p4_3[[sc]]$group_out_3st$emiss_sd$quant_2.5, 
                                                                    out_p4_4[[sc]]$group_out_3st$emiss_sd$quant_2.5)
  group_out_3st_emiss_sd[[sc + l_files_p2]]$quant_5      <-  Map(rbind,out_p4_1[[sc]]$group_out_3st$emiss_sd$quant_5, 
                                                                    out_p4_2[[sc]]$group_out_3st$emiss_sd$quant_5,
                                                                    out_p4_3[[sc]]$group_out_3st$emiss_sd$quant_5, 
                                                                    out_p4_4[[sc]]$group_out_3st$emiss_sd$quant_5)
  group_out_3st_emiss_sd[[sc + l_files_p2]]$quant_95     <-  Map(rbind,out_p4_1[[sc]]$group_out_3st$emiss_sd$quant_95, 
                                                                    out_p4_2[[sc]]$group_out_3st$emiss_sd$quant_95,
                                                                    out_p4_3[[sc]]$group_out_3st$emiss_sd$quant_95, 
                                                                    out_p4_4[[sc]]$group_out_3st$emiss_sd$quant_95)
  group_out_3st_emiss_sd[[sc + l_files_p2]]$quant_97.5   <-  Map(rbind,out_p4_1[[sc]]$group_out_3st$emiss_sd$quant_97.5, 
                                                                    out_p4_2[[sc]]$group_out_3st$emiss_sd$quant_97.5,
                                                                    out_p4_3[[sc]]$group_out_3st$emiss_sd$quant_97.5, 
                                                                    out_p4_4[[sc]]$group_out_3st$emiss_sd$quant_97.5)
  group_out_3st_gamma_prob[[sc + l_files_p2]]            <-  Map(rbind, out_p4_1[[sc]]$group_out_3st$trans_prob, 
                                                                     out_p4_2[[sc]]$group_out_3st$trans_prob,
                                                                     out_p4_3[[sc]]$group_out_3st$trans_prob, 
                                                                     out_p4_4[[sc]]$group_out_3st$trans_prob)
  group_out_3st_gamma_int[[sc + l_files_p2]]             <-  Map(rbind, out_p4_1[[sc]]$group_out_3st$trans_interc, 
                                                                     out_p4_2[[sc]]$group_out_3st$trans_interc,
                                                                     out_p4_3[[sc]]$group_out_3st$trans_interc, 
                                                                     out_p4_4[[sc]]$group_out_3st$trans_interc)
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

# extracted_results_3st <- readRDS("extracted_results_3st.RDS")
# selection_out <- extracted_results_3st$selection_out

# Model selection ####


AICc_mean <- data.frame(true_m = 3,
                        n = settings4[c(selected, selected_p2, selected_p4), 1], 
                        n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                        n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                        KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                        st1 = NA, 
                        st2 = NA, 
                        st3 = NA, 
                        st4 = NA)
min_AICc <- data.frame(true_m = 3,
                       n = settings4[c(selected, selected_p2, selected_p4), 1], 
                       n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                       n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                       KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                       prop_st1 = NA,
                       prop_st2 = NA,
                       prop_st3 = NA,
                       prop_st4 = NA)

AIC_mean <- data.frame(true_m = 3,
                       n = settings4[c(selected, selected_p2, selected_p4), 1], 
                       n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                       n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                       KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                       st1 = NA, 
                       st2 = NA, 
                       st3 = NA, 
                       st4 = NA)
min_AIC <- data.frame(true_m = 3,
                      n = settings4[c(selected, selected_p2, selected_p4), 1], 
                      n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                      n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                      KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                      prop_st1 = NA,
                      prop_st2 = NA,
                      prop_st3 = NA,
                      prop_st4 = NA)




for(sc in 1:length(c(selected, selected_p2, selected_p4))){
  AICc_mean[sc, 6:9] <- apply(selection_out[[sc]]$AICc_mean, 2, mean)
}

AICc_mean_melt <- melt(AICc_mean, id.vars = colnames(AICc_mean)[1:5], variable.name = "state")

ggplot(data = AICc_mean_melt, mapping = aes(x = state, y = value, group = as.factor(n_t), color = as.factor(n_t))) +
  geom_line(stat = "summary", fun = "mean") + ylab("Mean AICc") +
  facet_grid(rows = vars(as.factor(n)), cols = vars(as.factor(KL_div)))


for(sc in 1:length(c(selected, selected_p2, selected_p4))){
  min_AICc[sc, 6:9] <- table(factor(apply(selection_out[[sc]]$AICc_mean, 1, which.min), levels = 1:4)) / n_sim
}

gg_min_AICc <- melt(min_AICc, id.vars = colnames(min_AICc)[1:5], variable.name = "state")

ggplot(data = gg_min_AICc, mapping = aes(x = state, y = value, fill = as.factor(n_t))) +
  geom_bar(stat = "summary", fun = "mean", color="black", position=position_dodge()) +
  facet_grid(rows = vars(as.factor(n)), cols = vars(as.factor(KL_div))) +
  ylab("Proportion in state using minimum AICc") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

ggplot(data = gg_min_AICc[gg_min_AICc$n_dep == 8,], mapping = aes(x = state, y = value, fill = as.factor(n_t))) +
  geom_bar(stat = "summary", fun = "mean", color="black", position=position_dodge()) +
  facet_grid(rows = vars(as.factor(n)), cols = vars(as.factor(KL_div))) +
  ylab("Proportion in state using minimum AICc") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))


for(sc in 1:length(c(selected, selected_p2, selected_p4))){
  AIC_mean[sc, 6:9] <- apply(selection_out[[sc]]$AIC_mean, 2, mean)
}

for(sc in 1:length(c(selected, selected_p2, selected_p4))){
  min_AIC[sc, 6:9] <- table(factor(apply(selection_out[[sc]]$AIC_mean, 1, which.min), levels = 1:4)) / n_sim
}

gg_min_AIC <- melt(min_AIC, id.vars = colnames(min_AIC)[1:5], variable.name = "state")
gg_min_AIC$state <- factor(gg_min_AIC$state, labels = paste("state", 1:4))

ggplot(data = gg_min_AIC, mapping = aes(x = state, y = value, fill = as.factor(n_t))) +
  geom_bar(stat = "summary", fun = "mean", color="black", position=position_dodge()) +
  ylab("Proportion in state using minimum AIC") + xlab("Proportion in ...") +
  scale_fill_discrete(name = "Number of\n subjects") +
  facet_grid(rows = vars(as.factor(n)), cols = vars(as.factor(KL_div))) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

ggplot(data = gg_min_AIC[gg_min_AIC$n_dep == 8,], mapping = aes(x = state, y = value, fill = as.factor(n_t))) +
  geom_bar(stat = "summary", fun = "mean", color="black", position=position_dodge()) +
  ylab("Proportion in state using minimum AIC") + xlab("Proportion in ...") +
  scale_fill_discrete(name = "Number of\n subjects") +
  facet_grid(rows = vars(as.factor(n)), cols = vars(as.factor(KL_div))) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

model_selection_3st <- data.frame(true_m = 3,
                                  n = settings4[c(selected, selected_p2, selected_p4), 1], 
                                  n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                                  n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                                  KL_div = settings4[c(selected, selected_p2, selected_p4), 4],
                                  mean_AICc_st1 = AICc_mean$st1, 
                                  mean_AICc_st2 = AICc_mean$st2,
                                  mean_AICc_st3 = AICc_mean$st3,
                                  mean_AICc_st4 = AICc_mean$st4,
                                  prop_AICc_st1 = min_AICc$prop_st1, 
                                  prop_AICc_st2 = min_AICc$prop_st2, 
                                  prop_AICc_st3 = min_AICc$prop_st3, 
                                  prop_AICc_st4 = min_AICc$prop_st4, 
                                  mean_AIC_st1 = AIC_mean$st1, 
                                  mean_AIC_st2 = AIC_mean$st2,
                                  mean_AIC_st3 = AIC_mean$st3,
                                  mean_AIC_st4 = AIC_mean$st4,
                                  prop_AIC_st1 = min_AIC$prop_st1, 
                                  prop_AIC_st2 = min_AIC$prop_st2, 
                                  prop_AIC_st3 = min_AIC$prop_st3, 
                                  prop_AIC_st4 = min_AIC$prop_st4)

# Model performance ####

## Bias ####
### Emission distribution - MEAN ####
# group_out_3st_emiss_mean <- extracted_results_3st$group_out_3st_emiss_mean
# means_true <- extracted_results_3st$means_true

mean_abs_bias_emiss <- data.frame(true_m = 3,
                                  n = settings4[c(selected, selected_p2, selected_p4), 1], 
                                  n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                                  n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                                  KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                                  abs_bias = NA,
                                  abs_bias2 = NA)

mean_rel_bias_emiss <- data.frame(true_m = 3,
                                  n = settings4[c(selected, selected_p2, selected_p4), 1], 
                                  n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                                  n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                                  KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                                  rel_bias = NA,
                                  rel_bias2 = NA)

for(sc in 1:length(c(selected, selected_p2, selected_p4))){
  diff <- mapply('-', group_out_3st_emiss_mean[[sc]]$median, means_true[[sc]], SIMPLIFY = FALSE)
  abs.diff <- sapply(diff, 'abs')
  mean_abs_bias_emiss[sc, 6] <- median(abs.diff)
  mean_abs_bias_emiss[sc, 7] <- median(sapply(diff, mean))
  rel.diff <- mapply('/', diff, means_true[[sc]])
  mean_rel_bias_emiss[sc, 6] <- median(rel.diff)
  mean_rel_bias_emiss[sc, 7] <-  median(abs(rel.diff))
}

aggregate(mean_rel_bias_emiss, by = list(mean_rel_bias_emiss$KL_div), FUN = median)


ggplot(data = mean_abs_bias_emiss, mapping = aes(x = n_t, y = abs_bias, group = as.factor(n), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  ylab("mean absolute bias") +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep)))


ggplot(data = mean_rel_bias_emiss, mapping = aes(x = n_t, y = rel_bias, group = as.factor(n), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  ylab("median relative bias") +
  scale_color_discrete(name = "Number of\nsubjects") +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep))) + ylim(-0.5, 0.5) +
  geom_hline(yintercept=c(-0.1, 0.1), linetype="dashed", color = "grey") + xlab("number of observations per subject")


ggplot(data = mean_rel_bias_emiss, mapping = aes(x = n_t, y = rel_bias2, group = as.factor(n), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  ylab("median absolute relative bias") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep))) + ylim(0,1) +
  geom_hline(yintercept=0.1, linetype="dashed", color = "grey") + xlab("number of observations per subject")

# 9 panels 
ggplot(data = mean_rel_bias_emiss, mapping = aes(x = n_t, y = rel_bias2, group = interaction(as.factor(n), as.factor(n_dep)), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  geom_line() +
  ylab("median absolute relative bias") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep))) +
  geom_hline(yintercept=0.1, linetype="dashed", color = "grey") + xlab("number of observations per subject") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

# 3 panels 
ggplot(data = mean_rel_bias_emiss, mapping = aes(x = n_t, y = rel_bias2, group = interaction(as.factor(n), as.factor(n_dep)), color = as.factor(n_dep))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  geom_line(aes(linetype = as.factor(n))) +
  ylab("median absolute relative bias") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(cols = vars(as.factor(KL_div))) + ylim(0,1) +
  geom_hline(yintercept=0.1, linetype="dashed", color = "grey") + xlab("number of observations per subject")

### Emission distribution - SD ####
bias_emiss_SD <- data.frame(true_m = 3,
                                  n = settings4[c(selected, selected_p2, selected_p4), 1], 
                                  n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                                  n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                                  KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                                  mean_bias = NA,
                                  med_abs_rel_bias = NA)

for(sc in 1:length(c(selected, selected_p2, selected_p4))){
  diff <- mapply('-', group_out_3st_emiss_sd[[sc]]$median, 1, SIMPLIFY = FALSE)
  bias_emiss_SD[sc, 6] <- mean(unlist(diff))
  bias_emiss_SD[sc, 7] <-  median(mapply('/', diff, 1))
}

aggregate(bias_emiss_SD, by = list(bias_emiss_SD$KL_div), FUN = median)

ggplot(data = bias_emiss_SD, mapping = aes(x = n_t, y = med_abs_rel_bias, group = interaction(as.factor(n), as.factor(n_dep)), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  geom_line() +
  ylab("median absolute relative bias") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep))) + ylim(0,0.75) +
  geom_hline(yintercept=0.1, linetype="dashed", color = "grey") + xlab("number of observations per subject") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))


### Explaining bias in emission ####
#### Quantifying likely label switching problem occurences ####

rmse <- function(x){
  row_mean <- mean(x)
  rmse <- sqrt(sum((mean(x) - x)^2 / true_m))
  return(rmse)
}

proxy_labelsw <- vector("list", length = length(c(selected, selected_p2, selected_p4)))

for(sc in 1:length(c(selected, selected_p2, selected_p4))){
  proxy_labelsw[[sc]] <- lapply(group_out_3st_emiss_mean[[sc]]$median, apply, 1, rmse)
}

label_switch <- data.frame(true_m = 3,
                                  n = settings4[c(selected, selected_p2, selected_p4), 1], 
                                  n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                                  n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                                  KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                                  label_sw_prop = NA)


for(sc in 1:length(c(selected, selected_p2, selected_p4))){
  per_var <-  mapply("<", proxy_labelsw[[sc]], 0.20) * 1
  label_switch[sc, 6] <- sum(apply(per_var, 1, prod)) / n_sim
}


ggplot(data = label_switch, mapping = aes(x = n_t, y = label_sw_prop, group = interaction(as.factor(n), as.factor(n_dep)), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  geom_line() +
  ylab("Proportion of suspected label switching") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep))) + 
  xlab("number of observations per subject") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))


# give some examples of rmse distribution in appendix? 
# .2 seems to be good cutt off for ndep 2 an 4, 
# 65 >- 2 dep, KKLD of 5, 30 auvjects, 200. measurements
# 128 -> 4 dep, KLD of 5, 120 subjects, 400 measurments
# 61 -> 8 dep, KLD of 3, 30 subjects, 100 measurements
hist(unlist(proxy_labelsw[[65]]), breaks = 20)
group_out_3st_emiss_mean[[65]]$median

hist(unlist(proxy_labelsw[[128]]), breaks = 20)
group_out_3st_emiss_mean[[128]]$median

hist(unlist(proxy_labelsw[[61]]), breaks = 20)
group_out_3st_emiss_mean[[61]]$median

# can also calculate RMSE for true means

### Gamma, prob scale ####


true_gamma <-  matrix(c(0.7, 0.2, 0.1, 
                        0.1, 0.8, 0.1, 
                        0.1, 0.1, 0.8), ncol = true_m, byrow = TRUE)

diag_gamma <-  diag(true_gamma)
off_diag_gamma <- t(true_gamma)[lower.tri(true_gamma) | upper.tri(true_gamma)]
mean_bias_gamma <- data.frame(true_m = 3,
                                  n = settings4[c(selected, selected_p2, selected_p4), 1], 
                                  n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                                  n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                                  KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                                  bias_diag = NA,
                                  bias_off_diag = NA)

mean_rel_bias_gamma <- data.frame(true_m = 3,
                                  n = settings4[c(selected, selected_p2, selected_p4), 1], 
                                  n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                                  n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                                  KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                                  bias_diag = NA,
                                  bias_off_diag = NA)

for(sc in 1:length(c(selected, selected_p2, selected_p4))){
  mean_bias_gamma[sc, 6] <- mean(apply(group_out_3st_gamma_prob[[sc]]$median[,c(1, 5, 9)], 2, median) - diag_gamma)
  mean_bias_gamma[sc, 7] <- mean(apply(group_out_3st_gamma_prob[[sc]]$median[,-c(1, 5, 9)], 2, median) - off_diag_gamma)
  mean_rel_bias_gamma[sc, 6] <-  median((apply(group_out_3st_gamma_prob[[sc]]$median[,c(1, 5, 9)], 2, median) - diag_gamma) / diag_gamma)
  mean_rel_bias_gamma[sc, 7] <-  median((apply(group_out_3st_gamma_prob[[sc]]$median[,-c(1, 5, 9)], 2, median) - off_diag_gamma) / off_diag_gamma)
}

aggregate(mean_rel_bias_gamma, by = list(mean_rel_bias_gamma$KL_div), FUN = mean)

gg_mean_bias_gamma <- melt(mean_bias_gamma, id.vars = colnames(mean_bias_gamma)[1:5], variable.name = "diag")

ggplot(data = gg_mean_bias_gamma, mapping = aes(x = n_t, y = value, group = interaction(as.factor(n), as.factor(diag)), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  geom_line(aes(linetype = as.factor(diag))) +
  ylab("mean bias") +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep))) +
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_linetype_manual(name = "Transition type", 
                        values = c(1, 3),
                        labels = c("Self-transition", "Between state transition")) +
  # scale_linetype(name = "Transition type", values=c("Self-transition" = 'bias_dag',"Between state transition" = 'bias_off_diag')) +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  geom_hline(yintercept=0, linetype="solid", color = "grey") + xlab("number of observations per subject") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

ggplot(data = gg_mean_bias_gamma, mapping = aes(x = n_t, y = value, group = interaction(as.factor(n), as.factor(n_dep)), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  geom_line(aes(linetype = as.factor(n_dep))) +
  ylab("mean bias") +
  facet_grid(rows =vars(as.factor(diag)) , cols = vars(as.factor(KL_div))) +
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  geom_hline(yintercept=0, linetype="solid", color = "grey") + xlab("number of observations per subject")


gg_mean_rel_bias_gamma <- melt(mean_rel_bias_gamma, id.vars = colnames(mean_rel_bias_gamma)[1:5], variable.name = "diag")

ggplot(data = gg_mean_rel_bias_gamma, mapping = aes(x = n_t, y = value, group = as.factor(n), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  ylab("median relative bias") +
  geom_line(aes(linetype = as.factor(diag))) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep))) +
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_linetype_manual(name = "Transition type", 
                        values = c(1, 3),
                        labels = c("Self-transition", "Between state transition")) +
  # scale_linetype(name = "Transition type", values=c("Self-transition" = 'bias_dag',"Between state transition" = 'bias_off_diag')) +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  geom_hline(yintercept=0, linetype="solid", color = "grey") + xlab("number of observations per subject") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

ggplot(data = gg_mean_rel_bias_gamma, mapping = aes(x = n_t, y = value, group = interaction(as.factor(n), as.factor(n_dep)), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  geom_line(aes(linetype = as.factor(n_dep))) +
  ylab("median relative bias") +
  facet_grid(rows =vars(as.factor(diag)) , cols = vars(as.factor(KL_div))) +
  geom_hline(yintercept=c(-0.1, 0.1), linetype=c("dashed"), color = "grey") + 
  geom_hline(yintercept=c(0), linetype=c("solid"), color = "grey") + 
  xlab("number of observations per subject")



### Gamma, mnl intercept scale ####


## Empirical standard error / precision (sqrt(var(theta))) ####
### Emission distribution - MEAN ####
### NOT POSSIBLE, NOT THE SAME VALUE


### Emission distribution - SD ####
# group_out_3st_emiss_sd <- extracted_results_3st$group_out_3st_emiss_sd

SD_emp_se <- data.frame(true_m = 3,
                           n = settings4[c(selected, selected_p2, selected_p4), 1], 
                           n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                           n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                           KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                           emp_se = NA)

for(sc in 1:length(c(selected, selected_p2, selected_p4))){
  SD_emp_se[sc, 6] <- mean(sqrt(unlist(lapply(group_out_3st_emiss_sd[[sc]]$median, apply, 2, var))))
}

max(SD_emp_se$emp_se[SD_emp_se$KL_div == 7 & SD_emp_se$n_dep == 2])
max(SD_emp_se$emp_se[SD_emp_se$KL_div == 7 & SD_emp_se$n_dep == 4])
max(SD_emp_se$emp_se[SD_emp_se$KL_div == 7 & SD_emp_se$n_dep == 8])

ggplot(data = SD_emp_se, mapping = aes(x = n_t, y = emp_se, group = interaction(as.factor(n), as.factor(n_dep)), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  geom_line() +
  ylab("Emperical SE of emission SD") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep))) + 
  xlab("number of observations per subject") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))


### Gamma, prob scale ####
emp_se_gamma <- data.frame(true_m = 3,
                              n = settings4[c(selected, selected_p2, selected_p4), 1], 
                              n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                              n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                              KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                              emp_SE_diag = NA,
                              emp_SE_off_diag = NA)

for(sc in 1:length(c(selected, selected_p2, selected_p4))){
  emp_se_gamma[sc, 6] <- mean(sqrt(apply(group_out_3st_gamma_prob[[sc]]$median[,c(1, 5, 9)], 2, var)))
  emp_se_gamma[sc, 7] <-  mean(sqrt(apply(group_out_3st_gamma_prob[[sc]]$median[,-c(1, 5, 9)], 2, var)))
}

gg_emp_se_gamma <- melt(emp_se_gamma, id.vars = colnames(emp_se_gamma)[1:5], variable.name = "diag")

ggplot(data = gg_emp_se_gamma, mapping = aes(x = n_t, y = value, group = interaction(as.factor(n), as.factor(diag)), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  geom_line(aes(linetype = as.factor(diag))) +
  ylab("Emperical SE") +
  facet_grid(rows = vars(as.factor(n_dep)), cols = vars(as.factor(KL_div))) +
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_linetype_manual(name = "Transition type", 
                        values = c(1, 3),
                        labels = c("Self-transition", "Between state transition")) +
  # scale_linetype(name = "Transition type", values=c("Self-transition" = 'bias_dag',"Between state transition" = 'bias_off_diag')) +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  geom_hline(yintercept=0, linetype="solid", color = "grey") + xlab("number of observations per subject") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))


aggregate(emp_se_gamma, by = list(emp_se_gamma$n_t), FUN = mean)
aggregate(emp_se_gamma, by = list(emp_se_gamma$n), FUN = mean)

## Coverage ####
### Emission distribution - MEAN ####
Coverage_emiss <- data.frame(true_m = 3,
                             n = settings4[c(selected, selected_p2, selected_p4), 1], 
                             n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                             n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                             KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                             cov_90 = NA, 
                             cov_95 = NA)

for(sc in 1:length(c(selected, selected_p2, selected_p4))){
  lower_90 <- lapply(mapply('<', group_out_3st_emiss_mean[[sc]]$quant_5, means_true[[sc]], SIMPLIFY = FALSE), "*", 1)
  upper_90 <- lapply(mapply('>', group_out_3st_emiss_mean[[sc]]$quant_95, means_true[[sc]], SIMPLIFY = FALSE), "*", 1)
  Coverage_emiss[sc, 6] <- sum(mapply("*", lower_90, upper_90)) / (n_sim * Coverage_emiss$n_dep[sc] * true_m) * 100
  lower_95 <- lapply(mapply('<', group_out_3st_emiss_mean[[sc]]$quant_2.5, means_true[[sc]], SIMPLIFY = FALSE), "*", 1)
  upper_95 <- lapply(mapply('>', group_out_3st_emiss_mean[[sc]]$quant_97.5, means_true[[sc]], SIMPLIFY = FALSE), "*", 1)
  Coverage_emiss[sc, 7] <- sum(mapply("*", lower_95, upper_95)) / (n_sim * Coverage_emiss$n_dep[sc]  * true_m) * 100
}

ggplot(data = Coverage_emiss, mapping = aes(x = n_t, y = cov_95, group = as.factor(n), color = as.factor(n))) +
  geom_point() + 
  geom_line() +
  ylab("Coverage emission distribution") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep))) +
  xlab("number of observations per subject") +
  geom_hline(yintercept=c(90, 95), linetype="dashed", color = "grey") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

### Emission distribution - SD ####
Coverage_emiss_SD <- data.frame(true_m = 3,
                             n = settings4[c(selected, selected_p2, selected_p4), 1], 
                             n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                             n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                             KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                             cov_90 = NA, 
                             cov_95 = NA, 
                             cov_95_width = NA)

for(sc in 1:length(c(selected, selected_p2, selected_p4))){
  lower_90 <- lapply(mapply('<', group_out_3st_emiss_sd[[sc]]$quant_5, 1, SIMPLIFY = FALSE), "*", 1)
  upper_90 <- lapply(mapply('>', group_out_3st_emiss_sd[[sc]]$quant_95,1, SIMPLIFY = FALSE), "*", 1)
  Coverage_emiss_SD[sc, 6] <- sum(mapply("*", lower_90, upper_90)) / (n_sim * Coverage_emiss_SD$n_dep[sc] * true_m) * 100
  lower_95 <- lapply(mapply('<', group_out_3st_emiss_sd[[sc]]$quant_2.5, 1, SIMPLIFY = FALSE), "*", 1)
  upper_95 <- lapply(mapply('>', group_out_3st_emiss_sd[[sc]]$quant_97.5, 1, SIMPLIFY = FALSE), "*", 1)
  Coverage_emiss_SD[sc, 7] <- sum(mapply("*", lower_95, upper_95)) / (n_sim * Coverage_emiss_SD$n_dep[sc]  * true_m) * 100
  Coverage_emiss_SD[sc, 8] <- mean(unlist(group_out_3st_emiss_sd[[sc]]$quant_97.5)) - mean(unlist(group_out_3st_emiss_sd[[sc]]$quant_2.5))
}

aggregate(Coverage_emiss_SD, by = list(Coverage_emiss_SD$n), FUN = mean)
aggregate(Coverage_emiss_SD, by = list(Coverage_emiss_SD$n_t), FUN = mean)


ggplot(data = Coverage_emiss_SD, mapping = aes(x = n_t, y = cov_95, group = as.factor(n), color = as.factor(n))) +
  geom_point() + 
  geom_line() +
  ylab("Coverage emission distribution") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep))) +
  xlab("number of observations per subject") +
  geom_hline(yintercept=c(90, 95), linetype="dashed", color = "grey") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

## Emission performance table ####
emission_performance_3st <- data.frame(true_m = 3,
                                      n = settings4[c(selected, selected_p2, selected_p4), 1], 
                                      n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                                      n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                                      KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                                      MEAN_med_abs_rel_bias = mean_rel_bias_emiss$rel_bias2,
                                      label_switch_prop = label_switch$label_sw_prop,
                                      MEAN_cov_CrI_95 = Coverage_emiss$cov_95, 
                                      SD_med_abs_rel_bias = bias_emiss_SD$med_abs_rel_bias,
                                      SD_mean_bias = bias_emiss_SD$mean_bias, 
                                      SD_emperical_sd = SD_emp_se$emp_se, 
                                      SD_cov_CrI_95 = Coverage_emiss_SD$cov_95)




### Gamma, prob scale ####
true_gamma <-  matrix(c(0.7, 0.2, 0.1, 
                        0.1, 0.8, 0.1, 
                        0.1, 0.1, 0.8), ncol = true_m, byrow = TRUE)

diag_gamma <-  diag(true_gamma)
off_diag_gamma <- t(true_gamma)[lower.tri(true_gamma) | upper.tri(true_gamma)]

Coverage_gamma <- data.frame(true_m = 3,
                             n = settings4[c(selected, selected_p2, selected_p4), 1], 
                             n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                             n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                             KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                             cov_90_diag = NA, 
                             cov_90_off = NA,
                             cov_95_diag = NA,
                             cov_95_off = NA)



for(sc in 1:length(c(selected, selected_p2, selected_p4))){
  lower_90_diag <- (group_out_3st_gamma_prob[[sc]]$quant_5[,c(1, 5, 9)] < matrix(diag_gamma, nrow = n_sim, ncol = true_m, byrow = TRUE)) * 1
  upper_90_diag <-  (group_out_3st_gamma_prob[[sc]]$quant_95[,c(1, 5, 9)] > matrix(diag_gamma, nrow = n_sim, ncol = true_m, byrow = TRUE)) * 1
  lower_90_off <-  (group_out_3st_gamma_prob[[sc]]$quant_5[,-c(1, 5, 9)] < matrix(off_diag_gamma, nrow = n_sim, ncol = true_m*(true_m-1), byrow = TRUE)) * 1
  upper_90_off <-  (group_out_3st_gamma_prob[[sc]]$quant_95[,-c(1, 5, 9)] > matrix(off_diag_gamma, nrow = n_sim, ncol = true_m*(true_m-1), byrow = TRUE)) * 1
  
  Coverage_gamma[sc, 6] <- sum(lower_90_diag *upper_90_diag) / (n_sim * true_m) * 100
  Coverage_gamma[sc, 7] <- sum(lower_90_off *upper_90_off) / (n_sim * true_m*(true_m-1)) * 100
  
  lower_95_diag <-  (group_out_3st_gamma_prob[[sc]]$quant_2.5[,c(1, 5, 9)] < matrix(diag_gamma, nrow = n_sim, ncol = true_m, byrow = TRUE)) * 1
  upper_95_diag <- (group_out_3st_gamma_prob[[sc]]$quant_95[,c(1, 5, 9)] > matrix(diag_gamma, nrow = n_sim, ncol = true_m, byrow = TRUE)) * 1
  lower_95_off <-  (group_out_3st_gamma_prob[[sc]]$quant_2.5[,-c(1, 5, 9)] < matrix(off_diag_gamma, nrow = n_sim, ncol = true_m*(true_m-1), byrow = TRUE)) * 1
  upper_95_off <- (group_out_3st_gamma_prob[[sc]]$quant_95[,-c(1, 5, 9)] > matrix(off_diag_gamma, nrow = n_sim, ncol = true_m*(true_m-1), byrow = TRUE)) * 1
  Coverage_gamma[sc, 8] <- sum(lower_95_diag *upper_95_diag) / (n_sim * true_m) * 100
  Coverage_gamma[sc, 9] <- sum(lower_95_off *upper_95_off) / (n_sim * true_m*(true_m-1)) * 100
}



ggplot(data = Coverage_gamma, mapping = aes(x = n_t, y = cov_95_diag, group = as.factor(n), color = as.factor(n))) +
  geom_point() + 
  geom_line() +
  ylab("Coverage gamma (diagonal)") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep))) +
  xlab("number of observations per subject") +
  geom_hline(yintercept=c(90, 95), linetype="dashed", color = "grey") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

### Gamma, mnl intercept scale ####

## GAMMA performance table ####

gamma_performance_3st <- data.frame(true_m = 3,
                                    n = settings4[c(selected, selected_p2, selected_p4), 1], 
                                    n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                                    n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                                    KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                                    mean_bias_diag = mean_bias_gamma$bias_diag,
                                    mean_bias_off = mean_bias_gamma$bias_off_diag,
                                    emperical_sd_diag = emp_se_gamma$emp_SE_diag, 
                                    emperical_sd_off = emp_se_gamma$emp_SE_off_diag, 
                                    cov_CrI_95_diag = Coverage_gamma$cov_95_diag,
                                    cov_CrI_95_off = Coverage_gamma$cov_95_off)

## State Decoding ####
# inferred_states  <- extracted_results_3st$inferred_states
# true_states  <- extracted_results_3st$true_states
# true_state_probs  <- extracted_results_3st$true_state_probs

state_decoding <- data.frame(true_m = 3,
                           n = settings4[c(selected, selected_p2, selected_p4), 1], 
                           n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                           n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                           KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                           mean_prop_correct = NA,
                           mean_kappa = NA, 
                           mean_prop_over_20 = NA)

### Mean proportion of states classified correctly ####

for(sc in 1:length(c(selected, selected_p2, selected_p4))){
  state_decoding[sc, 6] <- sum(inferred_states[[sc]][,-1] == true_states[[sc]][,-1]) / (dim(true_states[[sc]][,-1])[1] * dim(true_states[[sc]][,-1])[2])
  inferred_totals <- table(inferred_states[[sc]][,-1]) / 100
  true_totals <- table(true_states[[sc]][,-1]) / 100
  matrix_true_inferred <- table(inferred_states[[sc]][,-1], true_states[[sc]][,-1]) / 100
  agreements <- sum(diag(matrix_true_inferred))
  change_agreements <- sum(inferred_totals * true_totals / sum(true_totals))
  state_decoding[sc, 7] <- (agreements - change_agreements) / ( sum(true_totals) - change_agreements)
  state_decoding[sc, 8] <- sum(true_state_probs[[sc]][,-1] > 0.20 )/ (dim(true_states[[sc]][,-1])[1] * dim(true_states[[sc]][,-1])[2])
}

mean(state_decoding$mean_prop_correct[state_decoding$KL_div > 3])
min(state_decoding$mean_prop_correct[state_decoding$KL_div > 3])
max(state_decoding$mean_prop_correct[state_decoding$KL_div > 3])

sum(state_decoding$mean_prop_correct[state_decoding$KL_div > 3] >= 0.8) / length(state_decoding$mean_prop_correct[state_decoding$KL_div > 3] >= 0.8)

mean(state_decoding$mean_kappa[state_decoding$KL_div > 5])
min(state_decoding$mean_kappa[state_decoding$KL_div > 5])
max(state_decoding$mean_kappa[state_decoding$KL_div > 5])
sum(state_decoding$mean_kappa[state_decoding$KL_div > 5] >= 0.8) / length(state_decoding$mean_kappa[state_decoding$KL_div > 5] >= 0.8)


ggplot(data = state_decoding, mapping = aes(x = n_t, y = mean_prop_correct, group = interaction(as.factor(n), as.factor(n_dep)), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  geom_line() +
  ylab("mean proportion correctly classified state") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep))) + 
  xlab("number of observations per subject") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))

ggplot(data = state_decoding, mapping = aes(x = n_t, y = mean_kappa, group = interaction(as.factor(n), as.factor(n_dep)), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_point() + 
  geom_line() +
  ylab("Cohen's kappa") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep))) + 
  xlab("number of observations per subject") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))


# Result tables ####

Result_table_3st <- list(model_selection_3st = model_selection_3st,
                         emission_performance_3st = emission_performance_3st,
                         gamma_performance_3st = gamma_performance_3st,
                         state_decoding_3st = state_decoding)

saveRDS(Result_table_3st, file = "result_tables/Result_table_3st.RDS")
# Convergence ####
