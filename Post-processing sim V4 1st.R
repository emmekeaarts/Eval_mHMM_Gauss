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
# 2 states
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

# Model selection ####


AICc_mean <- data.frame(true_m = 2,
                        n = settings4[c(selected, selected_p2, selected_p4), 1], 
                        n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                        n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                        KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                        st1 = NA, 
                        st2 = NA, 
                        st3 = NA, 
                        st4 = NA)
min_AICc <- data.frame(true_m = 2,
                       n = settings4[c(selected, selected_p2, selected_p4), 1], 
                       n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                       n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                       KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                       prop_st1 = NA,
                       prop_st2 = NA,
                       prop_st3 = NA,
                       prop_st4 = NA)

AIC_mean <- data.frame(true_m = 2,
                       n = settings4[c(selected, selected_p2, selected_p4), 1], 
                       n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                       n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                       KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                       st1 = NA, 
                       st2 = NA, 
                       st3 = NA, 
                       st4 = NA)
min_AIC <- data.frame(true_m = 2,
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

model_selection_1st <- data.frame(true_m = 2,
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


## MOdel performance ####
## Bias ####
### Emission distribution - MEAN ####
mean_abs_bias_emiss <- data.frame(true_m = 1,
                                  n = settings4[c(selected, selected_p2, selected_p4), 1], 
                                  n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                                  n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                                  KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                                  abs_bias = NA,
                                  abs_bias2 = NA)

mean_rel_bias_emiss <- data.frame(true_m = 1,
                                  n = settings4[c(selected, selected_p2, selected_p4), 1], 
                                  n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                                  n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                                  KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                                  rel_bias = NA,
                                  rel_bias2 = NA)

for(sc in 1:length(c(selected, selected_p2, selected_p4))){
  diff <- mapply('-', group_out_1st_emiss_mean[[sc]]$median, means_true[[sc]], SIMPLIFY = FALSE)
  abs.diff <- sapply(diff, 'abs')
  mean_abs_bias_emiss[sc, 6] <- median(abs.diff)
  mean_abs_bias_emiss[sc, 7] <- median(sapply(diff, mean))
  rel.diff <- mapply('/', diff, means_true[[sc]])
  mean_rel_bias_emiss[sc, 6] <- median(rel.diff)
  mean_rel_bias_emiss[sc, 7] <-  median(abs(rel.diff))
}

### Emission distribution - SD ####
bias_emiss_SD <- data.frame(true_m = 1,
                            n = settings4[c(selected, selected_p2, selected_p4), 1], 
                            n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                            n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                            KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                            mean_bias = NA,
                            med_abs_rel_bias = NA)

for(sc in 1:length(c(selected, selected_p2, selected_p4))){
  diff <- mapply('-', group_out_1st_emiss_sd[[sc]]$median, 1, SIMPLIFY = FALSE)
  bias_emiss_SD[sc, 6] <- mean(unlist(diff))
  bias_emiss_SD[sc, 7] <-  median(mapply('/', diff, 1))
}

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

## Empirical standard error / precision (sqrt(var(theta))) ####
### Emission distribution - MEAN ####
### NOT POSSIBLE, NOT THE SAME VALUE


### Emission distribution - SD ####
SD_emp_se <- data.frame(true_m = 1,
                        n = settings4[c(selected, selected_p2, selected_p4), 1], 
                        n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                        n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                        KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                        emp_se = NA)

for(sc in 1:length(c(selected, selected_p2, selected_p4))){
  SD_emp_se[sc, 6] <- mean(sqrt(unlist(lapply(group_out_1st_emiss_sd[[sc]]$median, apply, 2, var))))
}


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


## Coverage ####
### Emission distribution - MEAN ####
Coverage_emiss <- data.frame(true_m = 1,
                             n = settings4[c(selected, selected_p2, selected_p4), 1], 
                             n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                             n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                             KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                             cov_90 = NA, 
                             cov_95 = NA)

for(sc in 1:length(c(selected, selected_p2, selected_p4))){
  lower_90 <- lapply(mapply('<', group_out_1st_emiss_mean[[sc]]$quant_5, means_true[[sc]], SIMPLIFY = FALSE), "*", 1)
  upper_90 <- lapply(mapply('>', group_out_1st_emiss_mean[[sc]]$quant_95, means_true[[sc]], SIMPLIFY = FALSE), "*", 1)
  Coverage_emiss[sc, 6] <- sum(mapply("*", lower_90, upper_90)) / (n_sim * Coverage_emiss$n_dep[sc] * true_m) * 100
  lower_95 <- lapply(mapply('<', group_out_1st_emiss_mean[[sc]]$quant_2.5, means_true[[sc]], SIMPLIFY = FALSE), "*", 1)
  upper_95 <- lapply(mapply('>', group_out_1st_emiss_mean[[sc]]$quant_97.5, means_true[[sc]], SIMPLIFY = FALSE), "*", 1)
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
Coverage_emiss_SD <- data.frame(true_m = 1,
                                n = settings4[c(selected, selected_p2, selected_p4), 1], 
                                n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                                n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                                KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                                cov_90 = NA, 
                                cov_95 = NA)

for(sc in 1:length(c(selected, selected_p2, selected_p4))){
  lower_90 <- lapply(mapply('<', group_out_1st_emiss_sd[[sc]]$quant_5, 1, SIMPLIFY = FALSE), "*", 1)
  upper_90 <- lapply(mapply('>', group_out_1st_emiss_sd[[sc]]$quant_95,1, SIMPLIFY = FALSE), "*", 1)
  Coverage_emiss_SD[sc, 6] <- sum(mapply("*", lower_90, upper_90)) / (n_sim * Coverage_emiss$n_dep[sc] * true_m) * 100
  lower_95 <- lapply(mapply('<', group_out_1st_emiss_sd[[sc]]$quant_2.5, 1, SIMPLIFY = FALSE), "*", 1)
  upper_95 <- lapply(mapply('>', group_out_1st_emiss_sd[[sc]]$quant_97.5, 1, SIMPLIFY = FALSE), "*", 1)
  Coverage_emiss_SD[sc, 7] <- sum(mapply("*", lower_95, upper_95)) / (n_sim * Coverage_emiss$n_dep[sc]  * true_m) * 100
}

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
emission_performance_1st <- data.frame(true_m = 1,
                                       n = settings4[c(selected, selected_p2, selected_p4), 1], 
                                       n_t = settings4[c(selected, selected_p2, selected_p4), 2],
                                       n_dep = settings4[c(selected, selected_p2, selected_p4), 3], 
                                       KL_div = settings4[c(selected, selected_p2, selected_p4), 4], 
                                       MEAN_med_abs_rel_bias = mean_rel_bias_emiss$rel_bias2,
                                       MEAN_cov_CrI_95 = Coverage_emiss$cov_95, 
                                       SD_med_abs_rel_bias = bias_emiss_SD$med_abs_rel_bias,
                                       SD_mean_bias = bias_emiss_SD$mean_bias, 
                                       SD_emperical_sd = SD_emp_se$emp_se, 
                                       SD_cov_CrI_95 = Coverage_emiss_SD$cov_95)




Result_table_1st <- list(model_selection_1st = model_selection_1st,
                         emission_performance_1st = emission_performance_1st)

saveRDS(Result_table_1st, file = "result_tables/Result_table_1st.RDS")
