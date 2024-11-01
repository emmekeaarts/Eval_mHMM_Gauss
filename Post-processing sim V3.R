library(mHMMbayes)
library(ggplot2)
library(reshape2)

# Simulstion settings
settings3 <- read.csv("sim_scenarios_V3.csv")
true_m  <- 3

# scenario <- 1
# KL_div  <- settings[scenario, 4]
# n_t     <- settings[scenario, 2]
# n       <- settings[scenario, 1]
# n_dep   <- settings[scenario, 3]
# out_file <- paste0("out_n", n, "_nt", n_t, "_ndep", n_dep, "_", true_m, "st_KLD", KL_div)

## Reading in output files ####
# 3 states
selected <- sort(c(1:90, 101:106, 108))[-c(12, 38, 41, 44, 45, 48, 57, 65, 67, 73, 75, 80, 81)]


files <- paste0("out_n", settings3[selected, 1], "_nt", settings3[selected, 2], "_ndep", settings3[selected, 3], "_", true_m, "st_KLD", settings3[selected, 4])

# Map(rbind, L1, L2)
files_p2_1 <- 
files_p2_2 <- 
# out <- readRDS(paste0("out/", out_file, ".rds")) 

out <- vector("list", length(files))
names(out) <- files

for(sc in 1:length(files)){
  out[[sc]] <- readRDS(paste0("out_V3/", files[sc], ".rds")) 
}

n_sim <- out[[1]]$input$n_sim
J <- out[[1]]$input$J
burn_in <- out[[1]]$input$burn_in

# Model selection ####

out[[80]]$selection_out

AICc_mean <- data.frame(true_m = 3,
                        n = settings3[selected, 1], 
                        n_t = settings3[selected, 2],
                        n_dep = settings3[selected, 3], 
                        KL_div = settings3[selected, 4], 
                        st1 = NA, 
                        st2 = NA, 
                        st3 = NA, 
                        st4 = NA)
min_AICc <- data.frame(true_m = 3,
                       n = settings3[selected, 1], 
                       n_t = settings3[selected, 2],
                       n_dep = settings3[selected, 3], 
                       KL_div = settings3[selected, 4], 
                       prop_st1 = NA,
                       prop_st2 = NA,
                       prop_st3 = NA,
                       prop_st4 = NA)

AIC_mean <- data.frame(true_m = 3,
                       n = settings3[selected, 1], 
                       n_t = settings3[selected, 2],
                       n_dep = settings3[selected, 3], 
                       KL_div = settings3[selected, 4], 
                       st1 = NA, 
                       st2 = NA, 
                       st3 = NA, 
                       st4 = NA)
min_AIC <- data.frame(true_m = 3,
                      n = settings3[selected, 1], 
                      n_t = settings3[selected, 2],
                      n_dep = settings3[selected, 3], 
                      KL_div = settings3[selected, 4], 
                      prop_st1 = NA,
                      prop_st2 = NA,
                      prop_st3 = NA,
                      prop_st4 = NA)

for(sc in 1:length(selected)){
  AICc_mean[sc, 6:9] <- apply(out[[sc]]$selection_out$AICc_mean, 2, mean)
}

AICc_mean <- melt(AICc_mean, id.vars = colnames(AICc_mean)[1:5], variable.name = "state")

ggplot(data = AICc_mean, mapping = aes(x = state, y = value, group = as.factor(n_t), color = as.factor(n_t))) +
  geom_line(stat = "summary", fun = "mean") + ylab("Mean AICc") +
  facet_grid(rows = vars(as.factor(n)), cols = vars(as.factor(KL_div)))


for(sc in 1:length(selected)){
  min_AICc[sc, 6:9] <- table(factor(apply(out[[sc]]$selection_out$AICc_mean, 1, which.min), levels = 1:4)) / n_sim
}

gg_min_AICc <- melt(min_AICc, id.vars = colnames(min_AICc)[1:5], variable.name = "state")

ggplot(data = gg_min_AICc, mapping = aes(x = state, y = value, fill = as.factor(n_t))) +
  geom_bar(stat = "summary", fun = "mean", color="black", position=position_dodge()) +
  facet_grid(rows = vars(as.factor(n)), cols = vars(as.factor(KL_div))) +
  ylab("Proportion in state using minimum AICc") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))


for(sc in 1:length(selected)){
  AIC_mean[sc, 6:9] <- apply(out[[sc]]$selection_out$AIC_mean, 2, mean)
}

for(sc in 1:length(selected)){
  min_AIC[sc, 6:9] <- table(factor(apply(out[[sc]]$selection_out$AIC_mean, 1, which.min), levels = 1:4)) / n_sim
}

gg_min_AIC <- melt(min_AIC, id.vars = colnames(min_AIC)[1:5], variable.name = "state")

ggplot(data = gg_min_AIC, mapping = aes(x = state, y = value, fill = as.factor(n_t))) +
  geom_bar(stat = "summary", fun = "mean", color="black", position=position_dodge()) +
  ylab("Proportion in state using minimum AIC") +
  facet_grid(rows = vars(as.factor(n)), cols = vars(as.factor(KL_div))) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))



# Model performance ####

## Bias ####
### Emission distribution ####
mean_abs_bias_emiss <- data.frame(true_m = 3,
                                  n = settings3[selected, 1], 
                                  n_t = settings3[selected, 2],
                                  n_dep = settings3[selected, 3], 
                                  KL_div = settings3[selected, 4], 
                                  abs_bias = NA,
                                  abs_bias2 = NA)

mean_rel_bias_emiss <- data.frame(true_m = 3,
                                  n = settings3[selected, 1], 
                                  n_t = settings3[selected, 2],
                                  n_dep = settings3[selected, 3], 
                                  KL_div = settings3[selected, 4], 
                                  rel_bias = NA,
                                  rel_bias2 = NA)

for(sc in 1:length(selected)){
  diff <- mapply('-', out[[sc]]$group_out_3st$emiss_mean$median, out[[sc]]$means_true, SIMPLIFY = FALSE)
  abs.diff <- sapply(diff, 'abs')
  mean_abs_bias_emiss[sc, 6] <- median(abs.diff)
  mean_abs_bias_emiss[sc, 7] <- median(sapply(diff, mean))
  rel.diff <- mapply('/', diff, out[[sc]]$means_true)
  mean_rel_bias_emiss[sc, 6] <- median(rel.diff)
  mean_rel_bias_emiss[sc, 7] <-  median(abs(rel.diff))
}

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


## Empirical standard error ####


## Coverage ####

Coverage_emiss <- data.frame(true_m = 3,
                             n = settings3[selected, 1], 
                             n_t = settings3[selected, 2],
                             n_dep = settings3[selected, 3], 
                             KL_div = settings3[selected, 4], 
                             cov_90 = NA, 
                             cov_95 = NA)

for(sc in 1:length(selected)){
  lower_90 <- lapply(mapply('<', out[[sc]]$group_out_3st$emiss_mean$quant_5, out[[sc]]$means_true, SIMPLIFY = FALSE), "*", 1)
  upper_90 <- lapply(mapply('>', out[[sc]]$group_out_3st$emiss_mean$quant_95, out[[sc]]$means_true, SIMPLIFY = FALSE), "*", 1)
  Coverage_emiss[sc, 6] <- sum(mapply("*", lower_90, upper_90)) / (n_sim * out[[sc]]$input$n_dep * true_m) * 100
  lower_95 <- lapply(mapply('<', out[[sc]]$group_out_3st$emiss_mean$quant_2.5, out[[sc]]$means_true, SIMPLIFY = FALSE), "*", 1)
  upper_95 <- lapply(mapply('>', out[[sc]]$group_out_3st$emiss_mean$quant_97.5, out[[sc]]$means_true, SIMPLIFY = FALSE), "*", 1)
  Coverage_emiss[sc, 7] <- sum(mapply("*", lower_95, upper_95)) / (n_sim * out[[sc]]$input$n_dep * true_m) * 100
}

ggplot(data = Coverage_emiss, mapping = aes(x = n_t, y = cov_95, group = as.factor(n), color = as.factor(n))) +
  geom_point() + 
  ylab("Coverage emission distribution") + 
  scale_color_discrete(name = "Number of\nsubjects") +
  scale_x_continuous(trans='log2', breaks = c(50, 100, 200, 400, 800)) +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep))) +
  xlab("number of observations per subject") +
  geom_hline(yintercept=c(90, 95), linetype="dashed", color = "grey") 




## State Decoding ####

### Mean proportion of states classified correctly ####

### Cohen's kappa ####

### Mean proportion of correct states with prob => 0.20 ####


# Convergence ####
