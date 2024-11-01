library(mHMMbayes)
library(ggplot2)
library(reshape2)

# Simulstion settings
settings <- read.csv("sim_scenarios.csv")
true_m  <- 3

# scenario <- 1
# KL_div  <- settings[scenario, 4]
# n_t     <- settings[scenario, 2]
# n       <- settings[scenario, 1]
# n_dep   <- settings[scenario, 3]
# out_file <- paste0("out_n", n, "_nt", n_t, "_ndep", n_dep, "_", true_m, "st_KLD", KL_div)

## Reading in output files ####
selected <- sort(c(seq(1,39,3), seq(2,39,3), seq(43,69,3), seq(44,69,3), 76, 79, 74, 80, 
                   seq(82, 90, 3), seq(83, 90, 3)))
           
files <- paste0("out_n", settings[selected, 1], "_nt", settings[selected, 2], "_ndep", settings[selected, 3], "_", true_m, "st_KLD", settings[selected, 4])

# out <- readRDS(paste0("out/", out_file, ".rds")) 

out <- vector("list", length(files))
names(out) <- files

for(sc in 1:length(files)){
  out[[sc]] <- readRDS(paste0("out2/", files[sc], ".rds")) 
}

n_sim <- out[[1]]$input$n_sim
J <- out[[1]]$input$J
burn_in <- out[[1]]$input$burn_in

# Model selection ####

out[[1]]$selection_out

AICc_mean <- data.frame(true_m = 3,
                        n = settings[selected, 1], 
                        n_t = settings[selected, 2],
                        n_dep = settings[selected, 3], 
                        KL_div = settings[selected, 4], 
                        st1 = NA, 
                        st2 = NA, 
                        st3 = NA, 
                        st4 = NA)
min_AICc <- data.frame(true_m = 3,
                       n = settings[selected, 1], 
                       n_t = settings[selected, 2],
                       n_dep = settings[selected, 3], 
                       KL_div = settings[selected, 4], 
                       prop_st1 = NA,
                       prop_st2 = NA,
                       prop_st3 = NA,
                       prop_st4 = NA)

AIC_mean <- data.frame(true_m = 3,
                        n = settings[selected, 1], 
                        n_t = settings[selected, 2],
                        n_dep = settings[selected, 3], 
                        KL_div = settings[selected, 4], 
                        st1 = NA, 
                        st2 = NA, 
                        st3 = NA, 
                        st4 = NA)
min_AIC <- data.frame(true_m = 3,
                       n = settings[selected, 1], 
                       n_t = settings[selected, 2],
                       n_dep = settings[selected, 3], 
                       KL_div = settings[selected, 4], 
                       prop_st1 = NA,
                       prop_st2 = NA,
                       prop_st3 = NA,
                       prop_st4 = NA)

for(sc in 1:length(selected)){
  AICc_mean[sc, 6:9] <- apply(out[[sc]]$selection_out$AICc_mean, 2, mean)
}

AICc_mean <- melt(AICc_mean, id.vars = colnames(AICc_mean)[1:5], variable.name = "state")

ggplot(data = AICc_mean, mapping = aes(x = state, y = value, group = as.factor(n_t), color = as.factor(n_t))) +
  geom_line(stat = "summary", fun = "mean") + ylab("AICc") +
  facet_grid(rows = vars(as.factor(n)), cols = vars(as.factor(KL_div)))

  
for(sc in 1:length(selected)){
  min_AICc[sc, 6:9] <- table(factor(apply(out[[sc]]$selection_out$AICc_mean, 1, which.min), levels = 1:4)) / n_sim
}

gg_min_AICc <- melt(min_AICc, id.vars = colnames(min_AICc)[1:5], variable.name = "state")

ggplot(data = gg_min_AICc, mapping = aes(x = state, y = value, fill = as.factor(n_t))) +
  geom_bar(stat = "summary", fun = "mean", color="black", position=position_dodge()) +
  facet_grid(rows = vars(as.factor(n)), cols = vars(as.factor(KL_div)))


for(sc in 1:length(selected)){
  AIC_mean[sc, 6:9] <- apply(out[[sc]]$selection_out$AIC_mean, 2, mean)
}

for(sc in 1:length(selected)){
  min_AIC[sc, 6:9] <- table(factor(apply(out[[sc]]$selection_out$AIC_mean, 1, which.min), levels = 1:4)) / n_sim
}

gg_min_AIC <- melt(min_AIC, id.vars = colnames(min_AIC)[1:5], variable.name = "state")

ggplot(data = gg_min_AIC, mapping = aes(x = state, y = value, fill = as.factor(n_t))) +
  geom_bar(stat = "summary", fun = "mean", color="black", position=position_dodge()) +
  facet_grid(rows = vars(as.factor(n)), cols = vars(as.factor(KL_div)))


# Model performance ####

## Bias ####
### Emission distribution ####
mean_abs_bias_emiss <- data.frame(true_m = 3,
                      n = settings[selected, 1], 
                      n_t = settings[selected, 2],
                      n_dep = settings[selected, 3], 
                      KL_div = settings[selected, 4], 
                      abs_bias = NA,
                      abs_bias2 = NA)

mean_rel_bias_emiss <- data.frame(true_m = 3,
                              n = settings[selected, 1], 
                              n_t = settings[selected, 2],
                              n_dep = settings[selected, 3], 
                              KL_div = settings[selected, 4], 
                              rel_bias = NA,
                              rel_bias2 = NA)

for(sc in 1:length(selected)){
  diff <- mapply('-', out[[sc]]$group_out_3st$emiss_mean$median, out[[sc]]$means_true, SIMPLIFY = FALSE)
  abs.diff <- sapply(diff, 'abs')
  mean_abs_bias_emiss[sc, 6] <- mean(abs.diff)
  mean_abs_bias_emiss[sc, 7] <- mean(sapply(diff, mean))
  rel.diff <- mapply('/', diff, out[[sc]]$means_true)
  mean_rel_bias_emiss[sc, 6] <- mean(rel.diff)
  mean_rel_bias_emiss[sc, 7] <-  mean(abs(rel.diff))
}

ggplot(data = mean_abs_bias_emiss, mapping = aes(x = n_t, y = abs_bias, group = as.factor(n), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_line() + 
  ylab("mean absolute bias") +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep)))


ggplot(data = mean_rel_bias_emiss, mapping = aes(x = n_t, y = rel_bias, group = as.factor(n), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_line() + 
  ylab("mean relative bias") +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep)))

ggplot(data = mean_rel_bias_emiss, mapping = aes(x = n_t, y = rel_bias2, group = as.factor(n), color = as.factor(n))) +
  # geom_line(stat = "summary", fun = "mean") + 
  geom_line() + 
  ylab("mean absolute relative bias") +
  facet_grid(rows = vars(as.factor(KL_div)), cols = vars(as.factor(n_dep)))


out[[3]]$input

out[[52]]$group_out_3st$emiss_mean$median
out[[30]]$means_true
length(files)

abs.diff <- mapply('-', out[[30]]$group_out_3st$emiss_mean$median, out[[30]]$means_true, SIMPLIFY = FALSE)

rel.diff <-  mapply('/', diff, out[[30]]$means_true, SIMPLIFY = FALSE)


sapply(out[[4]]$means_true, mean)


## Empirical standard error ####


## Coverage ####

### 95 CrI ####

### 99 CrI ####


## State Decoding ####

### Mean proportion of states classified correctly ####

### Cohen's kappa ####

### Mean proportion of correct states with prob => 0.20 ####


# Convergence ####
