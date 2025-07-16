library(ggplot2)
library(cowplot)
library(combinat)


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


# ---------- Compute KLDs from Emmeke's Empirical Data ---


# ----- Bipolar (published) ------

## data ####
m <- 4
n_dep <- 12
bipolar_mean <- matrix(c( 8.9,  6.9, 16.2, 10.8, 33.0, 54.5,  9.6,  8.2, 38.9,  32.4,  30.5,  31.3,
                         20.2, 16.8, 25.8, 25.5, 34.1, 50.4, 33.8, 23.9, 44.5,  37.4,  36.9,  42.8,
                         49.6, 38.5, 45.1, 49.2, 58.3, 26.5, 35.4, 37.4, 31.8,  18.0,  20.9,  33.7,
                         41.4, 32.5, 31.3, 35.7, 55.3, 29.6, 18.5, 24.0, 22.9,  14.4,  13.3,  19.5), 
                       byrow = TRUE, nrow = m, ncol = n_dep)
colnames(bipolar_mean) <- c("down","dread","worry",
                            "inadequate", "tired", "content",
                            "agitated", "irritated", "switch_focus",
                            "extremely_well", "full_of_ideas", "thoughts_racing")
rownames(bipolar_mean) <- paste0("State_", 1:m)

bipolar_SD <-  matrix(c( 5.0,  4.9, 11.7,  6.7, 20.8, 12.9,  6.1,  5.7, 11.0,  12.1,  11.0,  10.1,
                        16.2, 19.5, 20.9, 19.6, 25.3, 17.6, 24.2, 22.3, 19.6,  20.7,  18.8,  21.1,
                        20.5, 23.7, 23.3, 22.2, 21.8, 17.0, 24.4, 23.9, 18.7,  12.2,  17.5,  21.3,
                        14.4, 19.6, 19.0, 18.7, 19.1, 14.9, 14.9, 17.8,  9.5,   7.9,   6.4,   9.3),
                        byrow = TRUE, nrow = m, ncol = n_dep)
colnames(bipolar_SD) <- c("down","dread","worry",
                            "inadequate", "tired", "content",
                            "agitated", "irritated", "switch_focus",
                            "extremely_well", "full_of_ideas", "thoughts_racing")
rownames(bipolar_SD) <- paste0("State_", 1:m)

## quick look ####
mu1 <- bipolar_mean[1,]
mu2 <- bipolar_mean[2,]
mu3 <- bipolar_mean[3,]
mu4 <- bipolar_mean[4,]

sig1 <- diag(n_dep) * bipolar_SD[1,]^2
sig2 <- diag(n_dep) * bipolar_SD[2,]^2
sig3 <- diag(n_dep) * bipolar_SD[3,]^2
sig4 <- diag(n_dep) * bipolar_SD[4,]^2


KLD(mu1=mu1, mu2=mu2, sig1=sig1, sig2=sig2)
KLD(mu1=mu1, mu2=mu3, sig1=sig1, sig2=sig3)
KLD(mu1=mu1, mu2=mu4, sig1=sig1, sig2=sig4)
KLD(mu1=mu2, mu2=mu3, sig1=sig2, sig2=sig3)
KLD(mu1=mu2, mu2=mu4, sig1=sig2, sig2=sig4)
KLD(mu1=mu3, mu2=mu4, sig1=sig3, sig2=sig4)

## more elaborate look ####
n_comp <- 6 
KLDs_bipolar <- data.frame(from_var = rep(c(rep(1,11), 7), each = n_comp),
                           to_var = rep(c(2:12,12), each = n_comp),
                           mu1 = c(rep(1,3), rep(2,2), 3),
                           mu2 = c(2, 3, 4, 3, 4, 4),
                           KLD = NA)

for(i in 1:length(KLDs_bipolar[,1])){
  vars <- KLDs_bipolar$from_var[i]:KLDs_bipolar$to_var[i]
  KLDs_bipolar$KLD[i] <- KLD(mu1 = bipolar_mean[KLDs_bipolar$mu1[i], vars], 
                             mu2 = bipolar_mean[KLDs_bipolar$mu2[i], vars], 
                             sig1 = (diag(length(vars)) * bipolar_SD[KLDs_bipolar$mu1[i], vars]^2), 
                             sig2 = (diag(length(vars)) * bipolar_SD[KLDs_bipolar$mu2[i], vars]^2))
}

KLDs_bipolar$var_from_to <- factor(paste0("var",KLDs_bipolar$from_var, "to", KLDs_bipolar$to_var), 
                                      levels = unique(paste0("var",KLDs_bipolar$from_var, "to", KLDs_bipolar$to_var)))

plot_mean_bipolar <- ggplot(data = KLDs_bipolar, aes(x = var_from_to, y = KLD)) +
  geom_bar(stat = "summary", fun = "mean", fill="steelblue") + 
  ggtitle("Mean KLDs for Bipolar") +
  ylab("Mean KLD over pairwise comparisons state 1 to 4") + ylim(0,28) +
  xlab("Selected variables") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
plot_mean_bipolar

plot_max_bipolar <- ggplot(data = KLDs_bipolar, aes(x = var_from_to, y = KLD)) +
  geom_bar(stat = "summary", fun = "max", fill="steelblue") + 
  ggtitle("Max KLDs for Bipolar") +
  ylab("Max KLD over pairwise comparisons state 1 to 4") + ylim(0,28) +
  xlab("Selected variables") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))
plot_max_bipolar

plot_grid(plot_mean_bipolar, plot_max_bipolar)



## Even more elaborate look ####
n_comp <- 6 

KLDs_bipolar2 <- data.frame(n_var = rep(1:12, each = n_comp),
                           mu1 = c(rep(1,3), rep(2,2), 3),
                           mu2 = c(2, 3, 4, 3, 4, 4),
                           KLD_mean = NA, 
                           KLD_max = NA)

for(i in 1:length(KLDs_bipolar2[,1])){
  if(KLDs_bipolar2$n_var[i] < n_dep){
    var_comb <- combn(1:n_dep, m = KLDs_bipolar2$n_var[i])
    n_comb <- dim(var_comb)[2]
    KLDs <- numeric(n_comb)
    
    for(j in 1:n_comb){
      KLDs[j] <- KLD(mu1 = bipolar_mean[KLDs_bipolar2$mu1[i], var_comb[,j]], 
                     mu2 = bipolar_mean[KLDs_bipolar2$mu2[i], var_comb[,j]], 
                     sig1 = (diag(dim(var_comb)[1]) * bipolar_SD[KLDs_bipolar2$mu1[i], var_comb[,j]]^2), 
                     sig2 = (diag(dim(var_comb)[1]) * bipolar_SD[KLDs_bipolar2$mu2[i], var_comb[,j]]^2))
    }
    KLDs_bipolar2$KLD_mean[i] <- mean(KLDs)
    KLDs_bipolar2$KLD_max[i] <- max(KLDs)
  } else {
    KLDs_bipolar2$KLD_mean[i] <- KLD(mu1 = bipolar_mean[KLDs_bipolar2$mu1[i], 1:n_dep], 
                                    mu2 = bipolar_mean[KLDs_bipolar2$mu2[i], 1:n_dep], 
                                    sig1 = (diag(n_dep) * bipolar_SD[KLDs_bipolar2$mu1[i], 1:n_dep]^2), 
                                    sig2 = (diag(n_dep) * bipolar_SD[KLDs_bipolar2$mu2[i], 1:n_dep]^2))
    KLDs_bipolar2$KLD_max[i] <- KLD(mu1 = bipolar_mean[KLDs_bipolar2$mu1[i], 1:n_dep], 
                                   mu2 = bipolar_mean[KLDs_bipolar2$mu2[i], 1:n_dep], 
                                   sig1 = (diag(n_dep) * bipolar_SD[KLDs_bipolar2$mu1[i], 1:n_dep]^2), 
                                   sig2 = (diag(n_dep) * bipolar_SD[KLDs_bipolar2$mu2[i], 1:n_dep]^2))
  }
  
}

KLDs_bipolar2$n_var <- factor(KLDs_bipolar2$n_var)
KLDs_bipolar2$state_comp <- factor(paste0("st", KLDs_bipolar2$mu1,"_vs_st", KLDs_bipolar2$mu2))


aggregate(KLDs_bipolar2, by = list(KLDs_bipolar2$n_var), FUN = mean)

### barplot visualization ####
plot_mean_mean_bipolar <- ggplot(data = KLDs_bipolar2, aes(x = n_var, y = KLD_mean)) +
  geom_bar(stat = "summary", fun = "mean", fill="steelblue") + 
  ggtitle("Mean mean KLDs for Bipolar") +
  ylab("Mean KLD over pairwise comparisons state 1 to 4") + ylim(0,28) +
  xlab("Number of selected variables from dataset") +
  theme_minimal() 
plot_mean_mean_bipolar

plot_max_mean_bipolar <- ggplot(data = KLDs_bipolar2, aes(x = n_var, y = KLD_mean)) +
  geom_bar(stat = "summary", fun = "max", fill="steelblue") + 
  ggtitle("Max mean KLDs for Bipolar") +
  ylab("Max KLD over pairwise comparisons state 1 to 4") + ylim(0,28) +
  xlab("Number of selected variables from dataset") +
  theme_minimal() 
plot_max_mean_bipolar

plot_max_bipolar <- ggplot(data = KLDs_bipolar2, aes(x = n_var, y = KLD_max)) +
  geom_bar(stat = "summary", fun = "max", fill="steelblue") + 
  ggtitle("Max KLDs for Bipolar") +
  ylab("Max KLD over pairwise comparisons state 1 to 4") + ylim(0,28) +
  xlab("Selected variables") +
  theme_minimal() 
plot_max_bipolar

plot_grid(plot_mean_mean_bipolar, plot_max_mean_bipolar)

### line visualization ####
plot_line_mean_bipolar <- ggplot(data = KLDs_bipolar2, aes(x = n_var, y = KLD_mean, group = state_comp, color = state_comp)) +
  geom_line() + 
  geom_point() + 
  ggtitle("Mean KLDs for Bipolar") +
  ylab("Mean KLD over variable combinations") + ylim(0,28) +
  xlab("Number of selected variables") +
  theme_minimal() + theme(legend.position="none")
plot_line_mean_bipolar

plot_line_max_bipolar <- ggplot(data = KLDs_bipolar2, aes(x = n_var, y = KLD_max, group = state_comp, color = state_comp)) +
  geom_line() + 
  geom_point() + 
  ggtitle("Max KLDs for Bipolar") +
  ylab("Max KLD over variable combinations") + ylim(0,28) +
  xlab("Number of selected variables") +
  theme_minimal() 
plot_line_max_bipolar

plot_grid(plot_line_mean_bipolar, plot_line_max_bipolar, rel_widths = c(0.8, 1))


#############
# ----- Suicidal crisis  ------
###########
m <- 4
n_dep <- 5
mu1 <- c(32, 44, 17, 42, 4)
mu2 <- c(52, 35, 30, 40, 48)
mu3 <- c(59, 33, 39, 42, 56)
mu4 <- c(72, 14, 49, 7, 66)

CAB_crisis_mean <- rbind(mu1, mu2, mu3, mu4)
colnames(CAB_crisis_mean) <- c("NegatiefAffect",
                               "Controle",
                               "Terugtrekken",
                               "BehoefteContact",
                               "Suicidaliteit")
rownames(CAB_crisis_mean) <- paste0("State_", 1:m)

sd1 <- c(13, 15, 16, 17,  2)
sd2 <- c(12, 14, 18, 21, 16)
sd3 <- c(15, 17, 26, 26, 20)
sd4 <- c(16,  9, 33,  4, 18)

CAB_crisis_SD <- rbind(sd1, sd2, sd3, sd4)
colnames(CAB_crisis_SD) <- c("NegatiefAffect",
                           "Controle",
                           "Terugtrekken",
                           "BehoefteContact",
                           "Suicidaliteit")
rownames(CAB_crisis_SD) <- paste0("State_", 1:m)

# more elobarte check 
n_comp <- 6 
KLDs_CAB_crisis <- data.frame(from_var = rep(c(rep(1,4)), each = n_comp),
                           to_var = rep(c(2:5), each = n_comp),
                           mu1 = c(rep(1,3), rep(2,2), 3),
                           mu2 = c(2, 3, 4, 3, 4, 4),
                           KLD = NA)

for(i in 1:length(KLDs_CAB_crisis[,1])){
  vars <- KLDs_CAB_crisis$from_var[i]:KLDs_CAB_crisis$to_var[i]
  KLDs_CAB_crisis$KLD[i] <- KLD(mu1 = CAB_crisis_mean[KLDs_CAB_crisis$mu1[i], vars], 
                             mu2 = CAB_crisis_mean[KLDs_CAB_crisis$mu2[i], vars], 
                             sig1 = (diag(length(vars)) * CAB_crisis_SD[KLDs_CAB_crisis$mu1[i], vars]^2), 
                             sig2 = (diag(length(vars)) * CAB_crisis_SD[KLDs_CAB_crisis$mu2[i], vars]^2))
}

KLDs_CAB_crisis$var_from_to <- factor(paste0("var",KLDs_CAB_crisis$from_var, "to", KLDs_CAB_crisis$to_var), 
                                   levels = unique(paste0("var",KLDs_CAB_crisis$from_var, "to", KLDs_CAB_crisis$to_var)))

plot_mean_CAB_crisis <- ggplot(data = KLDs_CAB_crisis, aes(x = var_from_to, y = KLD)) +
  geom_bar(stat = "summary", fun = "mean", fill="darkseagreen") + 
  ggtitle("Mean KLDs for CAB crisis") +
  ylab("Mean KLD over pairwise comparisons state 1 to 4") +
  xlab("Selected variables") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) 
plot_mean_CAB_crisis

plot_max_CAB_crisis <- ggplot(data = KLDs_CAB_crisis, aes(x = var_from_to, y = KLD)) +
  geom_bar(stat = "summary", fun = "max", fill="darkseagreen") + 
  ggtitle("Max KLDs for CAB crisis") +
  ylab("Max KLD over pairwise comparisons state 1 to 4") +
  xlab("Selected variables") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))
plot_max_CAB_crisis

plot_grid(plot_mean_CAB_crisis, plot_max_CAB_crisis)


plot_max_CAB_crisis <- ggplot(data = KLDs_CAB_crisis, aes(x = var_from_to, y = KLD)) +
  geom_bar(stat = "summary", fun = "max", fill="darkseagreen") + 
  ggtitle("Max KLDs for CAB crisis") +
  ylab("Max KLD over pairwise comparisons state 1 to 4") +
  xlab("Selected variables") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))
plot_max_CAB_crisis

plot_grid(plot_mean_CAB_crisis, plot_max_CAB_crisis)

## Even more elaborate look ####
n_comp <- 6 

KLDs_CAB_crisis2 <- data.frame(n_var = rep(1:n_dep, each = n_comp),
                             mu1 = c(rep(1,3), rep(2,2), 3),
                             mu2 = c(2, 3, 4, 3, 4, 4),
                             KLD_mean = NA, 
                             KLD_max = NA)

for(i in 1:length(KLDs_CAB_crisis2[,1])){
  if(KLDs_CAB_crisis2$n_var[i] < n_dep){
    var_comb <- combn(1:n_dep, m = KLDs_CAB_crisis2$n_var[i])
    n_comb <- dim(var_comb)[2]
    KLDs <- numeric(n_comb)
    
    for(j in 1:n_comb){
      KLDs[j] <- KLD(mu1 = CAB_crisis_mean[KLDs_CAB_crisis2$mu1[i], var_comb[,j]], 
                     mu2 = CAB_crisis_mean[KLDs_CAB_crisis2$mu2[i], var_comb[,j]], 
                     sig1 = (diag(dim(var_comb)[1]) * CAB_crisis_SD[KLDs_CAB_crisis2$mu1[i], var_comb[,j]]^2), 
                     sig2 = (diag(dim(var_comb)[1]) * CAB_crisis_SD[KLDs_CAB_crisis2$mu2[i], var_comb[,j]]^2))
    }
    KLDs_CAB_crisis2$KLD_mean[i] <- mean(KLDs)
    KLDs_CAB_crisis2$KLD_max[i] <- max(KLDs)
  } else {
    KLDs_CAB_crisis2$KLD_mean[i] <- KLD(mu1 = CAB_crisis_mean[KLDs_CAB_crisis2$mu1[i], 1:n_dep], 
                                      mu2 = CAB_crisis_mean[KLDs_CAB_crisis2$mu2[i], 1:n_dep], 
                                      sig1 = (diag(n_dep) * CAB_crisis_SD[KLDs_CAB_crisis2$mu1[i], 1:n_dep]^2), 
                                      sig2 = (diag(n_dep) * CAB_crisis_SD[KLDs_CAB_crisis2$mu2[i], 1:n_dep]^2))
    KLDs_CAB_crisis2$KLD_max[i] <- KLD(mu1 = CAB_crisis_mean[KLDs_CAB_crisis2$mu1[i], 1:n_dep], 
                                     mu2 = CAB_crisis_mean[KLDs_CAB_crisis2$mu2[i], 1:n_dep], 
                                     sig1 = (diag(n_dep) * CAB_crisis_SD[KLDs_CAB_crisis2$mu1[i], 1:n_dep]^2), 
                                     sig2 = (diag(n_dep) * CAB_crisis_SD[KLDs_CAB_crisis2$mu2[i], 1:n_dep]^2))
  }
  
}

KLDs_CAB_crisis2$n_var <- factor(KLDs_CAB_crisis2$n_var)
KLDs_CAB_crisis2$state_comp <- factor(paste0("st", KLDs_CAB_crisis2$mu1,"_vs_st", KLDs_CAB_crisis2$mu2))

### barplot visualization ####
plot_mean_mean_CAB_crisis <- ggplot(data = KLDs_CAB_crisis2, aes(x = n_var, y = KLD_mean)) +
  geom_bar(stat = "summary", fun = "mean", fill="steelblue") + 
  ggtitle("Mean mean KLDs for CAB crisis") +
  ylab("Mean KLD over pairwise comparisons state 1 to 4") + ylim(0,28) +
  xlab("Number of selected variables from dataset") +
  theme_minimal() 
plot_mean_mean_CAB_crisis

plot_max_mean_CAB_crisis <- ggplot(data = KLDs_CAB_crisis2, aes(x = n_var, y = KLD_mean)) +
  geom_bar(stat = "summary", fun = "max", fill="steelblue") + 
  ggtitle("Max mean KLDs for CAB crisis") +
  ylab("Max KLD over pairwise comparisons state 1 to 4") + ylim(0,28) +
  xlab("Number of selected variables from dataset") +
  theme_minimal() 
plot_max_mean_CAB_crisis

plot_max_CAB_crisis <- ggplot(data = KLDs_CAB_crisis2, aes(x = n_var, y = KLD_max)) +
  geom_bar(stat = "summary", fun = "max", fill="steelblue") + 
  ggtitle("Max KLDs for CAB crisis") +
  ylab("Max KLD over pairwise comparisons state 1 to 4") + ylim(0,28) +
  xlab("Selected variables") +
  theme_minimal() 
plot_max_CAB_crisis

plot_grid(plot_mean_mean_CAB_crisis, plot_max_mean_CAB_crisis)

### line visualization ####
plot_line_mean_CAB_crisis <- ggplot(data = KLDs_CAB_crisis2, aes(x = n_var, y = KLD_mean, group = state_comp, color = state_comp)) +
  geom_line() + 
  geom_point() + 
  ggtitle("Mean KLDs for CAB crisis") +
  ylab("Mean KLD over variable combinations") + 
  xlab("Number of selected variables") +
  theme_minimal() + theme(legend.position="none")
plot_line_mean_CAB_crisis

plot_line_max_CAB_crisis <- ggplot(data = KLDs_CAB_crisis2, aes(x = n_var, y = KLD_max, group = state_comp, color = state_comp)) +
  geom_line() + 
  geom_point() + 
  ggtitle("Max KLDs for CAB crisis") +
  ylab("Max KLD over variable combinations") + 
  xlab("Number of selected variables") +
  theme_minimal() 
plot_line_max_CAB_crisis

plot_grid(plot_line_mean_CAB_crisis, plot_line_max_CAB_crisis, rel_widths = c(0.8, 1))

#########
# ----- Psychotic substances (submitted) ------
########
m <- 4
n_dep <- 6
mu1 <- c(25, 17, 7, 7, 6, 5)
mu2 <- c(39, 35, 20, 20, 4, 4)
mu3 <- c(37, 31, 16, 19, 22, 12)
mu4 <- c(58, 57, 41, 40, 32, 29)

psychotic_mean <- rbind(mu1, mu2, mu3, mu4)
colnames(psychotic_mean) <- c("Stressed",
                              "Down",
                              "Suspicious",
                              "Paranoia",
                              "Broadcasting",
                              "External_Control" )
rownames(psychotic_mean) <- paste0("State_", 1:m)

sd1 <- c(4, 3, 2, 2, 2, 2)
sd2 <- c(4, 3, 2, 2, 2, 2)
sd3 <- c(4, 4, 4, 4, 4, 4)
sd4 <- c(5, 5, 5, 5, 5, 5)

psychotic_SD <- rbind(sd1, sd2, sd3, sd4)
colnames(psychotic_SD) <- c("Stressed",
                             "Down",
                             "Suspicious",
                             "Paranoia",
                             "Broadcasting",
                             "External_Control" )
rownames(psychotic_SD) <- paste0("State_", 1:m)


# more elobarte check 
n_comp <- 6 
KLDs_psychotic <- data.frame(from_var = rep(c(rep(1,5)), each = n_comp),
                             to_var = rep(c(2:6), each = n_comp),
                             mu1 = c(rep(1,3), rep(2,2), 3),
                             mu2 = c(2, 3, 4, 3, 4, 4),
                             KLD = NA)

for(i in 1:length(KLDs_psychotic[,1])){
  vars <- KLDs_psychotic$from_var[i]:KLDs_psychotic$to_var[i]
  KLDs_psychotic$KLD[i] <- KLD(mu1 = psychotic_mean[KLDs_psychotic$mu1[i], vars], 
                               mu2 = psychotic_mean[KLDs_psychotic$mu2[i], vars], 
                               sig1 = (diag(length(vars)) * psychotic_SD[KLDs_psychotic$mu1[i], vars]^2), 
                               sig2 = (diag(length(vars)) * psychotic_SD[KLDs_psychotic$mu2[i], vars]^2))
}

KLDs_psychotic$var_from_to <- factor(paste0("var",KLDs_psychotic$from_var, "to", KLDs_psychotic$to_var), 
                                     levels = unique(paste0("var",KLDs_psychotic$from_var, "to", KLDs_psychotic$to_var)))

plot_mean_psychotic <- ggplot(data = KLDs_psychotic, aes(x = var_from_to, y = KLD)) +
  geom_bar(stat = "summary", fun = "mean", fill="darkorchid4") + 
  ggtitle("Mean KLDs for psychotic") +
  ylab("Mean KLD over pairwise comparisons state 1 to 4") +
  xlab("Selected variables") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) 
plot_mean_psychotic

plot_max_psychotic <- ggplot(data = KLDs_psychotic, aes(x = var_from_to, y = KLD)) +
  geom_bar(stat = "summary", fun = "max", fill="darkorchid4") + 
  ggtitle("Max KLDs for psychotic") +
  ylab("Max KLD over pairwise comparisons state 1 to 4") +
  xlab("Selected variables") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))
plot_max_psychotic

plot_grid(plot_mean_psychotic, plot_max_psychotic)

## Even more elaborate look ####
n_comp <- 6 

KLDs_psychotic2 <- data.frame(n_var = rep(1:n_dep, each = n_comp),
                               mu1 = c(rep(1,3), rep(2,2), 3),
                               mu2 = c(2, 3, 4, 3, 4, 4),
                               KLD_mean = NA, 
                               KLD_max = NA)

for(i in 1:length(KLDs_psychotic2[,1])){
  if(KLDs_psychotic2$n_var[i] < n_dep){
    var_comb <- combn(1:n_dep, m = KLDs_psychotic2$n_var[i])
    n_comb <- dim(var_comb)[2]
    KLDs <- numeric(n_comb)
    
    for(j in 1:n_comb){
      KLDs[j] <- KLD(mu1 = psychotic_mean[KLDs_psychotic2$mu1[i], var_comb[,j]], 
                     mu2 = psychotic_mean[KLDs_psychotic2$mu2[i], var_comb[,j]], 
                     sig1 = (diag(dim(var_comb)[1]) * psychotic_SD[KLDs_psychotic2$mu1[i], var_comb[,j]]^2), 
                     sig2 = (diag(dim(var_comb)[1]) * psychotic_SD[KLDs_psychotic2$mu2[i], var_comb[,j]]^2))
    }
    KLDs_psychotic2$KLD_mean[i] <- mean(KLDs)
    KLDs_psychotic2$KLD_max[i] <- max(KLDs)
  } else {
    KLDs_psychotic2$KLD_mean[i] <- KLD(mu1 = psychotic_mean[KLDs_psychotic2$mu1[i], 1:n_dep], 
                                        mu2 = psychotic_mean[KLDs_psychotic2$mu2[i], 1:n_dep], 
                                        sig1 = (diag(n_dep) * psychotic_SD[KLDs_psychotic2$mu1[i], 1:n_dep]^2), 
                                        sig2 = (diag(n_dep) * psychotic_SD[KLDs_psychotic2$mu2[i], 1:n_dep]^2))
    KLDs_psychotic2$KLD_max[i] <- KLD(mu1 = psychotic_mean[KLDs_psychotic2$mu1[i], 1:n_dep], 
                                       mu2 = psychotic_mean[KLDs_psychotic2$mu2[i], 1:n_dep], 
                                       sig1 = (diag(n_dep) * psychotic_SD[KLDs_psychotic2$mu1[i], 1:n_dep]^2), 
                                       sig2 = (diag(n_dep) * psychotic_SD[KLDs_psychotic2$mu2[i], 1:n_dep]^2))
  }
  
}

KLDs_psychotic2$n_var <- factor(KLDs_psychotic2$n_var)
KLDs_psychotic2$state_comp <- factor(paste0("st", KLDs_psychotic2$mu1,"_vs_st", KLDs_psychotic2$mu2))

### barplot visualization ####
plot_mean_mean_psychotic <- ggplot(data = KLDs_psychotic2, aes(x = n_var, y = KLD_mean)) +
  geom_bar(stat = "summary", fun = "mean", fill="steelblue") + 
  ggtitle("Mean mean KLDs for psychotic") +
  ylab("Mean KLD over pairwise comparisons state 1 to 4") + ylim(0,28) +
  xlab("Number of selected variables from dataset") +
  theme_minimal() 
plot_mean_mean_psychotic

plot_max_mean_psychotic <- ggplot(data = KLDs_psychotic2, aes(x = n_var, y = KLD_mean)) +
  geom_bar(stat = "summary", fun = "max", fill="steelblue") + 
  ggtitle("Max mean KLDs for psychotic") +
  ylab("Max KLD over pairwise comparisons state 1 to 4") + ylim(0,28) +
  xlab("Number of selected variables from dataset") +
  theme_minimal() 
plot_max_mean_psychotic

plot_max_psychotic <- ggplot(data = KLDs_psychotic2, aes(x = n_var, y = KLD_max)) +
  geom_bar(stat = "summary", fun = "max", fill="steelblue") + 
  ggtitle("Max KLDs for psychotic") +
  ylab("Max KLD over pairwise comparisons state 1 to 4") + ylim(0,28) +
  xlab("Selected variables") +
  theme_minimal() 
plot_max_psychotic

plot_grid(plot_mean_mean_psychotic, plot_max_mean_psychotic)

### line visualization ####
plot_line_mean_psychotic <- ggplot(data = KLDs_psychotic2, aes(x = n_var, y = KLD_mean, group = state_comp, color = state_comp)) +
  geom_line() + 
  geom_point() + 
  ggtitle("Mean KLDs for psychotic") +
  ylab("Mean KLD over variable combinations") + 
  xlab("Number of selected variables") +
  theme_minimal() + theme(legend.position="none")
plot_line_mean_psychotic

plot_line_max_psychotic <- ggplot(data = KLDs_psychotic2, aes(x = n_var, y = KLD_max, group = state_comp, color = state_comp)) +
  geom_line() + 
  geom_point() + 
  ggtitle("Max KLDs for psychotic") +
  ylab("Max KLD over variable combinations") + 
  xlab("Number of selected variables") +
  theme_minimal() 
plot_line_max_psychotic

plot_grid(plot_line_mean_psychotic, plot_line_max_psychotic, rel_widths = c(0.8, 1))

aggregate(KLDs_psychotic2, by = list(KLDs_psychotic2$n_var), FUN = mean)

#########


# ----- Mindfulness (toy emp example) ------
mu1 <- c(68, 49, 64, 65, 4, 5, 9, 6)
mu2 <- c(61, 43, 54, 54, 17, 19, 25, 20)
mu3 <- c(36, 22, 34, 31, 34, 24, 47, 34)
d <- length(mu1)
sig1 <- diag(d) * c(15, 24, 21, 18, 4, 4, 8, 5)^2
sig2 <- diag(d) * c(16, 21, 22, 19, 17, 20, 19, 18)^2
sig3 <- diag(d) * c(19, 16, 22, 19, 28, 26, 25, 26)^2
# 3 Comparisons
KLD(mu1=mu1, mu2=mu2, sig1=sig1, sig2=sig2)
KLD(mu1=mu1, mu2=mu3, sig1=sig1, sig2=sig3)
KLD(mu1=mu2, mu2=mu3, sig1=sig2, sig2=sig3)

