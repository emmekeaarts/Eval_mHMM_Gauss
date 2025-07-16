library(mHMMbayes)
library(depmixS4)
library(combinat)

argv = commandArgs(trailingOnly=TRUE)
print(argv)
num_argv = as.numeric(argv)

settings <- read.csv("sim_scenarios_V4.csv")


# KL_div  <- 5
# n_t     <- 50
# n       <- 50
# n_dep   <- 4

KL_div  <- settings[num_argv, 4]
n_t     <- settings[num_argv, 2]
n       <- settings[num_argv, 1]
n_dep   <- settings[num_argv, 3]

true_m  <- 3
max_m   <- 4
n_sim   <- 25

J       <- 2000
burn_in <- 500

out_file <- paste0("out_n", n, "_nt", n_t, "_ndep", n_dep, "_", true_m, "st_KLD", KL_div)
out_file
input <- list(KL_div = KL_div, n_t = n_t, n = n, true_m = true_m, n_dep = n_dep, n_sim = n_sim, J = J, burn_in = burn_in)

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
    
    return(max(v_error))
    
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

means2 <- matrix(, nrow = n_sim, ncol = true_m * n_dep)
colnames(means2) <- paste0(rep(paste0("dep", 1:n_dep), each = true_m), "_state", 1:true_m)
means3 <- rep(list(matrix(, ncol = true_m, nrow = n_dep)), n_sim)
emiss_distr <- rep(list(rep(list(matrix(, ncol = 2, nrow = true_m)), n_dep)), n_sim)

if(n_dep == 8){
  init_range <- c(-2, 2)
} else {
  init_range <- c(-4, 4)
}

set.seed(129864 + num_argv)
for(s in 1:n_sim){
  means <- makeGGM_KLD(K = true_m, p = n_dep, target_KLD = KL_div, sig = 1, init_range = init_range)
  means3[[s]][] <- matrix(unlist(means$mu), ncol = true_m, nrow = n_dep)
  means2[s,] <- as.vector(t(means3[[s]]))
  for(q in 1:n_dep){
    emiss_distr[[s]][[q]][,1] <- means2[s,((q-1)*true_m + 1) : ((q-1)*true_m + true_m)]
    emiss_distr[[s]][[q]][,2] <- 1
  }
}
means_true <- rep(list(matrix(,ncol = true_m, nrow = n_sim, 
                              dimnames = list(1:n_sim, paste0("st", 1:true_m)))), n_dep)
for(q in 1:n_dep){
  for(i in 1:true_m){
    means_true[[q]][,i] <- unlist(sapply(sapply(emiss_distr, "[[", q, simplify = FALSE), "[[", i))
  }
}

# generating simulated dataset 
gamma <- matrix(c(0.7, 0.2, 0.1, 
                  0.1, 0.8, 0.1, 
                  0.1, 0.1, 0.8), ncol = true_m, byrow = TRUE)

data_cont <- vector(mode='list', length=n_sim)
depmix_means <- rep(list(vector(mode='list', length=n_sim)), max_m)

state_permutations <- matrix(unlist(permn(1:true_m)), byrow = TRUE, ncol = 3)
n_perm <- dim(state_permutations)[1]
SEE <- numeric(n_perm)

subj_par <- rep(list(list(gamma = matrix(, nrow = n, ncol = true_m * true_m, 
                                         dimnames = list(1:n, paste0("s", rep(1:true_m, each = true_m), "to_s", 1:true_m))), 
                          emiss_mean = rep(list(matrix(, nrow = n, ncol = true_m, 
                                                       dimnames = list(1:n, paste0("st", 1:true_m)))), n_dep))), n_sim)

true_states <- cbind(rep(1:n, each = n_t), matrix(,nrow = n_t*n, ncol = n_sim))
colnames(true_states) <- c("ID", paste0("sim", 1:n_sim))


set.seed(207956290 + num_argv)
for(s in 1:n_sim){
  # simulate n_sim datasets
  data_cont[[s]] <- sim_mHMM(n_t = n_t, n = n, data_distr = 'continuous', gen = list(m = true_m, n_dep = n_dep),
                             gamma = gamma, emiss_distr = emiss_distr[[s]], var_gamma = .05, var_emiss = c(rep(.2^2, n_dep)),
                             return_ind_par = TRUE)
  colnames(data_cont[[s]]$obs) <- c("subj", paste0("dep", 1:n_dep)) 
  
  # save true states
  true_states[,(s + 1)] <- data_cont[[s]]$states[,2]
  
  # save subject level parameters
  subj_par[[s]]$gamma[] <- matrix(unlist(lapply(data_cont[[s]]$subject_gamma, t)), byrow = TRUE, ncol = true_m * true_m) 
  for(q in 1:n_dep){
    for(i in 1:true_m){
      subj_par[[s]]$emiss_mean[[q]][,i] <- unlist(sapply(sapply(data_cont[[s]]$subject_emiss, "[[", q, simplify = FALSE), "[[", i))
    }
  }
  
  # fit depmix models to possible number of states to obtain reasonable starting values and prior means
  # only for m =>2
  for(i in 2:max_m){
    depmix_model <- depmix(mapply(as.formula, paste0("dep", 1:n_dep, " ~ 1")), 
                           data = data.frame(data_cont[[s]]$obs), 
                           nstates = i,
                           family = rep(list(gaussian()), n_dep))
    depmix_fit <- NULL
    attempt <- 1
    while(is.null(depmix_fit) && attempt <= 100){
      attempt <- attempt + 1
      try(depmix_fit <- fit(depmix_model), silent = TRUE)
    }
    depmix_means[[i]][[s]] <- matrix(t(summary(depmix_fit, "response")[,seq(1,n_dep * 2, 2)]), ncol = i, nrow = n_dep)
  }
  # re-arrange means in correct order only for true state model
  for(p in 1:n_perm){
    SEE[p] <- sum((means3[[s]] - depmix_means[[true_m]][[s]][,state_permutations[p,]])^2)
  }
  correct_order <- state_permutations[which.min(SEE),]
  depmix_means[[true_m]][[s]] <- depmix_means[[true_m]][[s]][,correct_order]
  
  # for one state solution, use variable means as starting values
  depmix_means[[1]][[s]] <- matrix(apply(data.frame(data_cont[[s]]$obs)[,2:(n_dep + 1)], 2, mean), ncol = 1)
}

means_true2 <-lapply(means_true,"[",c(20), TRUE)
data_cont2 <- list(data_cont[[20]])
depmix_means2 <- list(list(depmix_means[[true_m]][[20]]))
subj_par2 <- list(subj_par[[20]])

# Specify hyper-prior for the continuous emission distribution AND transition probs
conv_runs <- c(20)
n_conv_runs <- length(conv_runs)

manual_prior_emiss <- list(vector(mode='list', length=n_conv_runs))

manual_prior_gamma <- prior_gamma(m = 3, gamma_mu0 = list(prob_to_int(matrix(c(0.6, 0.2, 0.2), nrow = 1)),
                                                          prob_to_int(matrix(c(0.2, 0.6, 0.2), nrow = 1)),
                                                          prob_to_int(matrix(c(0.2, 0.2, 0.6), nrow = 1)))) 


for(s in 1:n_conv_runs){
  manual_prior_emiss[[s]] <- 
    prior_emiss_cont(gen = list(m = true_m, n_dep = n_dep),
                     emiss_mu0 = lapply(sapply(apply(
                       depmix_means2[[1]][[s]], 1, list), matrix), matrix, nrow = 1),
                     emiss_K0 = rep(list(1), n_dep),
                     emiss_V =  rep(list(rep(.1^2, true_m)), n_dep),
                     emiss_nu = rep(list(1), n_dep),
                     emiss_a0 = rep(list(rep(3, true_m)), n_dep),
                     emiss_b0 =  rep(list(rep(1, true_m)), n_dep))
}


# specify starting values for the models with 1 to max_m states
gamma_start2 <- diag(0.7, true_m)
gamma_start2[lower.tri(gamma_start2) | upper.tri(gamma_start2)] <- 
  (1 - 0.7) / (true_m - 1)
gamma_start3 <- diag(0.9, true_m)
gamma_start3[lower.tri(gamma_start3) | upper.tri(gamma_start3)] <- 
  (1 - 0.9) / (true_m - 1)
emiss_start2 <- list(vector(mode='list', length=n_conv_runs))
emiss_start3 <- list(vector(mode='list', length=n_conv_runs))

foo <- function(x){
  cbind(x, 1)
}

set.seed(25396  + num_argv)
for(s in 1:n_conv_runs){
  emiss_start2[[s]] <- lapply(sapply(apply(
    depmix_means2[[1]][[s]], 1, list), matrix), matrix, ncol = 1)
  emiss_start3[[s]] <- lapply(sapply(apply(
    depmix_means2[[1]][[s]], 1, list), matrix), matrix, ncol = 1)
  emiss_start2[[s]] <- 
    lapply(emiss_start2[[s]], function(x){x + runif(true_m, -0.2, 0.2)})
  emiss_start3[[s]] <- 
    lapply(emiss_start3[[s]], function(x){x + runif(true_m, -0.2, 0.2)})
  emiss_start2[[s]] <- lapply(emiss_start2[[s]], foo)
  emiss_start3[[s]] <- lapply(emiss_start3[[s]], foo)
}


group_out_3st_conv1 <- list(trans_prob = 
                              list(mean =  matrix(, nrow = n_conv_runs, ncol = true_m * true_m, 
                                                  dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m), "to_s", 1:true_m))),
                                   var = matrix(, nrow = n_conv_runs, ncol = true_m * true_m, 
                                                dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m), "to_s", 1:true_m))),
                                   mode = matrix(, nrow = n_conv_runs, ncol = true_m * true_m, 
                                                 dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m), "to_s", 1:true_m))),
                                   median = matrix(, nrow = n_conv_runs, ncol = true_m * true_m, 
                                                   dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m), "to_s", 1:true_m))), 
                                   quant_25 = matrix(, nrow = n_conv_runs, ncol = true_m * true_m, 
                                                     dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m), "to_s", 1:true_m))), 
                                   quant_75 = matrix(, nrow = n_conv_runs, ncol = true_m * true_m, 
                                                     dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m), "to_s", 1:true_m))), 
                                   quant_5 = matrix(, nrow = n_conv_runs, ncol = true_m * true_m, 
                                                    dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m), "to_s", 1:true_m))), 
                                   quant_95 = matrix(, nrow = n_conv_runs, ncol = true_m * true_m, 
                                                     dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m), "to_s", 1:true_m))), 
                                   quant_2.5 = matrix(, nrow = n_conv_runs, ncol = true_m * true_m, 
                                                      dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m), "to_s", 1:true_m))), 
                                   quant_97.5 = matrix(, nrow = n_conv_runs, ncol = true_m * true_m, 
                                                       dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m), "to_s", 1:true_m)))
                              ), 
                            trans_interc = 
                              list(mean =  matrix(, nrow = n_conv_runs, ncol = true_m * (true_m - 1), 
                                                  dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m - 1), "to_s", 2:true_m))),
                                   var = matrix(, nrow = n_conv_runs, ncol = true_m * (true_m - 1), 
                                                dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m - 1), "to_s", 2:true_m))),
                                   mode = matrix(, nrow = n_conv_runs, ncol = true_m * (true_m - 1), 
                                                 dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m - 1), "to_s", 2:true_m))),
                                   median = matrix(, nrow = n_conv_runs, ncol = true_m * (true_m - 1), 
                                                   dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m - 1), "to_s", 2:true_m))),
                                   quant_25 = matrix(, nrow = n_conv_runs, ncol = true_m * (true_m - 1), 
                                                     dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m - 1), "to_s", 2:true_m))),
                                   quant_75 = matrix(, nrow = n_conv_runs, ncol = true_m * (true_m - 1), 
                                                     dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m - 1), "to_s", 2:true_m))),
                                   quant_5 = matrix(, nrow = n_conv_runs, ncol = true_m * (true_m - 1), 
                                                    dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m - 1), "to_s", 2:true_m))),
                                   quant_95 = matrix(, nrow = n_conv_runs, ncol = true_m * (true_m - 1), 
                                                     dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m - 1), "to_s", 2:true_m))),
                                   quant_2.5 = matrix(, nrow = n_conv_runs, ncol = true_m * (true_m - 1), 
                                                      dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m - 1), "to_s", 2:true_m))),
                                   quant_97.5 = matrix(, nrow = n_conv_runs, ncol = true_m * (true_m - 1), 
                                                       dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m - 1), "to_s", 2:true_m)))
                              ), 
                            trans_betw_subj = 
                              list(mean =  matrix(, nrow = n_conv_runs, ncol = true_m * (true_m - 1), 
                                                  dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m - 1), "to_s", 2:true_m))),
                                   var = matrix(, nrow = n_conv_runs, ncol = true_m * (true_m - 1), 
                                                dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m - 1), "to_s", 2:true_m))),
                                   mode = matrix(, nrow = n_conv_runs, ncol = true_m * (true_m - 1), 
                                                 dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m - 1), "to_s", 2:true_m))),
                                   median = matrix(, nrow = n_conv_runs, ncol = true_m * (true_m - 1), 
                                                   dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m - 1), "to_s", 2:true_m))),
                                   quant_25 = matrix(, nrow = n_conv_runs, ncol = true_m * (true_m - 1), 
                                                     dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m - 1), "to_s", 2:true_m))),
                                   quant_75 = matrix(, nrow = n_conv_runs, ncol = true_m * (true_m - 1), 
                                                     dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m - 1), "to_s", 2:true_m))),
                                   quant_5 = matrix(, nrow = n_conv_runs, ncol = true_m * (true_m - 1), 
                                                    dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m - 1), "to_s", 2:true_m))),
                                   quant_95 = matrix(, nrow = n_conv_runs, ncol = true_m * (true_m - 1), 
                                                     dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m - 1), "to_s", 2:true_m))),
                                   quant_2.5 = matrix(, nrow = n_conv_runs, ncol = true_m * (true_m - 1), 
                                                      dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m - 1), "to_s", 2:true_m))),
                                   quant_97.5 = matrix(, nrow = n_conv_runs, ncol = true_m * (true_m - 1), 
                                                       dimnames = list(1:n_conv_runs,paste0("s", rep(1:true_m, each = true_m - 1), "to_s", 2:true_m)))
                              ), 
                            emiss_mean = 
                              list(mean = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                          dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep),
                                   var = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                         dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep), 
                                   mode = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                          dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep), 
                                   median = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                            dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep), 
                                   quant_25 = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                              dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep),
                                   quant_75 = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                              dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep),
                                   quant_5 = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                             dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep),
                                   quant_95 = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                              dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep),
                                   quant_2.5 = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                               dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep),
                                   quant_97.5 = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                                dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep)
                              ),
                            emiss_sd = 
                              list(mean = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                          dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep),
                                   var = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                         dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep), 
                                   mode = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                          dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep), 
                                   median = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                            dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep), 
                                   quant_25 = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                              dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep),
                                   quant_75 = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                              dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep),
                                   quant_5 = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                             dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep),
                                   quant_95 = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                              dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep),
                                   quant_2.5 = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                               dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep),
                                   quant_97.5 = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                                dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep)
                              ),
                            emiss_betw_subj = 
                              list(mean = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                          dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep),
                                   var = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                         dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep), 
                                   mode = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                          dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep), 
                                   median = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                            dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep), 
                                   quant_25 = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                              dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep),
                                   quant_75 = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                              dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep),
                                   quant_5 = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                             dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep),
                                   quant_95 = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                              dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep),
                                   quant_2.5 = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                               dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep),
                                   quant_97.5 = rep(list(matrix(,ncol = true_m, nrow = n_conv_runs, 
                                                                dimnames = list(1:n_conv_runs, paste0("st", 1:true_m)))), n_dep)
                              )
                            
)
group_out_3st_conv2 <- group_out_3st_conv1

subj_out_3st_conv1 <- list(subj_gamma_prob  = 
                             list(mean = vector("list", n_conv_runs),
                                  median = vector("list", n_conv_runs), 
                                  var = vector("list", n_conv_runs), 
                                  quant_2.5 = vector("list", n_conv_runs), 
                                  quant_97.5 = vector("list", n_conv_runs)),
                           subj_gamma_int =
                             list(mean = vector("list", n_conv_runs), 
                                  median = vector("list", n_conv_runs),
                                  var = vector("list", n_conv_runs), 
                                  quant_2.5 = vector("list", n_conv_runs), 
                                  quant_97.5 = vector("list", n_conv_runs)), 
                           subj_emiss = 
                             list(mean = vector("list", n_conv_runs), 
                                  median = vector("list", n_conv_runs),
                                  var = vector("list", n_conv_runs), 
                                  quant_2.5 = vector("list", n_conv_runs), 
                                  quant_97.5 = vector("list", n_conv_runs))
)
subj_out_3st_conv2 <- subj_out_3st_conv1


gamma_full_conv1 <- gamma_full_conv2 <- rep(list(vector("list", max_m)), n_conv_runs) 
emiss_mean_full_conv1 <- emiss_mean_full_conv2 <- rep(list(vector("list", max_m)), n_conv_runs) 
PD_full_conv1 <- PD_full_conv2 <- rep(list(vector("list", max_m)), n_conv_runs) 

# create fucntion to obtain mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

## Convergence run 1 ####
set.seed(3571 + num_argv)
for(s in 1:n_conv_runs){
  fit_out <- mHMM(s_data = data_cont2[[s]]$obs,
                  data_distr = 'continuous',
                  gen = list(m = true_m, n_dep = n_dep),
                  start_val = c(list(gamma_start2), emiss_start2[[s]]),
                  emiss_hyp_prior = manual_prior_emiss[[s]],
                  gamma_hyp_prior = manual_prior_gamma,
                  mcmc = list(J = J, burn_in = burn_in))
  
  
  gamma_full_conv1[[s]] <- fit_out$gamma_prob_bar
  emiss_mean_full_conv1[[s]] <- fit_out$emiss_mu_bar
  PD_full_conv1[[s]] <- fit_out$PD_subj
  
  input_m   <- fit_out$input
  n_subj  <- input_m$n_subj
  burn_in <- input_m$burn_in
  J       <- input_m$J
  m       <- input_m$m
  n_vary  <- input_m$n_vary
  n_dep   <- input_m$n_dep
  
  # save group level transition probability outcomes, probability scale 
  group_out_3st_conv1$trans_prob$mean[s,] <- apply(fit_out$gamma_prob_bar[((burn_in + 1): J),], 2, mean)
  group_out_3st_conv1$trans_prob$var[s,] <- apply(fit_out$gamma_prob_bar[((burn_in + 1): J),], 2, var)
  group_out_3st_conv1$trans_prob$mode[s,] <- apply(fit_out$gamma_prob_bar[((burn_in + 1): J),], 2, Mode)
  group_out_3st_conv1$trans_prob$median[s,] <- apply(fit_out$gamma_prob_bar[((burn_in + 1): J),], 2, median)
  quantiles <- apply(fit_out$gamma_prob_bar[((burn_in + 1): J),], 2, quantile, 
                     probs = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975))
  group_out_3st_conv1$trans_prob$quant_2.5[s,] <- quantiles[1,]
  group_out_3st_conv1$trans_prob$quant_5[s,] <- quantiles[2,] 
  group_out_3st_conv1$trans_prob$quant_25[s,] <- quantiles[3,]  
  group_out_3st_conv1$trans_prob$quant_75[s,] <- quantiles[4,]
  group_out_3st_conv1$trans_prob$quant_95[s,] <-  quantiles[5,]
  group_out_3st_conv1$trans_prob$quant_97.5[s,] <- quantiles[6,]
  
  # save group level transition probability outcomes, logit scale
  group_out_3st_conv1$trans_interc$mean[s,] <- apply(fit_out$gamma_int_bar[((burn_in + 1): J),], 2, mean)
  group_out_3st_conv1$trans_interc$var[s,] <- apply(fit_out$gamma_int_bar[((burn_in + 1): J),], 2, var)
  group_out_3st_conv1$trans_interc$mode[s,] <-  apply(fit_out$gamma_int_bar[((burn_in + 1): J),], 2, Mode)
  group_out_3st_conv1$trans_interc$median[s,] <- apply(fit_out$gamma_int_bar[((burn_in + 1): J),], 2, median)
  quantiles <- apply(fit_out$gamma_int_bar[((burn_in + 1): J),], 2, quantile, 
                     probs = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975))
  group_out_3st_conv1$trans_interc$quant_2.5[s,] <- quantiles[1,]
  group_out_3st_conv1$trans_interc$quant_5[s,] <- quantiles[2,]
  group_out_3st_conv1$trans_interc$quant_25[s,] <- quantiles[3,]  
  group_out_3st_conv1$trans_interc$quant_75[s,] <- quantiles[4,]
  group_out_3st_conv1$trans_interc$quant_95[s,] <-  quantiles[5,]
  group_out_3st_conv1$trans_interc$quant_97.5[s,] <- quantiles[6,]
  
  # save group level between subject variance on transition probability outcomes, logit scale
  group_out_3st_conv1$trans_betw_subj$mean[s,] <- apply(fit_out$gamma_V_int_bar[((burn_in + 1): J),], 2, mean)
  group_out_3st_conv1$trans_betw_subj$var[s,] <- apply(fit_out$gamma_V_int_bar[((burn_in + 1): J),], 2, var)
  group_out_3st_conv1$trans_betw_subj$mode[s,] <- apply(fit_out$gamma_V_int_bar[((burn_in + 1): J),], 2, Mode)
  group_out_3st_conv1$trans_betw_subj$median[s,] <- apply(fit_out$gamma_V_int_bar[((burn_in + 1): J),], 2, median)
  quantiles <- apply(fit_out$gamma_V_int_bar[((burn_in + 1): J),], 2, quantile, 
                     probs = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975))
  group_out_3st_conv1$trans_betw_subj$quant_2.5[s,] <- quantiles[1,]
  group_out_3st_conv1$trans_betw_subj$quant_5[s,] <- quantiles[2,]
  group_out_3st_conv1$trans_betw_subj$quant_25[s,] <- quantiles[3,]  
  group_out_3st_conv1$trans_betw_subj$quant_75[s,] <- quantiles[4,]
  group_out_3st_conv1$trans_betw_subj$quant_95[s,] <-  quantiles[5,]
  group_out_3st_conv1$trans_betw_subj$quant_97.5[s,] <- quantiles[6,]
  
  
  # save group level emission distribution parameters, for each dependent variable q
  for(q in 1:n_dep){
    # variable and state dependent emission means 
    group_out_3st_conv1$emiss_mean$mean[[q]][s,] <- apply(fit_out$emiss_mu_bar[[q]][((burn_in + 1): J),], 2, mean) 
    group_out_3st_conv1$emiss_mean$var[[q]][s,] <- apply(fit_out$emiss_mu_bar[[q]][((burn_in + 1): J),], 2, var) 
    group_out_3st_conv1$emiss_mean$mode[[q]][s,] <- apply(fit_out$emiss_mu_bar[[q]][((burn_in + 1): J),], 2, Mode) 
    group_out_3st_conv1$emiss_mean$median[[q]][s,] <- apply(fit_out$emiss_mu_bar[[q]][((burn_in + 1): J),], 2, median) 
    quantiles <- apply(fit_out$emiss_mu_bar[[q]][((burn_in + 1): J),], 2, quantile, 
                       probs = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975))
    group_out_3st_conv1$emiss_mean$quant_2.5[[q]][s,] <- quantiles[1,]
    group_out_3st_conv1$emiss_mean$quant_5[[q]][s,] <- quantiles[2,]
    group_out_3st_conv1$emiss_mean$quant_25[[q]][s,] <- quantiles[3,]  
    group_out_3st_conv1$emiss_mean$quant_75[[q]][s,] <- quantiles[4,]
    group_out_3st_conv1$emiss_mean$quant_95[[q]][s,] <-  quantiles[5,]
    group_out_3st_conv1$emiss_mean$quant_97.5[[q]][s,] <- quantiles[6,]
    
    # variable and state dependent emission SD's (fixed over subjects)
    group_out_3st_conv1$emiss_sd$mean[[q]][s,] <- apply(fit_out$emiss_sd_bar[[q]][((burn_in + 1): J),], 2, mean) 
    group_out_3st_conv1$emiss_sd$var[[q]][s,] <- apply(fit_out$emiss_sd_bar[[q]][((burn_in + 1): J),], 2, var) 
    group_out_3st_conv1$emiss_sd$mode[[q]][s,] <- apply(fit_out$emiss_sd_bar[[q]][((burn_in + 1): J),], 2, Mode) 
    group_out_3st_conv1$emiss_sd$median[[q]][s,] <- apply(fit_out$emiss_sd_bar[[q]][((burn_in + 1): J),], 2, median) 
    quantiles <- apply(fit_out$emiss_sd_bar[[q]][((burn_in + 1): J),], 2, quantile, 
                       probs = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975))
    group_out_3st_conv1$emiss_sd$quant_2.5[[q]][s,] <- quantiles[1,]
    group_out_3st_conv1$emiss_sd$quant_5[[q]][s,] <- quantiles[2,]
    group_out_3st_conv1$emiss_sd$quant_25[[q]][s,] <- quantiles[3,]  
    group_out_3st_conv1$emiss_sd$quant_75[[q]][s,] <- quantiles[4,]
    group_out_3st_conv1$emiss_sd$quant_95[[q]][s,] <-  quantiles[5,]
    group_out_3st_conv1$emiss_sd$quant_97.5[[q]][s,] <- quantiles[6,]
    
    # variable and state dependent emission between subject variance on means
    group_out_3st_conv1$emiss_betw_subj$mean[[q]][s,] <- apply(fit_out$emiss_varmu_bar[[q]][((burn_in + 1): J),], 2, mean) 
    group_out_3st_conv1$emiss_betw_subj$var[[q]][s,] <- apply(fit_out$emiss_varmu_bar[[q]][((burn_in + 1): J),], 2, var) 
    group_out_3st_conv1$emiss_betw_subj$mode[[q]][s,] <- apply(fit_out$emiss_varmu_bar[[q]][((burn_in + 1): J),], 2, Mode) 
    group_out_3st_conv1$emiss_betw_subj$median[[q]][s,] <- apply(fit_out$emiss_varmu_bar[[q]][((burn_in + 1): J),], 2, median) 
    quantiles <- apply(fit_out$emiss_varmu_bar[[q]][((burn_in + 1): J),], 2, quantile, 
                       probs = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975))
    group_out_3st_conv1$emiss_betw_subj$quant_2.5[[q]][s,] <- quantiles[1,]
    group_out_3st_conv1$emiss_betw_subj$quant_5[[q]][s,] <- quantiles[2,]
    group_out_3st_conv1$emiss_betw_subj$quant_25[[q]][s,] <- quantiles[3,]  
    group_out_3st_conv1$emiss_betw_subj$quant_75[[q]][s,] <- quantiles[4,]
    group_out_3st_conv1$emiss_betw_subj$quant_95[[q]][s,] <-  quantiles[5,]
    group_out_3st_conv1$emiss_betw_subj$quant_97.5[[q]][s,] <- quantiles[6,]
  }
  # save selection of subject level parameters
  subj_out_3st_conv1$subj_gamma_prob$mean[[s]] <- t(sapply(sapply(fit_out$PD_subj,"[[", 1, simplify = FALSE), apply, 2, mean))
  subj_out_3st_conv1$subj_gamma_prob$median[[s]] <- t(sapply(sapply(fit_out$PD_subj,"[[", 1, simplify = FALSE), apply, 2, median))
  subj_out_3st_conv1$subj_gamma_prob$var[[s]] <- t(sapply(sapply(fit_out$PD_subj,"[[", 1, simplify = FALSE), apply, 2, var))
  subj_out_3st_conv1$subj_gamma_prob$quant_2.5[[s]] <- t(sapply(sapply(fit_out$PD_subj,"[[", 1, simplify = FALSE), apply, 2, quantile, prob = c(0.025)))
  subj_out_3st_conv1$subj_gamma_prob$quant_97.5[[s]] <- t(sapply(sapply(fit_out$PD_subj,"[[", 1, simplify = FALSE), apply, 2, quantile, prob = c(0.975)))
  
  subj_out_3st_conv1$subj_gamma_int$mean[[s]] <- t(sapply(fit_out$gamma_int_subj, apply, 2, mean))
  subj_out_3st_conv1$subj_gamma_int$median[[s]] <- t(sapply(fit_out$gamma_int_subj, apply, 2, median))
  subj_out_3st_conv1$subj_gamma_int$var[[s]] <- t(sapply(fit_out$gamma_int_subj, apply, 2, var))
  subj_out_3st_conv1$subj_gamma_int$quant_2.5[[s]] <- t(sapply(fit_out$gamma_int_subj, apply, 2, quantile, prob = c(0.025)))
  subj_out_3st_conv1$subj_gamma_int$quant_97.5[[s]] <- t(sapply(fit_out$gamma_int_subj, apply, 2, quantile, prob = c(0.975)))
  
  subj_out_3st_conv1$subj_emiss$mean[[s]]  <- t(sapply(sapply(fit_out$PD_subj,"[[", 3, simplify = FALSE), apply, 2, mean))[,1:(true_m*n_dep)]
  subj_out_3st_conv1$subj_emiss$median[[s]]  <- t(sapply(sapply(fit_out$PD_subj,"[[", 3, simplify = FALSE), apply, 2, median))[,1:(true_m*n_dep)]
  subj_out_3st_conv1$subj_emiss$var[[s]]  <-t(sapply(sapply(fit_out$PD_subj,"[[", 3, simplify = FALSE), apply, 2, var))[,1:(true_m*n_dep)]
  subj_out_3st_conv1$subj_emiss$quant_2.5[[s]]  <-t(sapply(sapply(fit_out$PD_subj,"[[", 3, simplify = FALSE), apply, 2, quantile, prob = c(0.025)))[,1:(true_m*n_dep)]
  subj_out_3st_conv1$subj_emiss$quant_97.5[[s]]  <-t(sapply(sapply(fit_out$PD_subj,"[[", 3, simplify = FALSE), apply, 2, quantile, prob = c(0.975)))[,1:(true_m*n_dep)]
}

## Convergence run 2 ####

set.seed(9834 + num_argv)
for(s in 1:n_conv_runs){
  fit_out <- mHMM(s_data = data_cont2[[s]]$obs,
                  data_distr = 'continuous',
                  gen = list(m = true_m, n_dep = n_dep),
                  start_val = c(list(gamma_start3), emiss_start3[[s]]),
                  emiss_hyp_prior = manual_prior_emiss[[s]],
                  gamma_hyp_prior = manual_prior_gamma,
                  mcmc = list(J = J, burn_in = burn_in))
  
  
  gamma_full_conv2[[s]] <- fit_out$gamma_prob_bar
  emiss_mean_full_conv2[[s]] <- fit_out$emiss_mu_bar
  PD_full_conv2[[s]] <- fit_out$PD_subj
  
  input_m   <- fit_out$input
  n_subj  <- input_m$n_subj
  burn_in <- input_m$burn_in
  J       <- input_m$J
  m       <- input_m$m
  n_vary  <- input_m$n_vary
  n_dep   <- input_m$n_dep
  
  # save group level transition probability outcomes, probability scale 
  group_out_3st_conv2$trans_prob$mean[s,] <- apply(fit_out$gamma_prob_bar[((burn_in + 1): J),], 2, mean)
  group_out_3st_conv2$trans_prob$var[s,] <- apply(fit_out$gamma_prob_bar[((burn_in + 1): J),], 2, var)
  group_out_3st_conv2$trans_prob$mode[s,] <- apply(fit_out$gamma_prob_bar[((burn_in + 1): J),], 2, Mode)
  group_out_3st_conv2$trans_prob$median[s,] <- apply(fit_out$gamma_prob_bar[((burn_in + 1): J),], 2, median)
  quantiles <- apply(fit_out$gamma_prob_bar[((burn_in + 1): J),], 2, quantile, 
                     probs = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975))
  group_out_3st_conv2$trans_prob$quant_2.5[s,] <- quantiles[1,]
  group_out_3st_conv2$trans_prob$quant_5[s,] <- quantiles[2,] 
  group_out_3st_conv2$trans_prob$quant_25[s,] <- quantiles[3,]  
  group_out_3st_conv2$trans_prob$quant_75[s,] <- quantiles[4,]
  group_out_3st_conv2$trans_prob$quant_95[s,] <-  quantiles[5,]
  group_out_3st_conv2$trans_prob$quant_97.5[s,] <- quantiles[6,]
  
  # save group level transition probability outcomes, logit scale
  group_out_3st_conv2$trans_interc$mean[s,] <- apply(fit_out$gamma_int_bar[((burn_in + 1): J),], 2, mean)
  group_out_3st_conv2$trans_interc$var[s,] <- apply(fit_out$gamma_int_bar[((burn_in + 1): J),], 2, var)
  group_out_3st_conv2$trans_interc$mode[s,] <-  apply(fit_out$gamma_int_bar[((burn_in + 1): J),], 2, Mode)
  group_out_3st_conv2$trans_interc$median[s,] <- apply(fit_out$gamma_int_bar[((burn_in + 1): J),], 2, median)
  quantiles <- apply(fit_out$gamma_int_bar[((burn_in + 1): J),], 2, quantile, 
                     probs = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975))
  group_out_3st_conv2$trans_interc$quant_2.5[s,] <- quantiles[1,]
  group_out_3st_conv2$trans_interc$quant_5[s,] <- quantiles[2,]
  group_out_3st_conv2$trans_interc$quant_25[s,] <- quantiles[3,]  
  group_out_3st_conv2$trans_interc$quant_75[s,] <- quantiles[4,]
  group_out_3st_conv2$trans_interc$quant_95[s,] <-  quantiles[5,]
  group_out_3st_conv2$trans_interc$quant_97.5[s,] <- quantiles[6,]
  
  # save group level between subject variance on transition probability outcomes, logit scale
  group_out_3st_conv2$trans_betw_subj$mean[s,] <- apply(fit_out$gamma_V_int_bar[((burn_in + 1): J),], 2, mean)
  group_out_3st_conv2$trans_betw_subj$var[s,] <- apply(fit_out$gamma_V_int_bar[((burn_in + 1): J),], 2, var)
  group_out_3st_conv2$trans_betw_subj$mode[s,] <- apply(fit_out$gamma_V_int_bar[((burn_in + 1): J),], 2, Mode)
  group_out_3st_conv2$trans_betw_subj$median[s,] <- apply(fit_out$gamma_V_int_bar[((burn_in + 1): J),], 2, median)
  quantiles <- apply(fit_out$gamma_V_int_bar[((burn_in + 1): J),], 2, quantile, 
                     probs = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975))
  group_out_3st_conv2$trans_betw_subj$quant_2.5[s,] <- quantiles[1,]
  group_out_3st_conv2$trans_betw_subj$quant_5[s,] <- quantiles[2,]
  group_out_3st_conv2$trans_betw_subj$quant_25[s,] <- quantiles[3,]  
  group_out_3st_conv2$trans_betw_subj$quant_75[s,] <- quantiles[4,]
  group_out_3st_conv2$trans_betw_subj$quant_95[s,] <-  quantiles[5,]
  group_out_3st_conv2$trans_betw_subj$quant_97.5[s,] <- quantiles[6,]
  
  
  # save group level emission distribution parameters, for each dependent variable q
  for(q in 1:n_dep){
    # variable and state dependent emission means 
    group_out_3st_conv2$emiss_mean$mean[[q]][s,] <- apply(fit_out$emiss_mu_bar[[q]][((burn_in + 1): J),], 2, mean) 
    group_out_3st_conv2$emiss_mean$var[[q]][s,] <- apply(fit_out$emiss_mu_bar[[q]][((burn_in + 1): J),], 2, var) 
    group_out_3st_conv2$emiss_mean$mode[[q]][s,] <- apply(fit_out$emiss_mu_bar[[q]][((burn_in + 1): J),], 2, Mode) 
    group_out_3st_conv2$emiss_mean$median[[q]][s,] <- apply(fit_out$emiss_mu_bar[[q]][((burn_in + 1): J),], 2, median) 
    quantiles <- apply(fit_out$emiss_mu_bar[[q]][((burn_in + 1): J),], 2, quantile, 
                       probs = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975))
    group_out_3st_conv2$emiss_mean$quant_2.5[[q]][s,] <- quantiles[1,]
    group_out_3st_conv2$emiss_mean$quant_5[[q]][s,] <- quantiles[2,]
    group_out_3st_conv2$emiss_mean$quant_25[[q]][s,] <- quantiles[3,]  
    group_out_3st_conv2$emiss_mean$quant_75[[q]][s,] <- quantiles[4,]
    group_out_3st_conv2$emiss_mean$quant_95[[q]][s,] <-  quantiles[5,]
    group_out_3st_conv2$emiss_mean$quant_97.5[[q]][s,] <- quantiles[6,]
    
    # variable and state dependent emission SD's (fixed over subjects)
    group_out_3st_conv2$emiss_sd$mean[[q]][s,] <- apply(fit_out$emiss_sd_bar[[q]][((burn_in + 1): J),], 2, mean) 
    group_out_3st_conv2$emiss_sd$var[[q]][s,] <- apply(fit_out$emiss_sd_bar[[q]][((burn_in + 1): J),], 2, var) 
    group_out_3st_conv2$emiss_sd$mode[[q]][s,] <- apply(fit_out$emiss_sd_bar[[q]][((burn_in + 1): J),], 2, Mode) 
    group_out_3st_conv2$emiss_sd$median[[q]][s,] <- apply(fit_out$emiss_sd_bar[[q]][((burn_in + 1): J),], 2, median) 
    quantiles <- apply(fit_out$emiss_sd_bar[[q]][((burn_in + 1): J),], 2, quantile, 
                       probs = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975))
    group_out_3st_conv2$emiss_sd$quant_2.5[[q]][s,] <- quantiles[1,]
    group_out_3st_conv2$emiss_sd$quant_5[[q]][s,] <- quantiles[2,]
    group_out_3st_conv2$emiss_sd$quant_25[[q]][s,] <- quantiles[3,]  
    group_out_3st_conv2$emiss_sd$quant_75[[q]][s,] <- quantiles[4,]
    group_out_3st_conv2$emiss_sd$quant_95[[q]][s,] <-  quantiles[5,]
    group_out_3st_conv2$emiss_sd$quant_97.5[[q]][s,] <- quantiles[6,]
    
    # variable and state dependent emission between subject variance on means
    group_out_3st_conv2$emiss_betw_subj$mean[[q]][s,] <- apply(fit_out$emiss_varmu_bar[[q]][((burn_in + 1): J),], 2, mean) 
    group_out_3st_conv2$emiss_betw_subj$var[[q]][s,] <- apply(fit_out$emiss_varmu_bar[[q]][((burn_in + 1): J),], 2, var) 
    group_out_3st_conv2$emiss_betw_subj$mode[[q]][s,] <- apply(fit_out$emiss_varmu_bar[[q]][((burn_in + 1): J),], 2, Mode) 
    group_out_3st_conv2$emiss_betw_subj$median[[q]][s,] <- apply(fit_out$emiss_varmu_bar[[q]][((burn_in + 1): J),], 2, median) 
    quantiles <- apply(fit_out$emiss_varmu_bar[[q]][((burn_in + 1): J),], 2, quantile, 
                       probs = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975))
    group_out_3st_conv2$emiss_betw_subj$quant_2.5[[q]][s,] <- quantiles[1,]
    group_out_3st_conv2$emiss_betw_subj$quant_5[[q]][s,] <- quantiles[2,]
    group_out_3st_conv2$emiss_betw_subj$quant_25[[q]][s,] <- quantiles[3,]  
    group_out_3st_conv2$emiss_betw_subj$quant_75[[q]][s,] <- quantiles[4,]
    group_out_3st_conv2$emiss_betw_subj$quant_95[[q]][s,] <-  quantiles[5,]
    group_out_3st_conv2$emiss_betw_subj$quant_97.5[[q]][s,] <- quantiles[6,]
  }
  # save selection of subject level parameters
  subj_out_3st_conv2$subj_gamma_prob$mean[[s]] <- t(sapply(sapply(fit_out$PD_subj,"[[", 1, simplify = FALSE), apply, 2, mean))
  subj_out_3st_conv2$subj_gamma_prob$median[[s]] <- t(sapply(sapply(fit_out$PD_subj,"[[", 1, simplify = FALSE), apply, 2, median))
  subj_out_3st_conv2$subj_gamma_prob$var[[s]] <- t(sapply(sapply(fit_out$PD_subj,"[[", 1, simplify = FALSE), apply, 2, var))
  subj_out_3st_conv2$subj_gamma_prob$quant_2.5[[s]] <- t(sapply(sapply(fit_out$PD_subj,"[[", 1, simplify = FALSE), apply, 2, quantile, prob = c(0.025)))
  subj_out_3st_conv2$subj_gamma_prob$quant_97.5[[s]] <- t(sapply(sapply(fit_out$PD_subj,"[[", 1, simplify = FALSE), apply, 2, quantile, prob = c(0.975)))
  
  subj_out_3st_conv2$subj_gamma_int$mean[[s]] <- t(sapply(fit_out$gamma_int_subj, apply, 2, mean))
  subj_out_3st_conv2$subj_gamma_int$median[[s]] <- t(sapply(fit_out$gamma_int_subj, apply, 2, median))
  subj_out_3st_conv2$subj_gamma_int$var[[s]] <- t(sapply(fit_out$gamma_int_subj, apply, 2, var))
  subj_out_3st_conv2$subj_gamma_int$quant_2.5[[s]] <- t(sapply(fit_out$gamma_int_subj, apply, 2, quantile, prob = c(0.025)))
  subj_out_3st_conv2$subj_gamma_int$quant_97.5[[s]] <- t(sapply(fit_out$gamma_int_subj, apply, 2, quantile, prob = c(0.975)))
  
  subj_out_3st_conv2$subj_emiss$mean[[s]]  <- t(sapply(sapply(fit_out$PD_subj,"[[", 3, simplify = FALSE), apply, 2, mean))[,1:(true_m*n_dep)]
  subj_out_3st_conv2$subj_emiss$median[[s]]  <- t(sapply(sapply(fit_out$PD_subj,"[[", 3, simplify = FALSE), apply, 2, median))[,1:(true_m*n_dep)]
  subj_out_3st_conv2$subj_emiss$var[[s]]  <-t(sapply(sapply(fit_out$PD_subj,"[[", 3, simplify = FALSE), apply, 2, var))[,1:(true_m*n_dep)]
  subj_out_3st_conv2$subj_emiss$quant_2.5[[s]]  <-t(sapply(sapply(fit_out$PD_subj,"[[", 3, simplify = FALSE), apply, 2, quantile, prob = c(0.025)))[,1:(true_m*n_dep)]
  subj_out_3st_conv2$subj_emiss$quant_97.5[[s]]  <-t(sapply(sapply(fit_out$PD_subj,"[[", 3, simplify = FALSE), apply, 2, quantile, prob = c(0.975)))[,1:(true_m*n_dep)]
}

## Save results ####

out_conv1 = list(gamma_full_conv1 = gamma_full_conv1,
                 emiss_mean_full_conv1 = emiss_mean_full_conv1,
                 PD_full_conv1 = PD_full_conv1,
                 group_out_3st_conv1 = group_out_3st_conv1, 
                 subj_out_3st_conv1 = subj_out_3st_conv1, 
                 input = input, 
                 subj_par_conv1 = subj_par2, 
                 means_true_conv1 = means_true2)

out_conv2 = list(gamma_full_conv2 = gamma_full_conv2,
                 emiss_mean_full_conv2 = emiss_mean_full_conv2,
                 PD_full_conv2 = PD_full_conv2,
                 group_out_3st_conv2 = group_out_3st_conv2, 
                 subj_out_3st_conv2 = subj_out_3st_conv2, 
                 input = input, 
                 subj_par_conv2 = subj_par2, 
                 means_true_conv2 = means_true2)

saveRDS(out_conv1, file = paste0("out_V4/", out_file, "_conv1_1.rds"))
saveRDS(out_conv2, file = paste0("out_V4/", out_file, "_conv2_1.rds"))



