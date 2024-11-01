n_dep <- 4
m <- 3
KL_div <- 15

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
    
    return(sum(v_error))
    
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

means <- makeGGM_KLD(K = m, p = n_dep, target_KLD = KL_div, sig = 1, init_range = c(-5, 5))
means2 <- matrix(unlist(means$mu), ncol = m, nrow = n_dep)
rownames(means2) <- paste0("dep_", 1:n_dep)
colnames(means2) <- paste0("state_", 1:m)

emiss_distr <- rep(list(matrix(, ncol = 2, nrow = m)), n_dep)
for(q in 1:n_dep){
  emiss_distr[[q]][,1] <- means2[q,]
  emiss_distr[[q]][,2] <- 1
}

x <- seq(-5, 5, 0.1)

library(ggplot2)

dens <- data.frame(dep = factor(rep(1:n_dep, each = (length(x) * m))),
                   state = factor(rep(rep(1:m, each = length(x), n_dep))), 
                   x = rep(x, n_dep*m),
                   dens = NA)

for(q in 1:n_dep){
  for (i in 1:m){
    dens$dens[dens$dep == q & dens$state == i] <- dnorm(x, mean = emiss_distr[[q]][i,1], sd = 1)
  }
}

ggplot(data = dens, mapping = aes(x = x, y = dens, grouping = state, col = state)) +
  geom_line() +
  facet_wrap(vars(dep), ncol = 1)


emiss_distr <- list(matrix(c( -2.1, 1,
                              -1.1, 1), nrow = m, byrow = TRUE),
                    matrix(c(3.8, 1,
                             2.1, 1), nrow = m, byrow = TRUE))


plot_ndep2_KL2 <- ggplot(data = dens, mapping = aes(x = x, y = dens, grouping = state, col = state)) +
  geom_line() +
  facet_wrap(vars(dep), ncol = 1) +
  theme(legend.position = "none")

plot_ndep4_KL2 <- ggplot(data = dens, mapping = aes(x = x, y = dens, grouping = state, col = state)) +
  geom_line() +
  facet_wrap(vars(dep), ncol = 1) +
  theme(legend.position = "none")

plot_ndep8_KL2 <- ggplot(data = dens, mapping = aes(x = x, y = dens, grouping = state, col = state)) +
  geom_line() +
  facet_wrap(vars(dep), ncol = 2)

plot_ndep2_KL35 <- ggplot(data = dens, mapping = aes(x = x, y = dens, grouping = state, col = state)) +
  geom_line() +
  facet_wrap(vars(dep), ncol = 1) +
  theme(legend.position = "none")

plot_ndep4_KL35 <- ggplot(data = dens, mapping = aes(x = x, y = dens, grouping = state, col = state)) +
  geom_line() +
  facet_wrap(vars(dep), ncol = 1) +
  theme(legend.position = "none")

plot_ndep8_KL35 <- ggplot(data = dens, mapping = aes(x = x, y = dens, grouping = state, col = state)) +
  geom_line() +
  facet_wrap(vars(dep), ncol = 2)

plot_ndep2_KL5 <- ggplot(data = dens, mapping = aes(x = x, y = dens, grouping = state, col = state)) +
  geom_line() +
  facet_wrap(vars(dep), ncol = 1) +
  theme(legend.position = "none")

plot_ndep4_KL5 <- ggplot(data = dens, mapping = aes(x = x, y = dens, grouping = state, col = state)) +
  geom_line() +
  facet_wrap(vars(dep), ncol = 1) +
  theme(legend.position = "none")

plot_ndep8_KL5 <- ggplot(data = dens, mapping = aes(x = x, y = dens, grouping = state, col = state)) +
  geom_line() +
  facet_wrap(vars(dep), ncol = 2)
plot_ndep8_KL5


library(cowplot)

plot_grid(plot_ndep2_KL2, plot_ndep4_KL2, plot_ndep8_KL2, nrow = 1, rel_widths = c(1, 1, 2))

plot_grid(plot_ndep2_KL35, plot_ndep4_KL35, plot_ndep8_KL35, nrow = 1, rel_widths = c(1, 1, 2))

plot_grid(plot_ndep2_KL5, plot_ndep4_KL5, plot_ndep8_KL5, nrow = 1, rel_widths = c(1, 1, 2))
