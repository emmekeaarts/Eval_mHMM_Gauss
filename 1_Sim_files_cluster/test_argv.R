argv = commandArgs(trailingOnly=TRUE)
print(argv)
num_argv = as.numeric(argv)
print(num_argv)

settings <- read.csv("sim_scenarios.csv")

KL_div  <- settings[num_argv, 4]
n_t     <- settings[num_argv, 2]
n       <- settings[num_argv, 1]
n_dep   <- settings[num_argv, 3]

true_m  <- 3
max_m   <- 4
n_sim   <- 20

J       <- 20
burn_in <- 10

out_file <- paste0("out_n", n, "_nt", n_t, "_ndep", n_dep, "_", true_m, "st_KLD", KL_div)
print(out_file)