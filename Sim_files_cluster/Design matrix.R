# create design matrix, specify: 
# KL_div, n_t, n, n_dep, 


KLD <- c(2, 3.5, 5)
n <- rep(c(30, 60, 120, 240), times = c(5, 5, 4, 3))
n_t <- c(rep(c(50,100, 200, 400, 800),2), 50, 100, 200, 400, 50, 100, 200)

sim_scenarios <- rbind(data.frame(n = n, 
                            n_t = n_t, 
                            ndep = 2, 
                            KLD = 2),
                 data.frame(n = n, 
                            n_t = n_t, 
                            ndep = 4, 
                            KLD = 2),
                 data.frame(n = n, 
                            n_t = n_t, 
                            ndep = 8, 
                            KLD = 2),
                 data.frame(n = n, 
                            n_t = n_t, 
                            ndep = 2, 
                            KLD = 3.5),
                 data.frame(n = n, 
                            n_t = n_t, 
                            ndep = 4, 
                            KLD = 3.5),
                 data.frame(n = n, 
                            n_t = n_t, 
                            ndep = 8, 
                            KLD = 3.5),
                 data.frame(n = n, 
                            n_t = n_t, 
                            ndep = 2, 
                            KLD = 5),
                 data.frame(n = n, 
                            n_t = n_t, 
                            ndep = 4, 
                            KLD = 5),
                 data.frame(n = n, 
                            n_t = n_t, 
                            ndep = 8, 
                            KLD = 5)
                 )


sim_scenarios <- sim_scenarios[order(sim_scenarios$n, sim_scenarios$n_t),]
table(sim_scenarios$n)

write.csv(sim_scenarios, file = "sim_scenarios.csv", row.names = FALSE)



# This file needs to be iterated over 3 times, once for true states 1, true states 2, and true states 3. 
# use 32 cores in node, need 3 nodes to together cover 30 and 60 pp (each node set to iterating over 30 scenarios), 
# and the remaining 4 of 120 (iterate over 32 cores) will be added to 240 (totalling to 31 cores)

# note: the R scripts differs over the iterations of the true number of states, hence, this parameter is not included here


########## KLD2 ##########
# Using different KLD values, more close to what we found in real-world settings
##########################

KLD2 <- c(4, 8, 12)
n <- rep(c(30, 60, 120, 240), times = c(5, 5, 4, 3))
n_t <- c(rep(c(50,100, 200, 400, 800),2), 50, 100, 200, 400, 50, 100, 200)

sim_scenarios_V2 <- rbind(data.frame(n = n, 
                                  n_t = n_t, 
                                  ndep = 2, 
                                  KLD = 4),
                       data.frame(n = n, 
                                  n_t = n_t, 
                                  ndep = 4, 
                                  KLD = 4),
                       data.frame(n = n, 
                                  n_t = n_t, 
                                  ndep = 8, 
                                  KLD = 4),
                       data.frame(n = n, 
                                  n_t = n_t, 
                                  ndep = 2, 
                                  KLD = 8),
                       data.frame(n = n, 
                                  n_t = n_t, 
                                  ndep = 4, 
                                  KLD = 8),
                       data.frame(n = n, 
                                  n_t = n_t, 
                                  ndep = 8, 
                                  KLD = 8),
                       data.frame(n = n, 
                                  n_t = n_t, 
                                  ndep = 2, 
                                  KLD = 12),
                       data.frame(n = n, 
                                  n_t = n_t, 
                                  ndep = 4, 
                                  KLD = 12),
                       data.frame(n = n, 
                                  n_t = n_t, 
                                  ndep = 8, 
                                  KLD = 12)
)


sim_scenarios_V2 <- sim_scenarios_V2[order(sim_scenarios_V2$n, sim_scenarios_V2$n_t),]
table(sim_scenarios_V2$n)

write.csv(sim_scenarios_V2, file = "sim_scenarios_V2.csv", row.names = FALSE)


########## KLD3 ##########
# final version on KLD settings and added group size number, more close to what we found in real-world settings
##########################

KLD2 <- c(5, 10, 15)
n <- rep(c(15, 30, 60, 120, 240), times = c(5, 5, 5, 4, 3))
n_t <- c(rep(c(50,100, 200, 400, 800), 3), 50, 100, 200, 400, 50, 100, 200)

sim_scenarios_V3 <- rbind(data.frame(n = n, 
                                     n_t = n_t, 
                                     ndep = 2, 
                                     KLD = 5),
                          data.frame(n = n, 
                                     n_t = n_t, 
                                     ndep = 4, 
                                     KLD = 5),
                          data.frame(n = n, 
                                     n_t = n_t, 
                                     ndep = 8, 
                                     KLD = 5),
                          data.frame(n = n, 
                                     n_t = n_t, 
                                     ndep = 2, 
                                     KLD = 10),
                          data.frame(n = n, 
                                     n_t = n_t, 
                                     ndep = 4, 
                                     KLD = 10),
                          data.frame(n = n, 
                                     n_t = n_t, 
                                     ndep = 8, 
                                     KLD = 10),
                          data.frame(n = n, 
                                     n_t = n_t, 
                                     ndep = 2, 
                                     KLD = 15),
                          data.frame(n = n, 
                                     n_t = n_t, 
                                     ndep = 4, 
                                     KLD = 15),
                          data.frame(n = n, 
                                     n_t = n_t, 
                                     ndep = 8, 
                                     KLD = 15)
)


sim_scenarios_V3 <- sim_scenarios_V3[order(sim_scenarios_V3$n, sim_scenarios_V3$n_t, sim_scenarios_V3$ndep),]
table(sim_scenarios_V3$n)

sim_scenarios <- sim_scenarios_V3
rownames(sim_scenarios) <- 1:dim(sim_scenarios_V3)[1]

sim_scenarios3 <- rbind(sim_scenarios[1:27,],
                       sim_scenarios[46:48,],
                       sim_scenarios[28:45,],
                       sim_scenarios[49:81,],
                       sim_scenarios[91:99,],
                       sim_scenarios[82:90,],
                       sim_scenarios[100:117,],
                       sim_scenarios[136:138,],
                       sim_scenarios[118:126,],
                       sim_scenarios[139:153,],
                       sim_scenarios[127:135,],
                       sim_scenarios[154:162,],
                       sim_scenarios[172:177,],
                       sim_scenarios[163:171,],
                       sim_scenarios[178:198,]
                       )
rownames(sim_scenarios3) <- 1:dim(sim_scenarios3)[1]

write.csv(sim_scenarios3, file = "sim_scenarios_V3.csv", row.names = FALSE)

####### And again different values #########
n <- rep(c(15, 30, 60, 120), times = c(5, 5, 5, 4))
n_t <- c(rep(c(50,100, 200, 400, 800), 3), 50, 100, 200, 400)

sim_scenarios_V4 <- rbind(data.frame(n = n, 
                                     n_t = n_t, 
                                     ndep = 2, 
                                     KLD = 3),
                          data.frame(n = n, 
                                     n_t = n_t, 
                                     ndep = 4, 
                                     KLD = 3),
                          data.frame(n = n, 
                                     n_t = n_t, 
                                     ndep = 8, 
                                     KLD = 3),
                          data.frame(n = n, 
                                     n_t = n_t, 
                                     ndep = 2, 
                                     KLD = 5),
                          data.frame(n = n, 
                                     n_t = n_t, 
                                     ndep = 4, 
                                     KLD = 5),
                          data.frame(n = n, 
                                     n_t = n_t, 
                                     ndep = 8, 
                                     KLD = 5),
                          data.frame(n = n, 
                                     n_t = n_t, 
                                     ndep = 2, 
                                     KLD = 7),
                          data.frame(n = n, 
                                     n_t = n_t, 
                                     ndep = 4, 
                                     KLD = 7),
                          data.frame(n = n, 
                                     n_t = n_t, 
                                     ndep = 8, 
                                     KLD = 7)
)


sim_scenarios_V4 <- sim_scenarios_V4[order(sim_scenarios_V4$n, sim_scenarios_V4$n_t, sim_scenarios_V4$ndep),]
table(sim_scenarios_V4$n)

sim_scenarios <- sim_scenarios_V4
rownames(sim_scenarios) <- 1:dim(sim_scenarios_V4)[1]

sim_scenarios4 <- rbind(sim_scenarios[1:27,],
                        sim_scenarios[46:48,],
                        sim_scenarios[28:45,],
                        sim_scenarios[49:81,],
                        sim_scenarios[91:99,],
                        sim_scenarios[82:90,],
                        sim_scenarios[100:117,],
                        sim_scenarios[136:138,],
                        sim_scenarios[118:126,],
                        sim_scenarios[139:153,],
                        sim_scenarios[127:135,],
                        sim_scenarios[154:171,]
)
rownames(sim_scenarios4) <- 1:dim(sim_scenarios4)[1]

write.csv(sim_scenarios4, file = "sim_scenarios_V4.csv", row.names = FALSE)
