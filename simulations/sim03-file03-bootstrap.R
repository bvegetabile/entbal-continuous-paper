source('simulation-functions/load_all_funcs_and_libs.R')

# ------------------------------------------------------------------------
# Simulation parameters ------
source('simulations/simulation-parameters-01.R')

# ------------------------------------------------------------------------

# #--- Visualizing one span solution
n_moments <- 3
set.seed(20200109)
n_obs <- 200
t1 <- system.time(testdat <- simrun_lite(X=33, n_obs = n_obs, MX1 = MX1, MX2 = MX2, MX3 = MX3, nm=n_moments))
print(t1)

# # -- Bootstrap Simulation
NB <- 100
NS <- n_sims <- 1000
resmat <- matrix(NA, nrow = NS, ncol = 100)
bscovers <- matrix(NA, nrow = NS, ncol = 100)
simests <- matrix(NA, nrow = NS, ncol = 100)

haseffect <- TRUE

for(b in 1:NS){
  simdat <- simrun_lite(b+1, n_obs = 200, 
                        MX1 = MX1, MX2 = MX2, MX3 = MX3,
                        amin = amin, amax = amax, 
                        nm=n_moments, a_effect = haseffect,
                        lowspan = 0.25)  
  simests[b,] <- simdat[[1]]$ests_loess
  options(cores = 20)
  cl <- makeCluster(20)
  
  registerDoParallel(cl)
  bstime <- system.time(bsres <- foreach(
    i = 1:NB,
    .combine = rbind,
    .verbose = F) %dopar% bsest_smalln(i, simdat$dat, A_test, nm = n_moments, sval = simdat$eb1$spanres))
  stopCluster(cl)
  print(bstime)
  cat(paste('Sim:', b, '\n'))
  message(paste('Sim:',b, ', time diff', round(bstime[3],3))) 
  lower <- simdat[[1]]$ests_loess - 1.96 * apply(bsres, 2, sd)
  upper <- simdat[[1]]$ests_loess + 1.96 * apply(bsres, 2, sd)
  
  covers <- (upper > truth) & (lower < truth)
  bscovers[b, ] <- covers
  resmat[b,] <- apply(bsres, 2, sd)
}

outro <- list('bs-cov' = bscovers, 'bs-se' = resmat, 'ests' = simests, 'truth' = truth)

saveRDS(outro, file = 'output/sim03-bootstrap-results.rds')  
