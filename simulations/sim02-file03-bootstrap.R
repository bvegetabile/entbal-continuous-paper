source('simulation-functions/load_all_funcs_and_libs.R')

# ------------------------------------------------------------------------
# Simulation parameters ------
source('simulations/simulation-parameters-02.R')

# ------------------------------------------------------------------------

# #--- Visualizing one span solution
set.seed(20200109)
n_obs <- 1000
t1 <- system.time(testdat <- simrun(X=33, n_obs = n_obs, MX1 = MX1, MX2 = MX2, MX3 = MX3, lite_version = T, a_effect = F))
print(t1)

# # -- Bootstrap Simulation
NB <- 100
NS <- n_sims <- 1000
resmat <- matrix(NA, nrow = NS, ncol = 100)
bscovers <- matrix(NA, nrow = NS, ncol = 100)
simests <- matrix(NA, nrow = NS, ncol = 100)
n_moments <- 3
haseffect <- FALSE
for(b in 1:NS){
  simdat <- simrun(b+1, n_obs = 1000, 
                   MX1 = MX1, MX2 = MX2, MX3 = MX3,
                   amin = amin, amax = amax, 
                   lite_version = T, a_effect = haseffect)  
  simests[b,] <- simdat[[n_moments]]$ests_loess
  options(cores = 20)
  cl <- makeCluster(20)
  
  registerDoParallel(cl)
  bstime <- system.time(bsres <- foreach(
    i = 1:NB,
    .combine = rbind,
    .verbose = F) %dopar% bsest(i, simdat$dat, A_test, nm = n_moments))
  stopCluster(cl)
  print(bstime)
  cat(paste('Sim:', b, '\n'))
  message(paste('Sim:',b, ', time diff', round(bstime[3],3))) 
  lower <- simdat[[n_moments]]$ests_loess - 1.96 * apply(bsres, 2, sd)
  upper <- simdat[[n_moments]]$ests_loess + 1.96 * apply(bsres, 2, sd)
  
  covers <- (upper > truth) & (lower < truth)
  bscovers[b, ] <- covers
  resmat[b,] <- apply(bsres, 2, sd)
}

outro <- list('bs-cov' = bscovers, 'bs-se' = resmat, 'ests' = simests, 'truth' = truth)

if(haseffect){
  saveRDS(outro, file = 'output/sim01-bootstrap-results.rds')  
} else {
  saveRDS(outro, file = 'output/sim02-bootstrap-results.rds')  
}
