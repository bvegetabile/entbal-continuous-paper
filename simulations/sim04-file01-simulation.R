source('simulation-functions/load_all_funcs_and_libs.R')

# ------------------------------------------------------------------------
# Simulation parameters ------
source('simulations/simulation-parameters-01.R')

# ------------------------------------------------------------------------

# #--- Visualizing one span solution
set.seed(20200109)
n_obs <- 1000
t1 <- system.time(testdat <- simrun(X=33, n_obs = n_obs, MX1 = MX1, MX2 = MX2, MX3 = MX3, lite_version = T))
