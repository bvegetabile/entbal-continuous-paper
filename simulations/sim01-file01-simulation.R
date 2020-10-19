source('simulation-functions/load_all_funcs_and_libs.R')

# ------------------------------------------------------------------------
# Simulation parameters ------
source('simulations/simulation-parameters-01.R')

# ------------------------------------------------------------------------

# #--- Visualizing one span solution
set.seed(20200109)
n_obs <- 1000
t1 <- system.time(testdat <- simrun(X=1, n_obs = n_obs, MX1 = MX1, MX2 = MX2, MX3 = MX3, lite_version = T))
system.time(mod1 <- loess_split_sample(testdat$dat$Y, testdat$dat$A, testdat$dat[,3:5], 
                                       n_moments = 2, spans = seq(0.05,1,0.05), plot_curves = T))
system.time(mod2 <- bspline_split_sample(testdat$dat$Y, testdat$dat$A, testdat$dat[,3:5], 
                                         n_moments = 2,  plot_curves = T, dfs = seq(3,10,1)))
outro <- model_estimation(testdat$dat, testdat$ebmod2$wts, testpts = A_test)
plot(A_test, truth, type = 'l', lwd = 3)
lines(A_test, outro$ests_loess, lwd = 5, col = 'blue')
lines(A_test, outro$ests_splmod, lwd = 5, col = 'red')


# ------------------------------------------------------------------------
# Simulation ------
n_sims <- 20
options(cores = 20)
cl <- makeCluster(20)
registerDoParallel(cl)
system.time(simout <- foreach(i = 1:n_sims, 
                              .packages = c('CBPS', 'splines'),
                              # .errorhandling = 'remove',
                              .verbose = T) 
            %dopar% simrun(i, n_obs = 1000,
                           MX1 = MX1, MX2 = MX2, MX3 = MX3, 
                           amin = amin, amax = amax,
                           a_effect = T, lite_version = F))
stopCluster(cl)

# Storing simulation results 

saveRDS(simout, file = 'output/sim1-output1-effect-estimation.RDS')
