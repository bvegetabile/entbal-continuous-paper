source('simulation-functions/load_all_funcs_and_libs.R')

# ------------------------------------------------------------------------
# Simulation parameters ------
source('simulations/simulation-parameters-02.R')

# #--- Visualizing one span solution
set.seed(20200109)
n_obs <- 1000
t1 <- system.time(testdat <- simrun(X=33, n_obs = n_obs, MX1 = MX1, MX2 = MX2, MX3 = MX3, lite_version = T, a_effect = FALSE))
print(t1)

system.time(mod1 <- loess_split_sample(testdat$dat$Y, testdat$dat$A, testdat$dat[,3:7], 
                                       n_moments = 3, spans = seq(0.05,1,0.05), plot_curves = T))
system.time(mod2 <- bspline_split_sample(testdat$dat$Y, testdat$dat$A, testdat$dat[,3:7], 
                                         n_moments = 3,  plot_curves = T, dfs = seq(3,10,1)))

outro <- model_estimation(testdat$dat, testdat$eb2$wts, testpts = A_test)
naive <- model_estimation(testdat$dat, testpts = A_test)
plot(testdat$dat$A, testdat$dat$Y, pch = 19, col = rgb(0,0,0,0.25))
lines(A_test, truth, type = 'l', lwd = 3, col = 'red')
lines(A_test, testdat$truth, type = 'l', lwd = 3, col = 'red')
lines(A_test, outro$ests_loess, lwd = 5, col = 'blue')
lines(A_test, outro$ests_splmod, lwd = 5, col = 'darkgreen')
lines(A_test, outro$ests_lmmod2, lwd = 5, col = 'orange')
lines(A_test, naive$ests_lmmod2, lwd = 5, col = 'lightblue')


# ------------------------------------------------------------------------
# Simulation ------
n_sims <- 1000
options(cores = 20)
cl <- makeCluster(20)
registerDoParallel(cl)
system.time(simout <- foreach(i = 1:n_sims, 
                              .packages = c('CBPS', 'splines', 'gbm'),
                              # .errorhandling = 'remove',
                              .verbose = T) 
            %dopar% simrun(i, n_obs = 1000,
                           MX1 = MX1, MX2 = MX2, MX3 = MX3, 
                           amin = amin, amax = amax,
                           a_effect = F, lite_version = F))
stopCluster(cl)

# Storing simulation results 

saveRDS(simout, file = 'output/sim2-output1-effect-estimation.RDS')
