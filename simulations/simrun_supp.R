source('simulation-functions/load_all_funcs_and_libs.R')
source('ebcont/makeCA.R')
# ------------------------------------------------------------------------
# Simulation parameters ------
source('simulations/simulation-parameters-01.R')

# ------------------------------------------------------------------------

update_array <- function(mod, which_m, which_s, arry, TRUTH){
  arry[which_m, 1, which_s, ] <- mod$ests_loess - TRUTH
  arry[which_m, 2, which_s, ] <- mod$ests_lmmod1 - TRUTH
  arry[which_m, 3, which_s, ] <- mod$ests_lmmod2 - TRUTH
  arry[which_m, 4, which_s, ] <- mod$ests_splmod - TRUTH
  arry
}

verbose = F

set.seed(2020 + 1)
n_pts <- 100
A_test <- seq(1.5, 45, length.out = n_pts)
truth <- - 0.15 * A_test^2 + A_test * (2 + MX1^2 + MX2^2) - 15
truth <- truth / 50

n_sims <- 25
n_obs <- 1000

hm_array <- array(NA, c(4,4,n_sims,n_pts))
fm_array <- array(NA, c(4,4,n_sims,n_pts))

for(s in 1:n_sims){
  X1 <- rnorm(n_obs, mean = MX1, sd = 1)
  X2 <- rnorm(n_obs, mean = MX2, sd = 1)
  X3 <- rnorm(n_obs, mean = 0, sd = 1)
  X4 <- rnorm(n_obs, mean = MX2, sd = 1)
  X5 <- rbinom(n_obs, 1, prob = MX3)
  
  Z1 <- exp(X1 / 2)
  Z2 <- (X2 / (1 + exp(X1))) + 10
  Z3 <- (X1 * X3 / 25) + 0.6
  Z4 <- (X4 - MX2)**2 
  Z5 <- X5
  
  
  muA <- 5 * abs(X1) + 6 * abs(X2) + 3 * abs(X5) + abs(X4)
  
  A <- rchisq(n_obs, df = 3, ncp = muA)
  Y <- -0.15 * A^2 + A * (X1^2 + X2^2) - 15 + rnorm(n_obs, sd = 1) + (X1 + 3)^2 + 2 *(X2 - 25)^2 + X3 - (1 + MX1^2 + 6 * MX1 + 9) - 2 * (1 + MX2^2 - 50 * MX2 + 625)
  Y <- Y / 50
  
  datz <- data.frame(Y, A, 'X1' = Z1, 'X2' = Z2, 'X3' = Z3, 'X4' = Z4, 'X5' = Z5)
  
  # Weights Estimation -----------------------
  
  # Entopy Balancing
  CZ1 <- makeC2(datz[,3:7], datz$A, n_moments = 1)
  CZ2 <- makeC2(datz[,3:7], datz$A, n_moments = 2)
  CZ3 <- makeC2(datz[,3:7], datz$A, n_moments = 3)
  CZ4 <- makeC2(datz[,3:7], datz$A, n_moments = 4)
  
  CA1 <- makeCA(datz[,3:7], datz$A, n_moments = 1)
  CA2 <- makeCA(datz[,3:7], datz$A, n_moments = 2)
  CA3 <- makeCA(datz[,3:7], datz$A, n_moments = 3)
  CA4 <- makeCA(datz[,3:7], datz$A, n_moments = 4)
  
  ebz1 <- entbal_fit(CZ1, rep(0,ncol(CZ1)), 
                     n_moments = 1, 
                     verbose = verbose, opt_constraints = c(-250,250),
                     bal_tol = 1e-8)
  
  ebz2 <- entbal_fit(CZ2, rep(0,ncol(CZ2)), 
                     n_moments = 1, 
                     verbose = verbose, opt_constraints = c(-250,250),
                     bal_tol = 1e-8)
  
  ebz3 <- entbal_fit(CZ3, rep(0,ncol(CZ3)), 
                     n_moments = 1, 
                     verbose = verbose, opt_constraints = c(-250,250),
                     bal_tol = 1e-8)
  
  ebz4 <- entbal_fit(CZ4, rep(0,ncol(CZ4)), 
                     n_moments = 1, 
                     verbose = verbose, opt_constraints = c(-250,250),
                     bal_tol = 1e-8)
  
  eba1 <- entbal_fit(CA1, rep(0,ncol(CA1)), 
                     n_moments = 1, 
                     verbose = verbose, opt_constraints = c(-250,250),
                     bal_tol = 1e-8)
  
  eba2 <- entbal_fit(CA2, rep(0,ncol(CA2)), 
                     n_moments = 1, 
                     verbose = verbose, opt_constraints = c(-250,250),
                     bal_tol = 1e-8)
  
  eba3 <- entbal_fit(CA3, rep(0,ncol(CA3)), 
                     n_moments = 1, 
                     verbose = verbose, opt_constraints = c(-250,250),
                     bal_tol = 1e-8)
  
  eba4 <- entbal_fit(CA4, rep(0,ncol(CA4)), 
                     n_moments = 1, 
                     verbose = verbose, opt_constraints = c(-250,250),
                     bal_tol = 1e-8)
  
  # Proposed Method
  
  mod_eb1 <- model_estimation(datz, ebz1$wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test, nm = 1)
  mod_eb2 <- model_estimation(datz, ebz2$wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test, nm = 2)
  mod_eb3 <- model_estimation(datz, ebz3$wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test, nm = 3)
  mod_eb4 <- model_estimation(datz, ebz4$wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test, nm = 4)
  
  mod_eb1a <- model_estimation(datz, eba1$wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test, nm = 1)
  mod_eb2a <- model_estimation(datz, eba2$wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test, nm = 2)
  mod_eb3a <- model_estimation(datz, eba3$wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test, nm = 3)
  mod_eb4a <- model_estimation(datz, eba4$wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test, nm = 4)

  # Update results  
  hm_array <- update_array(mod_eb1, 1, s, hm_array, truth)
  hm_array <- update_array(mod_eb2, 2, s, hm_array, truth)
  hm_array <- update_array(mod_eb3, 3, s, hm_array, truth)
  hm_array <- update_array(mod_eb4, 4, s, hm_array, truth)

  fm_array <- update_array(mod_eb1a, 1, s, fm_array, truth)
  fm_array <- update_array(mod_eb2a, 2, s, fm_array, truth)
  fm_array <- update_array(mod_eb3a, 3, s, fm_array, truth)
  fm_array <- update_array(mod_eb4a, 4, s, fm_array, truth)
  
  message(paste(s,'.', sep=''), appendLF = F)
  if(!(s %% 10)) message('', appendLF = T)
}

# saveRDS(hm_array, 'output/hm_array.RDS')
# saveRDS(fm_array, 'output/fm_array.RDS')

hm_array <- readRDS('output/hm_array.RDS')
fm_array <- readRDS('output/fm_array.RDS')


get_model <- function(i, print_mdl = T){
  if(i == 1) mdl <- 'LOESS'
  if(i == 2) mdl <- 'LMMOD1'
  if(i == 3) mdl <- 'LMMOD2'
  if(i == 4) mdl <- 'BSPLINE'  
  if(print_mdl) message(mdl)
  invisible(mdl)
}


restab <- matrix(NA, nrow = 16, ncol = 7)
counter <- 0
for(mdl in 1:4){
  mdl_name <- get_model(mdl)
  for(mnts in 1:4){
    counter <- counter + 1
    avg_bias1 <- mean(apply(hm_array[mnts, mdl, , ],2,mean))
    avg_mse1 <- mean(apply(hm_array[mnts, mdl, , ]^2,2,mean))
    avg_bias2 <- mean(apply(fm_array[mnts, mdl, , ],2,mean))
    avg_mse2 <- mean(apply(fm_array[mnts, mdl, , ]^2,2,mean))
    ratio_bias <- avg_bias2 / avg_bias1
    restab[counter, ] <- c(mdl_name, 
                           mnts, 
                           round(avg_bias1, 3), 
                           round(avg_mse1, 3),
                           round(avg_bias2, 3), 
                           round(avg_mse2, 3),
                           round(ratio_bias,3))
  }
}
colnames(restab) <- c('Model', 'Moments', 'HM-AvgBias', 'HM-AvgMSE', 'FM-AvgBias', 'FM-AvgMSE', 'Ratio-Bias')
knitr::kable(restab)
xtable::xtable(restab[,1:6], )
