source('simulation-functions/load_all_funcs_and_libs.R')

# ------------------------------------------------------------------------
# Simulation parameters ------
source('simulations/simulation-parameters-01.R')
# ------------------------------------------------------------------------------

update_results <- function(mod, bias_array, which_s, which_m){
  bias_array[which_m, 1, which_s, ] <- mod$ests_loess - truth
  bias_array[which_m, 2,  which_s, ] <- mod$ests_lmmod1 - truth
  bias_array[which_m, 3,  which_s, ] <- mod$ests_lmmod2 - truth
  bias_array[which_m, 4, which_s, ] <- mod$ests_splmod - truth
  bias_array
}

MX1 <- -0.5
MX2 <- 1
MX3 <- 0.3

amin <- 1.5
amax <- 45
n_obs <- 1000

A_test <- seq(amin, amax, length.out = 100)
truth <- - 0.15 * A_test^2 + A_test * (2 + MX1^2 + MX2^2) - 15 #+ 5 * (1 + MX1^2 + 6 * MX1 + 9) + 15 * (1 + MX2^2 + 6 * MX2 + 9) + MX3
truth <- truth / 50

# ------------------------------------------------------------------------------
verbose = F
n_sims <- 100
biasres1 <- array(NA, c(4, 4, n_sims, n_obs))
biasres2 <- array(NA, c(4, 4, n_sims, n_obs))
for(s in 1:n_sims){
  set.seed(2020 + s)  
  
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
  
  Y <- - 0.15 * A^2 + A * (X1^2 + X2^2) - 15 + (X1+3)^2 + 2 * (X2-25)^2 + X3 + rnorm(n_obs, sd = 1) - ((1 + MX1^2 + 6 * MX1 + 9) + 2 * (1 + MX2^2 - 50 * MX2 + 625))
  Y <- Y / 50
  truth <- - 0.15 * A_test^2 + A_test * (2 + MX1^2 + MX2^2) - 15 #+ 5 * (1 + MX1^2 + 6 * MX1 + 9) + 15 * (1 + MX2^2 + 6 * MX2 + 9) + MX3
  truth <- truth / 50
  
  
  datz <- data.frame(Y, A, 'X1' = Z1, 'X2' = Z2, 'X3' = Z3, 'X4' = Z4, 'X5' = Z5)
  

  
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
  
  mod_eb1 <- model_estimation(datz, ebz1$wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test)
  mod_eb2 <- model_estimation(datz, ebz2$wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test)
  mod_eb3 <- model_estimation(datz, ebz3$wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test)
  mod_eb4 <- model_estimation(datz, ebz4$wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test)
  
  mod_eb1a <- model_estimation(datz, eba1$wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test)
  mod_eb2a <- model_estimation(datz, eba2$wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test)
  mod_eb3a <- model_estimation(datz, eba3$wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test)
  mod_eb4a <- model_estimation(datz, eba4$wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test)
  
  biasres1 <- update_results(mod_eb1, biasres1, s, 1)
  biasres1 <- update_results(mod_eb2, biasres1, s, 2)
  biasres1 <- update_results(mod_eb3, biasres1, s, 3)
  biasres1 <- update_results(mod_eb4, biasres1, s, 4)
  
  biasres2 <- update_results(mod_eb1a, biasres2, s, 1)
  biasres2 <- update_results(mod_eb2a, biasres2, s, 2)
  biasres2 <- update_results(mod_eb3a, biasres2, s, 3)
  biasres2 <- update_results(mod_eb4a, biasres2, s, 4)
  message(paste(s,'.', sep=''), appendLF = F)
  if(!s %% 10) message('', appendLF = T)
}

dim(biasres1[1,1,,])

print_method <- function(x, printval = T){
  outval <- NA
  if(x == 1) outval <- 'LOESS'
  if(x == 2) outval <- 'LMMOD1'
  if(x == 3) outval <- 'LMMOD2'
  if(x == 4) outval <- 'BSLINE'
  if(printval) print(outval)
  invisible(outval)
}

for(i in 1:4){
  # methods
  print_method(i)
  for(j in 1:4){
    # moments
    outro <- c(
      mean(apply(biasres1[j,i,,], 1, mean, na.rm = T)),
      mean(apply(biasres2[j,i,,], 1, mean, na.rm = T))  
    )
    print(round(outro,4))
  }
}



