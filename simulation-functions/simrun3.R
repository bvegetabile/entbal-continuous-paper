simrun_lite <- function(X, 
                   n_obs = 1000, 
                   npts = 100, 
                   verbose = F,
                   MX1 = -0.5,
                   MX2 = 1,
                   MX3 = 0.3,
                   amin = 0,
                   amax = 45,
                   a_effect = T,
                   nm = 3,
                   lowspan = 0.05){
  # Data generation --------------------------
  
  set.seed(2020 + X)
  
  A_test <- seq(amin, amax, length.out = npts)
  
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
  
  if(a_effect){
    Y <- - 0.15 * A^2 + A * (X1^2 + X2^2) - 15 + (X1+3)^2 + 2 * (X2-25)^2 + X3 + rnorm(n_obs, sd = 1) - ((1 + MX1^2 + 6 * MX1 + 9) + 2 * (1 + MX2^2 - 50 * MX2 + 625))
    Y <- Y / 50
    truth <- - 0.15 * A_test^2 + A_test * (2 + MX1^2 + MX2^2) - 15 #+ 5 * (1 + MX1^2 + 6 * MX1 + 9) + 15 * (1 + MX2^2 + 6 * MX2 + 9) + MX3
    truth <- truth / 50
    # truth <- - (A_test - 10) * (A_test - 10) / 5 + 5* A_test * (2 + MX1^2 + MX2^2) / 1 - 15 * (MX1 - MX2) + MX3 + (1 + MX2^2 - 40 * MX2 + 400)
  } else {
    Y <- X1 + X1^2 + X2 + X2^2 + X1 * X2 + X5 + rnorm(n_obs, sd = 1)
    truth <- rep(MX1 + (MX1^2 + 1) + MX2 + (MX2^2 + 1) + MX1 * MX2 + MX3 , 100)
  }
  
  datz <- data.frame(Y, A, 'X1' = Z1, 'X2' = Z2, 'X3' = Z3, 'X4' = Z4, 'X5' = Z5)
  
  # Weights Estimation -----------------------
  
  # Entopy Balancing
  CZM <- makeC2(datz[,3:7], datz$A, n_moments = nm)
  
  ebz1 <- entbal_fit(CZM, rep(0,ncol(CZM)), 
                     n_moments = 1, 
                     verbose = verbose, opt_constraints = c(-250,250),
                     bal_tol = 1e-8)

  # Proposed Method
  
  mod_eb1 <- model_estimation_lite(datz, ebz1$wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test, nm = nm, ls = lowspan)

  return(  
    list('eb1' = mod_eb1,
         'ebmod1' = ebz1,
         'dat' = datz,
         'truth' = truth)
  )
}