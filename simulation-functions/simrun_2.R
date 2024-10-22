simrun <- function(X, 
                   n_obs = 1000, 
                   npts = 100, 
                   verbose = F,
                   MX1 = -0.5,
                   MX2 = 1,
                   MX3 = 0.3,
                   amin = 0,
                   amax = 45,
                   a_effect = T,
                   lite_version = F){
  # Data generation --------------------------
  
  set.seed(2020 + X)
  
  A_test <- seq(amin, amax, length.out = npts)
  
  X1 <- rnorm(n_obs, mean = MX1, sd = 1)
  X2 <- rnorm(n_obs, mean = MX2, sd = 1)
  X3 <- rbinom(n_obs, 1, prob = MX3)
  
  muA <- 5 * abs(X1) + 6 * abs(X2) + 3 * abs(X3)
  
  A <- rchisq(n_obs, df = 3, ncp = muA)
  
  if(a_effect){
    Y <- - (A + 5) * (A - 5) / 300 + A * (X1^2 + X2^2) / 25 + X1 + X2 + X3 + rnorm(n_obs, sd = 1)
    truth <- - (A_test + 5) * (A_test - 5) / 300 + A_test * (2 + MX1^2 + MX2^2) / 25 + MX1 + MX2 + MX3
  } else {
    Y <- X1 + X1^2 + X2 + X2^2 + X1 * X2 + X3 + rnorm(n_obs, sd = 1)
    truth <- rep(MX1 + (MX1^2 + 1) + MX2 + (MX2^2 + 1) + MX1 * MX2 + MX3 , 100)
  }
  
  datz <- data.frame(Y, A, X1, X2, X3)
  
  # Weights Estimation -----------------------
  
  # Entopy Balancing
  CZ1 <- makeC2(datz[,3:5], datz$A, n_moments = 1)
  CZ2 <- makeC2(datz[,3:5], datz$A, n_moments = 2)
  CZ3 <- makeC2(datz[,3:5], datz$A, n_moments = 3)
  CZ4 <- makeC2(datz[,3:5], datz$A, n_moments = 4)
  
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
  
  datz$W1 <- ebz1$wts
  datz$W2 <- ebz2$wts
  datz$W3 <- ebz1$wts
  datz$W4 <- ebz2$wts
  
  # Proposed Method
  
  mod_eb1 <- model_estimation(datz, ebz1$wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test)
  mod_eb2 <- model_estimation(datz, ebz2$wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test)
  mod_eb3 <- model_estimation(datz, ebz3$wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test)
  mod_eb4 <- model_estimation(datz, ebz4$wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test)
  
  if(!lite_version){
    # linear model weights
    ps_lm_mod <- lm(A ~ X1 + X2 + X3, data = datz)
    ps_lm <- predict(ps_lm_mod)
    lm_wts <- dnorm(A, mean = mean(A), sd = sd(A)) / dnorm(A, mean = ps_lm, sd = sd(ps_lm_mod$residuals))
    lm_wts <- lm_wts / sum(lm_wts)
    # CBPS 
    
    cbpsmod <- CBPS::CBPS(A ~ X1 + X2 + X3, data = datz, print.level = 0, method = 'exact')
    
    junkmsg <- capture.output(npcbpsmod1 <- CBPS::npCBPS(A ~ X1 + X2 + X3, 
                                                         data = datz, 
                                                         corprior = 1e-8, print.level = 0))
    
    junkmsg <- capture.output(npcbpsmod2 <- CBPS::npCBPS(A ~ poly(X1,2) + poly(X2,2) + X3, 
                                                         data = datz, 
                                                         corprior = 1e-8, print.level = 0))
    
    junkmsg <- capture.output(npcbpsmod3 <- CBPS::npCBPS(A ~ poly(X1,3) + poly(X2,3) + X3, 
                                                         data = datz, 
                                                         corprior = 1e-8, print.level = 0))
    
    junkmsg <- capture.output(npcbpsmod4 <- CBPS::npCBPS(A ~ poly(X1,4) + poly(X2,4) + X3, 
                                                         data = datz, 
                                                         corprior = 1e-8, print.level = 0))
    
    mod_unwtd <- model_estimation(datz, poly_reg = T, wts = rep(1/n_obs, n_obs),  p = ifelse(a_effect, 2, 1), A_test,
                                  which_model = 'none')
    mod_lm <- model_estimation(datz, lm_wts, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test,
                               which_model = 'linear')
    mod_cpbs <- model_estimation(datz, cbpsmod$weights, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test, which_model = 'CBPS')
    mod_npcbps1 <- model_estimation(datz, npcbpsmod1$weights, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test, which_model = 'npCBPS')
    mod_npcbps2 <- model_estimation(datz, npcbpsmod2$weights, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test, which_model = 'npCBPS')
    mod_npcbps3 <- model_estimation(datz, npcbpsmod3$weights, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test, which_model = 'npCBPS')
    mod_npcbps4 <- model_estimation(datz, npcbpsmod4$weights, poly_reg = T, p = ifelse(a_effect, 2, 1), A_test, which_model = 'npCBPS')
    
    return(
      list(
        'unwtd' = mod_unwtd,
        'eb1' = mod_eb1,
        'eb2' = mod_eb2,
        'eb3' = mod_eb3,
        'eb4' = mod_eb4,
        'lm' = mod_lm,
        'cbps' = mod_cpbs,
        'npcbps1' = mod_npcbps1,
        'npcbps2' = mod_npcbps2,
        'npcbps3' = mod_npcbps3,
        'npcbps4' = mod_npcbps4,
        'dat' = datz,
        'truth' = truth)
    )
  } else {
    return(  
      list('eb1' = mod_eb1,
           'eb2' = mod_eb2,
           'eb3' = mod_eb3,
           'eb4' = mod_eb4,
           'ebmod1' = ebz1,
           'ebmod2' = ebz2,
           'ebmod3' = ebz3,
           'ebmod4' = ebz4,
           'dat' = datz,
           'truth' = truth)
    )
  }
}