# Given a simulation data struction and a set of weights:
# - Estimate using LOESS
# - Estimate linear regression 
# - Estimate polynomial regression 

model_estimation <- function(dat, 
                             wts=rep(1/nrow(dat), 
                                     nrow(dat)),
                             poly_reg = T, p = 1,
                             testpts, nm = 2,
                             which_model = 'EB',
                             lite=F){
  
  span <- loess_split_sample(dat$Y, dat$A, dat[,3:7],
                             spans = seq(0.05,1,0.05),
                             n_moments = nm,
                             plot_curves = F,
                             model = which_model)
  loessmod <- loess(Y ~ A,
                    weights = wts,
                    data = dat,
                    degree = 1,
                    span = span$best_span,
                    control = loess.control(surface = 'direct'))
  
  dfs <- bspline_split_sample(dat$Y, dat$A, dat[,3:7], n_moments = nm,
                              plot_curves = F, model = which_model)    
  
  lmmod1 <- lm(Y ~ A, data = dat, weights = wts)
  lmmod2 <- lm(Y ~ poly(A,2), data = dat, weights = wts)
  arange <- c(0, max(dat$A))
  splmod <- lm(Y ~ bSpline(A, df = dfs$best_df, Boundary.knots = arange), data = dat, weights = wts)
  
  return(
    list(
      'wts' = wts,
      'spanres' = span$best_span,
      'dfres' = dfs$best_df,
      'ests_loess' = predict(loessmod, newdata = data.frame(A = testpts)),
      'ests_lmmod1' = predict(lmmod1, newdata = data.frame(A = testpts)),
      'ests_lmmod2' = predict(lmmod2, newdata = data.frame(A = testpts)),
      'ests_splmod' = predict(splmod, newdata = data.frame(A = testpts)),
      'ess' = (sum(wts))^2 / sum(wts^2),
      'ksstats' = apply(dat[,2:7], 2, function(x) check_ks_stat(x, wts)),
      # 'regbal' = apply(dat[3:7], 2, function(x) lm_ps(x, dat$A, wts = wts)$ests[2,]),
      'cors' = apply(dat[,3:7], 2, function(x) wcor(wts, x, dat[,2]))
    )
  )
}

model_estimation_lite <- function(dat, 
                                  wts=rep(1/nrow(dat), 
                                          nrow(dat)),
                                  poly_reg = T, p = 1,
                                  testpts, nm = 2,
                                  which_model = 'EB',
                                  ls = 0.05,
                                  lite=F){
  
  span <- loess_split_sample(dat$Y, dat$A, dat[,3:7],
                             spans = seq(ls,1,0.05),
                             n_moments = nm,
                             plot_curves = F,
                             model = which_model)
  loessmod <- loess(Y ~ A,
                    weights = wts,
                    data = dat,
                    degree = 1,
                    span = span$best_span,
                    control = loess.control(surface = 'direct'))
  
  return(
    list(
      'wts' = wts,
      'spanres' = span$best_span,
      'ests_loess' = predict(loessmod, newdata = data.frame(A = testpts))
    )
  )
}



# Bootstrap estimation ----

bsest <- function(X, simdat, testpts, nm = 3) {
  bsamp <- sample(1:nrow(simdat), nrow(simdat), replace = T)
  datbs <- simdat[bsamp, ]
  
  C <- makeC2(datbs[,3:7], datbs$A, n_moments = nm)
  eb <- entbal_fit(C, rep(0,ncol(C)),
                   n_moments = 1, 
                   verbose = F,
                   bal_tol = 1e-8)
  datbs$W <- eb$wts
  
  bspan <- loess_split_sample(datbs$Y, datbs$A, datbs[,3:7],
                              spans = seq(0.05,1,0.05),
                              n_moments = nm,
                              plot_curves = F,
                              model = 'EB')$best_span
  
  loessmod <- loess(Y ~ A, 
                    weights = W,
                    data = datbs, 
                    degree = 1, 
                    span = bspan,
                    control = loess.control(surface = 'direct'))
  
  pls <- predict(loessmod, newdata = data.frame(A = testpts))
  pls
}

# Bootstrap estimation ----

bsest_smalln <- function(X, simdat, testpts, nm = 3, sval = 0.75) {
  bsamp <- sample(1:nrow(simdat), nrow(simdat), replace = T)
  datbs <- simdat[bsamp, ]
  
  C <- makeC2(datbs[,3:7], datbs$A, n_moments = nm)
  eb <- entbal_fit(C, rep(0,ncol(C)),
                   n_moments = 1, 
                   verbose = F,
                   bal_tol = 1e-5,
                   opt_constraints = c(-50,50))
  datbs$W <- eb$wts
  
  # bspan <- 0.75
  # bspan <- loess_split_sample(datbs$Y, datbs$A, datbs[,3:7],
  #                             spans = seq(0.5,1,0.05),
  #                             n_moments = nm,
  #                             plot_curves = F,
  #                             model = 'EB')$best_span
  bspan <- sval
  
  loessmod <- loess(Y ~ A, 
                    weights = W,
                    data = datbs, 
                    degree = 1, 
                    span = bspan,
                    control = loess.control(surface = 'direct'))
  
  pls <- predict(loessmod, newdata = data.frame(A = testpts))
  pls
}


# Simulation Analysis ----------

plot_sim <- function(estsmat, atest, truth, 
                     aquants, 
                     ylimits = c(-10,10), ylines = 2,
                     maintitle = '') {
  plot(atest, truth,
       ylim = ylimits,
       type = 'l', col = rgb(1,0,0,0.5), lwd =3,
       xlab = 'Treatment',
       ylab = 'Ests.',
       main = maintitle,
       axes = F); axis(1); axis(2, las = 2);
  abline(h = seq(ylimits[1],ylimits[2],ylines), col = rgb(0,0,0,0.15))
  abline(v = aquants, lty = c(2,3,3,2), col = rgb(0,0,0,0.5))
  lines(A_test, apply(estsmat,2,mean), lwd = 3, col = rgb(0,0,0.75,0.5))
  lines(A_test, apply(estsmat,2,quantile, 0.025), lty = 1, lwd = 2, col = rgb(0,0,0.75,0.25))
  lines(A_test, apply(estsmat,2,quantile, 0.975), lty = 1, lwd = 2, col = rgb(0,0,0.75,0.25))
  lines(A_test, apply(estsmat,2,quantile, 0.025), lty = 3, lwd = 2, col = rgb(0,0,0,0.75))
  lines(A_test, apply(estsmat,2,quantile, 0.975), lty = 3, lwd = 2, col = rgb(0,0,0,0.75))
}
