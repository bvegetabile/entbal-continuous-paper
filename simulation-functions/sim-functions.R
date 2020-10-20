# Given a simulation data struction and a set of weights:
# - Estimate using LOESS
# - Estimate linear regression 
# - Estimate polynomial regression 

model_estimation <- function(dat, 
                             wts=rep(1/nrow(dat), 
                                     nrow(dat)),
                             poly_reg = T, p = 1,
                             testpts, nm = 2,
                             which_model = 'EB'){
  
  span <- loess_split_sample(dat$Y, dat$A, dat[,3:7],
                             spans = seq(0.05,1,0.05),
                             n_moments = nm,
                             plot_curves = F,
                             model = which_model)
  dfs <- bspline_split_sample(dat$Y, dat$A, dat[,3:7], n_moments = nm,
                              plot_curves = F, model = which_model)
  loessmod <- loess(Y ~ A,
                    weights = wts,
                    data = dat,
                    degree = 1,
                    span = span$best_span,
                    control = loess.control(surface = 'direct'))
  lmmod1 <- lm(Y ~ A, data = dat, weights = wts)
  lmmod2 <- lm(Y ~ poly(A,2), data = dat, weights = wts)
  arange <- c(0, max(dat$A))
  splmod <- lm(Y ~ bs(A, df = dfs$best_df, Boundary.knots = arange), data = dat, weights = wts)
  
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
      'regbal' = apply(dat[3:7], 2, function(x) lm_ps(x, dat$A, wts = wts)$ests[2,]),
      'cors' = apply(dat[,3:7], 2, function(x) wcor(wts, x, dat[,2]))
    )
  )
}


# Bootstrap estimation ----

bsest <- function(X, simdat, testpts) {
  bsamp <- sample(1:nrow(simdat), nrow(simdat), replace = T)
  datbs <- simdat[bsamp, ]
  
  C <- makeC2(datbs[,3:7], datbs$A, n_moments = 3)
  eb <- entbal::entbal_fit(C, rep(0,ncol(C)),
                           n_moments = 1, 
                           verbose = F,
                           bal_tol = 1e-8)
  datbs$W <- eb$wts
  # 
  
  bspan <- loess_cv(datbs$Y,
                    datbs$A,
                    eb$wts,
                    nfolds = 2,
                    spans = seq(0.1, 1, 0.1),
                    plot_curves = F)$best_span
  
  loessmod <- loess(Y ~ A, 
                    weights = W,
                    data = datbs, 
                    degree = 1, 
                    span = bspan,
                    control = loess.control(surface = 'direct'))
  
  pls <- predict(loessmod, newdata = data.frame(A = testpts))
  pls
}


# LOESS Cross-Validation
loess_cv <- function(outcome, 
                     treatment, 
                     weights,
                     nfolds = 5,
                     spans = seq(0.01, 0.5, length.out = 25),
                     plot_curves = T){
  
  n_obs <- length(outcome)
  
  dat <- data.frame(Y = outcome, 
                    A = treatment, 
                    W = weights)
  
  folds <- sample(rep(1:nfolds, length.out = n_obs), n_obs, replace = F)
  
  resmat <- array(NA, c(length(spans), nfolds, 3))
  for(s in 1:length(spans)){
    for(f in 1:nfolds){
      span <- spans[s]
      traindat <- dat[folds != f, ]
      testdat <- dat[folds == f, ]
      
      trainmod <- loess(Y ~ A,
                        weights = W, 
                        data = traindat, 
                        degree = 1, 
                        span = span,
                        control = loess.control(surface = 'direct'))
      
      trainpreds <- predict(trainmod, newdata = traindat)
      testpreds <- predict(trainmod, newdata = data.frame(A=testdat$A))
      
      trainerr <- wmse(traindat$Y, trainpreds, traindat$W)
      testerr <- wmse(testdat$Y, testpreds, testdat$W)
      
      resmat[s,f,] <- c(span, trainerr, testerr)
    }
  }
  
  traintest_errs <- apply(resmat, 1, colMeans)
  best_span <- traintest_errs[1, which(traintest_errs[3,] == min(traintest_errs[3,]))]
  
  if(plot_curves){
    par(mfrow=c(1,1))
    plot(traintest_errs[1,], traintest_errs[2,],
         ylim = range(traintest_errs[2:3, ]),
         type = 'l', col = rgb(0.75,0,0,0.5),
         lwd = 3)
    lines(traintest_errs[1,], traintest_errs[3,],
          col = rgb(0,0,0.75,0.5),
          lwd = 3)
  }
  
  
  list('best_span' = best_span,
       'cv_errs' = traintest_errs)
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
