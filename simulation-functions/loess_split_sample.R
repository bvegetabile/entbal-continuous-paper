
# LOESS Cross-Validation
loess_split_sample <- function(outcome, 
                               treatment, 
                               covariates,
                               n_moments, 
                               model = 'EB',
                               spans = seq(0.01, 0.5, length.out = 25),
                               plot_curves = T){
  nfolds <- 2
  n_obs <- length(outcome)
  
  dat <- data.frame(Y = outcome, 
                    A = treatment, 
                    covariates)
  dat$W <- NA
  
  folds <- sample(rep(1:2, length.out = n_obs), n_obs, replace = F)
  
  resmat <- array(NA, c(length(spans), nfolds, 3))
  for(f in 1:nfolds){
    CZ <- makeC2(covariates[folds == f, ], dat$A[folds == f], n_moments = n_moments)
    if(model == 'EB'){
      ebmod <- entbal_fit(CZ, rep(0,ncol(CZ)), 
                          n_moments = 1, 
                          verbose = FALSE, opt_constraints = c(-250,250),
                          bal_tol = 1e-8)  
      dat$W[folds == f] <- ebmod$wts
    } else if (model == 'npCBPS') {
      if(n_moments == 1){
        junkmsg <- capture.output(npcbpsmod <- CBPS::npCBPS(A ~ X1 + X2 + X3 + X4 + X5, 
                                                            data = dat[folds == f, ], 
                                                            corprior = 1e-8, print.level = 0))  
      } else {
        junkmsg <- capture.output(npcbpsmod <- CBPS::npCBPS(A ~ poly(X1,n_moments) + poly(X2,n_moments) + poly(X3,n_moments) + poly(X4,n_moments) + X5, 
                                                            data = dat[folds == f, ], 
                                                            corprior = 1e-8, print.level = 0))
      }
      dat$W[folds == f] <- npcbpsmod$weights
    } else if (model == 'CBPS'){
      cbpsmod <- CBPS::CBPS(A ~ X1 + X2 + X3 + X4 + X5, data = dat[folds==f,], print.level = 0, method = 'exact')
      dat$W[folds == f] <- cbpsmod$weights
    } else if (model == 'none'){
      dat$W[folds == f] <- 1 / length(dat$W[folds == f])
    } else if ( model == 'linear'){
      Af <- dat$A[folds == f]
      ps_lm_mod <- lm(A ~ X1 + X2 + X3 + X4 + X5, data = dat[folds == f,])
      ps_lm <- predict(ps_lm_mod)
      lm_wts <- dnorm(Af, mean = mean(Af), sd = sd(Af)) / dnorm(Af, mean = ps_lm, sd = sd(ps_lm_mod$residuals))
      lm_wts <- lm_wts / sum(lm_wts)
      dat$W[folds == f] <- lm_wts
    }
  }
  for(f in 1:nfolds){
    for(s in 1:length(spans)){
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
