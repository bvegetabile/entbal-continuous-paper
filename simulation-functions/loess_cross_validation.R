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