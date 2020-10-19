
# Weighted Least Squares with Variance Correction
lm_ps <- function(Y, x, wts, true_val = NULL){
  X <- cbind(1,x)
  W <- diag(as.vector(wts))
  invXtWX <- solve(t(X) %*% W %*% X)
  hatmat <- invXtWX %*% t(X) %*% W  
  betas <- hatmat %*% Y
  pseudodf <- sum(wts)^2 / sum(wts^2)
  
  Yhat <- X %*% betas
  resids <- Y - Yhat
  
  # varmat <- sighat * invXtWX %*% t(X) %*% W %*% W %*% X %*% invXtWX
  varmat <- invXtWX %*% t(X) %*% W %*% diag(as.vector(resids)^2) %*% W %*% X %*% invXtWX
  std_errs <- sqrt(diag(varmat))
  
  low_int <- betas - 1.96 * std_errs
  upp_int <- betas + 1.96 * std_errs
  
  res <- cbind(betas, std_errs, betas/std_errs, 
               2 * (1 - pt(abs(betas/std_errs), df = pseudodf - 1)),
               low_int, upp_int)
  colnames(res) <- c('coef', 'stderrs', 't-value', 'p-value', 'low95', 'upp95')
  
  if(!is.null(true_val)){
    cover_test <- res[2,3] < true_val & res[2,4] > true_val
    return(list('ests' = data.frame(res),
                'covers' = cover_test))
  } else{
    return(list('ests' = res))
  }
}