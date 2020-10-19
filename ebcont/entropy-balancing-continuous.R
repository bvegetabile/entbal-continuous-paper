makeC2 <- function(XD,A, n_moments = 3){
  NC <- ncol(XD)
  NR <- nrow(XD)
  outmat <- matrix(NA, nrow = NR, ncol = 0)
  #covariance 
  for(c in 1:NC){
    nuniq <- length(unique(XD[,c]))
    if(nuniq <= 1) {
      message(paste('Column: ', 
                    colnames(XD)[c], 
                    ', has <= 1 unique value and was not included'))
    } else if(nuniq == 2) {
      outmat <- cbind(outmat, scale(XD[,c]) * scale(A))
    } else {
      for(p in 1:n_moments){
        outmat <- cbind(outmat, scale(XD[,c]**p) * scale(A))
      }  
    }
  }
  outmat <- cbind(outmat, poly(A, n_moments))
  for(c in 1:NC){
    nuniq <- length(unique(XD[,c]))
    if(nuniq == 2 ) {
      outmat <- cbind(outmat, scale(XD[,c])) 
    } else {
      outmat <- cbind(outmat, poly(XD[,c], n_moments))   
    }
    
  }
  
  outmat
}

# Weight Construction for Entropy Balancing
entbal_wts <- function(Q, C, Z){
  norm_c <- Q %*% exp( - C %*% Z )
  Q * exp( - C %*% Z ) / c(norm_c)
}

# Stand Alone implementation of Entropy Balancing
entbal_fit <- function(C, targets,
                       n_moments = 2,
                       max_iters = 1000,
                       verbose = 0,
                       optim_method = 'L-BFGS-B',
                       bal_tol = 0.0005,
                       opt_constraints = c(-100, 100)){
  n_obs <- nrow(C)
  Q <- rep(1/n_obs, n_obs)
  M <- targets
  n_targets <- length(M)
  
  loss_func0 <- function(f){
    loss <- log(t(Q) %*% exp( - C %*% f )) + t(M) %*% f
    return(loss)
  }
  
  grad_func0 <- function(f){
    W <- entbal_wts(Q, C, f)
    grad <- M - t(C) %*% W
    return(grad)
  }
  
  f_init <- solve(t(C) %*% C + diag(ncol(C))) %*% M
  
  if(optim_method == 'L-BFGS-B'){
    opt_val <- optim(par = f_init,
                     fn = loss_func0,
                     gr = grad_func0,
                     method = optim_method,
                     lower = opt_constraints[1],
                     upper = opt_constraints[2],
                     control = list(trace = verbose,
                                    maxit = max_iters,
                                    lmm = 5,
                                    pgtol = bal_tol))
  } else if (optim_method == 'BFGS') {
    opt_val <- optim(par = f_init,
                     fn = loss_func0,
                     gr = grad_func0,
                     method = optim_method,
                     control = list(trace = verbose,
                                    maxit = max_iters))
  } else {
    stop('Unknown optimization method: Only L-BFGS-B and BFGS supported at this point')
  }
  
  return(list(optim_obj = opt_val,
              f = opt_val$par,
              wts = entbal_wts(Q, C, opt_val$par)))
  
}
