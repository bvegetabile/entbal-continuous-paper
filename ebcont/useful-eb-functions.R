# weighted functions --------------

# Weighted Mean
wmean <- function(wts, x){
  sum(x * wts) / sum(wts)
}

# Weighted Variance
wvar <- function(wts, x){
  wm <- wmean(wts, x)
  wv <- sum(wts * (x - wm)^2) / sum(wts)
  wv
}

# Weighted Correlation
wcor <- function(wts, x, y){
  wmx <- wmean(wts, x)
  wmy <- wmean(wts, y)
  wvx <- wvar(wts, x)
  wvy <- wvar(wts, y)
  topval <- sum(wts * (x - wmx) * (y - wmy)) / sum(wts)
  topval / sqrt(wvy * wvx)
}

# Compare the weighted KS-Statistic for X and weighted X
check_ks_stat <- function(X, wts){
  B1 <- ecdf(X)
  B2 <- entbal::wtd_ecdf(X, wts)
  XTEST <- seq(min(X), max(X), length.out = 250)
  max(abs(B1(XTEST) - B2(XTEST)))
}

# Weighted Covariance
wcov <- function(wts, x, y){
  wmx <- wmean(wts, x)
  wmy <- wmean(wts, y)
  wvx <- wvar(wts, x)
  wvy <- wvar(wts, y)
  topval <- sum(wts * (x - wmx) * (y - wmy)) / sum(wts)
  topval
}

# Weighted Mean Squared Error
wmse <- function(y, yhat, w) {
  sum(w * (y - yhat)^2) / sum(w)
}

# Other Functions ------------

# Mean Squared Error
mse <- function(y,yhat) mean((y-yhat)^2)
