
n_obs <- 200

source('simulations/simulation-parameters-01.R')

get_quantsA <- function(x){
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
  quantile(A, c(0.01, 0.05, 0.95, 0.99))
}

apply(sapply(1:10000, get_quantsA), 1, mean)
