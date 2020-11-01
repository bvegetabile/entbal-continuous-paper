library(gbm)

bal_by_trees <- function(model, treatment, covs, trees = 1000, print_bal = FALSE){
  preds <- predict(model, newdata = covs, n.trees = trees)
  sd_resids <- sd(treatment - preds)
  nums <- dnorm(treatment, mean = mean(treatment), sd = sd(treatment))
  denoms <- dnorm(treatment, mean = preds, sd = sd_resids)
  wts <- nums / denoms
  dim_bal <- abs(apply(covs, 2, function(x) wcor(wts, x, treatment)))
  if(print_bal) print(dim_bal)
  bal <- sum(dim_bal)
  bal 
}

get_weights <- function(model, treatment, covs, trees){
  preds <- predict(model, newdata = covs, n.trees = trees)
  sd_resids <- sd(treatment - preds)
  nums <- dnorm(treatment, mean = mean(treatment), sd = sd(treatment))
  denoms <- dnorm(treatment, mean = preds, sd = sd_resids)
  wts <- nums / denoms
  wts
}

gbm_weights <- function(treatment, 
                        covariates,
                        nt = 1000){
  mod <- gbm.fit(covariates, treatment, 
                 distribution = 'gaussian', n.trees = nt, 
                 shrinkage = 0.01, interaction.depth = 3,
                 verbose = F)
  tree_seq <- seq(25, nt, 25)
  system.time(bal_path <- sapply(seq(25, nt, 25), function(x) bal_by_trees(mod, treatment, covariates, x, print_bal = FALSE)))
  best_trees <- tree_seq[which(bal_path == min(bal_path))]
  wts <- get_weights(mod, treatment, covariates, best_trees)
  list('wts' = wts,
       'trees' = best_trees)
}

# WTS <- gbm_weights(testdat$dat$A, testdat$dat[,3:7])
# 
# abs(apply(testdat$dat[,3:7], 2, function(x) wcor(WTS$wts, x, testdat$dat$A)))
