source('simulation-functions/load_all_funcs_and_libs.R')

# ------------------------------------------------------------------------
# Simulation parameters ------
source('simulations/simulation-parameters-02.R')

# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Reading Bootstrop Simulations ------
bsdata <- readRDS('output/sim02-bootstrap-results.rds')
bscovers <- bsdata$`bs-cov`
bsstderr <- bsdata$`bs-se`

simout <- readRDS('output/sim2-output1-effect-estimation.RDS')

ests <- array(NA, dim = c(12, length(simout), length(A_test)))
for(s in 1:length(simout)){
  out <- simout[[s]]
  ests[1,s,] <- out$unwtd$ests_loess
  ests[2,s,] <- out$eb1$ests_loess
  ests[3,s,] <- out$eb2$ests_loess
  ests[4,s,] <- out$eb3$ests_loess
  ests[5,s,] <- out$eb4$ests_loess
  ests[6,s,] <- out$lm$ests_loess
  ests[7,s,] <- out$cbps$ests_loess
  ests[8,s,] <- out$npcbps1$ests_loess
  ests[9,s,] <- out$npcbps2$ests_loess
  ests[10,s,] <- out$npcbps3$ests_loess
  ests[11,s,] <- out$npcbps4$ests_loess
  ests[12,s,] <- out$gbm$ests_loess
  
  if(! s %% 10) message(paste(s,':', sep =''), appendLF = F)
  if(! s %% 100) message('', appendLF = T)
}


# One BS --------

NB = 100
options(cores = 20)
#-----------------------------------------------------------------
cl <- makeCluster(20)
registerDoParallel(cl)
simres <- simrun(1, n_obs = 1000,
                 MX1 = MX1, MX2 = MX2, MX3 = MX3,
                 amin = amin, amax = amax,
                 a_effect = F,
                 lite_version = T)
bstime <- system.time(bsres <- foreach(
  i = 1:NB,
  .combine = rbind,
  .verbose = F) %dopar% bsest(i, simres$dat, A_test, nm = 3)
)
stopCluster(cl)
print(bstime)



# plot(A_test, apply(bscovers, 2, mean), ylim = c(0,1))
# abline(h = 0.95)
# # 
# # 
# plot(A_test, apply(bsstderr,2,mean)/apply(ests[4,,],2,sd))
# abline(h=1)
# # 
# # pdf('figures/bs-plot-one.pdf', height = 3.5, width = 12)
# par(mfrow=c(1,3))
# plot(A_test, truth, type = 'l', lwd = 3, col = rgb(0.75,0,0,0.75), ylim = c(-6, 5),
#      ylab = 'Outcome', xlab = 'Exposure Level',
#      main = 'Bootstrap Confidence Interval', axes = F)
# axis(1); axis(2, las = 2); box();
# abline(v = Aquants, lty = c(2,3,3,2), col = rgb(0,0,0,0.5))
# lines(A_test, simres$eb3$ests_loess, lwd = 3, col = rgb(0,0,0.75,0.5))
# lines(A_test, simres$eb3$ests_loess + 1.96 * apply(bsres, 2, sd), lwd=3, col = rgb(0,0,0.75,0.25))
# lines(A_test, simres$eb3$ests_loess - 1.96 * apply(bsres, 2, sd), lwd=3, col = rgb(0,0,0.75,0.25))
# lines(A_test, simres$eb3$ests_loess + 1.96 * apply(bsres, 2, sd), lty = 3, lwd=3, col = rgb(0,0,0,0.75))
# lines(A_test, simres$eb3$ests_loess - 1.96 * apply(bsres, 2, sd), lty = 3, lwd=3, col = rgb(0,0,0,0.75))
# # 
# plot(A_test, apply(bsstderr,2,mean)/apply(ests[4,,],2,sd),
#      type = 'l', lwd = 3,
#      col = rgb(0,0,0,0.5),
#      ylab = 'Ratio of S.E.',
#      xlab = 'Exposure Level',
#      main = 'Comparing Standard Errors',
#      axes = F)
# abline(h = 1, lty = 3)
# abline(v = Aquants, lty = c(2,3,3,2), col = rgb(0,0,0,0.5))
# axis(1); axis(2, las = 2); box();
# # 
# plot(A_test, apply(bscovers, 2, mean), ylim = c(0.8,1),
#      type = 'l', lwd = 3,
#      col = rgb(0,0,0,0.5),
#      ylab = 'Percent Coverage',
#      xlab = 'Exposure Level',
#      main = 'Confidence Interval Coverage',
#      axes = F)
# abline(h = 0.95, lty = 3)
# abline(v = Aquants, lty = c(2,3,3,2), col = rgb(0,0,0,0.5))
# axis(1); axis(2, las = 2); box();
# # 
# dev.off()




which_est <- 1


pdf('paper-figures/ne-bs-plot-one.pdf', height = 6, width = 12)
par(mfrow=c(2,2), oma = c(0,0,0,0), mar = c(4.25,4,2,1))
plot(A_test, bsdata$truth, type = 'l', lwd = 3, col = rgb(0.75,0,0,0.75), ylim = c(0, 10),
     ylab = 'Outcome', xlab = 'Exposure Level',
     main = 'Example Bootstrap Confidence Interval', axes = F)
axis(1); axis(2, las = 2); box();
abline(v = Aquants, lty = c(2,3,3,2), col = rgb(0,0,0,0.5))
lines(A_test, bsdata$ests[which_est,], lwd = 3, col = rgb(0,0,0.75,0.5))
lines(A_test, bsdata$ests[which_est,] + 2 * bsdata$`bs-se`[which_est,], lwd=3, col = rgb(0,0,0.75,0.25))
lines(A_test, bsdata$ests[which_est,] - 2 * bsdata$`bs-se`[which_est,], lwd=3, col = rgb(0,0,0.75,0.25))
lines(A_test, bsdata$ests[which_est,] + 2 * bsdata$`bs-se`[which_est,], lty = 3, lwd=3, col = rgb(0,0,0,0.75))
lines(A_test, bsdata$ests[which_est,] - 2 * bsdata$`bs-se`[which_est,], lty = 3, lwd=3, col = rgb(0,0,0,0.75))
legend('topleft', c('Truth', 'Estimated'), 
       col = c( rgb(0.75,0,0,0.75), rgb(0,0,0.75,0.5)), lwd = 3., cex = 1.25, bg='white')

# 
plot(A_test, apply(bsdata$`bs-cov`, 2, mean), ylim = c(0.5,1),
     type = 'l', lwd = 3,
     col = rgb(0,0,0,0.5),
     ylab = 'Percent Coverage',
     xlab = 'Exposure Level',
     main = 'Confidence Interval Coverage',
     axes = F)
abline(h = 0.95, lty = 3, col = 'red')
abline(h = c(0,1), lty = 3)
abline(v = Aquants, lty = c(2,3,3,2), col = rgb(0,0,0,0.5))
axis(1); axis(2, las = 2); box();

ymax <- max(c(apply(bsdata$`bs-se`,2,mean), apply(bsdata$ests,2,sd)))

plot(A_test, apply(bsdata$`bs-se`,2,mean),
     type = 'l', lwd = 3,
     ylim = c(0,ymax),
     col = rgb(0,0,0.75,0.5),
     ylab = 'Standard Error',
     xlab = 'Exposure Level',
     main = 'Comparing Standard Errors - Magnitude',
     axes = F)
lines(A_test, apply(bsdata$ests,2,sd), 
      lwd = 3, col = rgb(0.75,0,0,0.75))
abline(h = 1, lty = 3)
abline(v = Aquants, lty = c(2,3,3,2), col = rgb(0,0,0,0.5))
legend('topleft', c('Empirical SE', 'Bootstrap SE'), 
       col = c( rgb(0.75,0,0,0.75), rgb(0,0,0.75,0.5)), lwd = 3, cex = 1.25, bg='white')
axis(1); axis(2, las = 2); box();

# 
plot(A_test, apply(bsdata$`bs-se`,2,mean)/apply(bsdata$ests,2,sd),
     type = 'l', lwd = 3,
     col = rgb(0,0,0,0.5),
     ylab = 'Ratio of S.E.',
     xlab = 'Exposure Level',
     main = 'Comparing Standard Errors - Ratio',
     axes = F)
abline(h = 1, lty = 3)
abline(v = Aquants, lty = c(2,3,3,2), col = rgb(0,0,0,0.5))
axis(1); axis(2, las = 2); box()

dev.off()

