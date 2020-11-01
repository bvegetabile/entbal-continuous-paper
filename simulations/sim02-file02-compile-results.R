source('simulation-functions/load_all_funcs_and_libs.R')

# ------------------------------------------------------------------------
# Simulation parameters ------
source('simulations/simulation-parameters-02.R')
# ------------------------------------------------------------------------

# Reading in results ---------------

simout <- readRDS('output/sim2-output1-effect-estimation.RDS')

truth <- rep(MX1 + (MX1^2 + 1) + MX2 + (MX2^2 + 1) + MX1 * MX2 + MX3, 100)

ests_loess <- simplify2array(lapply(simout, function(x) cbind(x$unwtd$ests_loess,
                                                              x$eb1$ests_loess,
                                                              x$eb2$ests_loess,
                                                              x$eb3$ests_loess,
                                                              x$eb4$ests_loess,
                                                              x$lm$ests_loess,
                                                              x$cbps$ests_loess,
                                                              x$npcbps1$ests_loess,
                                                              x$npcbps2$ests_loess,
                                                              x$npcbps3$ests_loess,
                                                              x$npcbps4$ests_loess,
                                                              x$gbm$ests_loess)))
ests_lmmod1 <- simplify2array(lapply(simout, function(x) cbind(x$unwtd$ests_lmmod1,
                                                               x$eb1$ests_lmmod1,
                                                               x$eb2$ests_lmmod1,
                                                               x$eb3$ests_lmmod1,
                                                               x$eb4$ests_lmmod1,
                                                               x$lm$ests_lmmod1,
                                                               x$cbps$ests_lmmod1,
                                                               x$npcbps1$ests_lmmod1,
                                                               x$npcbps2$ests_lmmod1,
                                                               x$npcbps3$ests_lmmod1,
                                                               x$npcbps4$ests_lmmod1,
                                                               x$gbm$ests_lmmod1)))

ests_lmmod2 <- simplify2array(lapply(simout, function(x) cbind(x$unwtd$ests_lmmod2,
                                                               x$eb1$ests_lmmod2,
                                                               x$eb2$ests_lmmod2,
                                                               x$eb3$ests_lmmod2,
                                                               x$eb4$ests_lmmod2,
                                                               x$lm$ests_lmmod2,
                                                               x$cbps$ests_lmmod2,
                                                               x$npcbps1$ests_lmmod2,
                                                               x$npcbps2$ests_lmmod2,
                                                               x$npcbps3$ests_lmmod2,
                                                               x$npcbps4$ests_lmmod2,
                                                               x$gbm$ests_lmmod2)))

ests_splmod <- simplify2array(lapply(simout, function(x) cbind(x$unwtd$ests_splmod,
                                                               x$eb1$ests_splmod,
                                                               x$eb2$ests_splmod,
                                                               x$eb3$ests_splmod,
                                                               x$eb4$ests_splmod,
                                                               x$lm$ests_splmod,
                                                               x$cbps$ests_splmod,
                                                               x$npcbps1$ests_splmod,
                                                               x$npcbps2$ests_splmod,
                                                               x$npcbps3$ests_splmod,
                                                               x$npcbps4$ests_splmod,
                                                               x$gbm$ests_splmod)))


cors <- simplify2array(lapply(simout, function(x) cbind(x$unwtd$cors,
                                                        x$eb1$cors,
                                                        x$eb2$cors,
                                                        x$eb3$cors,
                                                        x$eb4$cors,
                                                        x$lm$cors,
                                                        x$cbps$cors,
                                                        x$npcbps1$cors,
                                                        x$npcbps2$cors,
                                                        x$npcbps3$cors,
                                                        x$npcbps4$cors,
                                                        x$gbm$cors)))

ksstats <- simplify2array(lapply(simout, function(x) cbind(x$unwtd$ksstats,
                                                           x$eb1$ksstats,
                                                           x$eb2$ksstats,
                                                           x$eb3$ksstats,
                                                           x$eb4$ksstats,
                                                           x$lm$ksstats,
                                                           x$cbps$ksstats,
                                                           x$npcbps1$ksstats,
                                                           x$npcbps2$ksstats,
                                                           x$npcbps3$ksstats,
                                                           x$npcbps4$ksstats,
                                                           x$gbm$ksstats)))

ess <- simplify2array(lapply(simout, function(x) cbind(x$unwtd$ess,
                                                       x$eb1$ess,
                                                       x$eb2$ess,
                                                       x$eb3$ess,
                                                       x$eb4$ess,
                                                       x$lm$ess,
                                                       x$cbps$ess,
                                                       x$npcbps1$ess,
                                                       x$npcbps2$ess,
                                                       x$npcbps3$ess,
                                                       x$npcbps4$ess,
                                                       x$gbm$ess)))


cortab_mean <- t(apply(abs(cors),2,rowMeans))
cortab_max <- t(apply(abs(cors),2,function(x) apply(x, 1, max)))
kstab_mean <- t(apply(abs(ksstats),2,rowMeans))
kstab_max <- t(apply(abs(ksstats),2,function(x) apply(x, 1, max)))
ess_tab <- t(apply(abs(ess),2,rowMeans))
ess_max <- t(apply(abs(ess),2,function(x) apply(x, 1, max)))

bias1 <- apply(apply(ests_lmmod1, 2, function(x) apply(x, 2, function(x) x - truth)),2,mean)
bias2 <- apply(apply(ests_lmmod2, 2, function(x) apply(x, 2, function(x) x - truth)),2,mean)
bias3 <- apply(apply(ests_loess, 2, function(x) apply(x, 2, function(x) x - truth)),2,mean)
bias4 <- apply(apply(ests_splmod, 2, function(x) apply(x, 2, function(x) x - truth)),2,mean)

bias_quant1 <- apply(apply(ests_lmmod1, 2, function(x) apply(x, 2, function(x) x - truth)),2,quantile, c(0.05,0.25, 0.75,0.95))
bias_quant2 <- apply(apply(ests_lmmod2, 2, function(x) apply(x, 2, function(x) x - truth)),2,quantile, c(0.05,0.25, 0.75,0.95))
bias_quant3 <- apply(apply(ests_loess, 2, function(x) apply(x, 2, function(x) x - truth)),2,quantile, c(0.05,0.25, 0.75,0.95))
bias_quant4 <- apply(apply(ests_splmod, 2, function(x) apply(x, 2, function(x) x - truth)),2,quantile, c(0.05,0.25, 0.75,0.95))

absbias1 <- apply(apply(ests_lmmod1, 2, function(x) apply(x, 2, function(x) abs(x - truth))),2,mean)
absbias2 <- apply(apply(ests_lmmod2, 2, function(x) apply(x, 2, function(x) abs(x - truth))),2,mean)
absbias3 <- apply(apply(ests_loess, 2, function(x) apply(x, 2, function(x) abs(x - truth))),2,mean)
absbias4 <- apply(apply(ests_splmod, 2, function(x) apply(x, 2, function(x) abs(x - truth))),2,mean)

mse1 <- apply(apply(ests_lmmod1, 2, function(x) apply(x, 2, function(x) mse(x, truth))),2,mean)
mse2 <- apply(apply(ests_lmmod2, 2, function(x) apply(x, 2, function(x) mse(x, truth))),2,mean)
mse3 <- apply(apply(ests_loess, 2, function(x) apply(x, 2, function(x) mse(x, truth))),2,mean)
mse4 <- apply(apply(ests_splmod, 2, function(x) apply(x, 2, function(x) mse(x, truth))),2,mean)

msetab <- cbind(mse1, mse2, mse3, mse4)

set.seed(20200108)

out <- simout[[1]]
A50 <- seq(amin, 65, length.out = 100)
truth50 <- rep(MX1 + (MX1^2 + 1) + MX2 + (MX2^2 + 1) + MX1 * MX2 + MX3, 100)
aspan <- 0.45
smther <- loess(Y ~ A, data = out$dat, span = aspan)

# Visualizations --------------------------------------------------

pdf('paper-figures/ne-fig1.pdf', height = 4, width = 9)
par(mfrow=c(1,2), oma = c(0,0,0,0))
hist(out$dat$A,
     breaks = seq(0, 65, 2.5),
     xlab = 'Treatment',
     col = rgb(0,0,0,0.25),
     main = '',
     axes = F)
axis(1); axis(2, las = 2);
abline(v = Aquants, lty = c(2,3,3,2))
legend('topright', c('5th & 95th Quantile', '1st & 99th Quantile'), 
       lty = c(3, 2), bg = 'white', box.col = 'white')
plot(out$dat$A, out$dat$Y,
     pch = 19, col = rgb(0,0,0,0.25),
     xlim = c(0,65),
     xlab = 'Treatment',
     ylab = 'Outcome',
     axes = F)
axis(1); axis(2, las = 2);
lines(A50, predict(lm(Y ~ poly(A,2), data = out$dat), newdata = data.frame(A=A50)),
      col = 'green', lwd = 3)
lines(A50, predict(smther, newdata = data.frame(A=A50)),
      col = 'blue', lwd = 3)
lines(A50, truth50, lwd = 3, col = 'red')
dev.off()


pdf('paper-figures/ne-eb-response-loess.pdf', height = 6, width = 12)
par(mfrow=c(3,4), oma = c(0,0,0,0), mar = c(4,4,1,1)+0.1)
plot_sim(t(ests_loess[,1,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Unweighted')
plot_sim(t(ests_loess[,6,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Linear Model')
plot_sim(t(ests_loess[,7,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'CBPS - Parametric')
plot_sim(t(ests_loess[,12,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'GBM')

plot_sim(t(ests_loess[,2,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (1)')
plot_sim(t(ests_loess[,3,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (2)')
plot_sim(t(ests_loess[,4,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (3)')
plot_sim(t(ests_loess[,5,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (4)')

plot_sim(t(ests_loess[,8,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (1)')
plot_sim(t(ests_loess[,9,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (2)')
plot_sim(t(ests_loess[,10,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (3)')
plot_sim(t(ests_loess[,11,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (4)')
dev.off()


pdf('paper-figures/ne-eb-response-lmmod1.pdf', height = 6, width = 12)
par(mfrow=c(3,4), oma = c(0,0,0,0), mar = c(4,4,1,1)+0.1)
plot_sim(t(ests_lmmod1[,1,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Unweighted')
plot_sim(t(ests_lmmod1[,6,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Linear Model')
plot_sim(t(ests_lmmod1[,7,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'CBPS - Parametric')
plot_sim(t(ests_lmmod1[,12,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'GBM')

plot_sim(t(ests_lmmod1[,2,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (1)')
plot_sim(t(ests_lmmod1[,3,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (2)')
plot_sim(t(ests_lmmod1[,4,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (3)')
plot_sim(t(ests_lmmod1[,5,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (4)')

plot_sim(t(ests_lmmod1[,8,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (1)')
plot_sim(t(ests_lmmod1[,9,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (2)')
plot_sim(t(ests_lmmod1[,10,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (3)')
plot_sim(t(ests_lmmod1[,11,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (4)')
dev.off()


pdf('paper-figures/ne-eb-response-lmmod2.pdf', height = 6, width = 12)
par(mfrow=c(3,4), oma = c(0,0,0,0), mar = c(4,4,1,1)+0.1)
plot_sim(t(ests_lmmod2[,1,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Unweighted')
plot_sim(t(ests_lmmod2[,6,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Linear Model')
plot_sim(t(ests_lmmod2[,7,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'CBPS - Parametric')
plot_sim(t(ests_lmmod2[,12,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'GBM')

plot_sim(t(ests_lmmod2[,2,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (1)')
plot_sim(t(ests_lmmod2[,3,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (2)')
plot_sim(t(ests_lmmod2[,4,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (3)')
plot_sim(t(ests_lmmod2[,5,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (4)')

plot_sim(t(ests_lmmod2[,8,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (1)')
plot_sim(t(ests_lmmod2[,9,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (2)')
plot_sim(t(ests_lmmod2[,10,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (3)')
plot_sim(t(ests_lmmod2[,11,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (4)')
dev.off()


pdf('paper-figures/ne-eb-response-splmod.pdf', height = 6, width = 12)
par(mfrow=c(3,4), oma = c(0,0,0,0), mar = c(4,4,1,1)+0.1)
plot_sim(t(ests_splmod[,1,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Unweighted')
plot_sim(t(ests_splmod[,6,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Linear Model')
plot_sim(t(ests_splmod[,7,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'CBPS - Parametric')
plot_sim(t(ests_splmod[,12,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'GBM')

plot_sim(t(ests_splmod[,2,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (1)')
plot_sim(t(ests_splmod[,3,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (2)')
plot_sim(t(ests_splmod[,4,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (3)')
plot_sim(t(ests_splmod[,5,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (4)')

plot_sim(t(ests_splmod[,8,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (1)')
plot_sim(t(ests_splmod[,9,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (2)')
plot_sim(t(ests_splmod[,10,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (3)')
plot_sim(t(ests_splmod[,11,]), atest = A_test, truth = truth, aquants = Aquants, 
         ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (4)')
dev.off()



# 
# pdf('paper-figures/eb-response-lmmod1.pdf', height = 4, width = 12)
# par(mfrow=c(2,5), oma = c(0,0,0,0), mar = c(4,4,2,2)+0.1)
# plot_sim(t(ests_lmmod1[,1,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'Unweighted')
# plot_sim(t(ests_lmmod1[,2,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (1)')
# plot_sim(t(ests_lmmod1[,3,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (2)')
# plot_sim(t(ests_lmmod1[,4,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (3)')
# plot_sim(t(ests_lmmod1[,5,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (4)')
# 
# plot_sim(t(ests_lmmod1[,6,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'Linear Model')
# plot_sim(t(ests_lmmod1[,8,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (1)')
# plot_sim(t(ests_lmmod1[,9,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (2)')
# plot_sim(t(ests_lmmod1[,10,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (3)')
# plot_sim(t(ests_lmmod1[,11,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (4)')
# dev.off()
# 
# pdf('paper-figures/eb-response-lmmod2.pdf', height = 4, width = 12)
# par(mfrow=c(2,5), oma = c(0,0,0,0), mar = c(4,4,2,2)+0.1)
# plot_sim(t(ests_lmmod2[,1,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'Unweighted')
# plot_sim(t(ests_lmmod2[,2,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (1)')
# plot_sim(t(ests_lmmod2[,3,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (2)')
# plot_sim(t(ests_lmmod2[,4,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (3)')
# plot_sim(t(ests_lmmod2[,5,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (4)')
# 
# plot_sim(t(ests_lmmod2[,6,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'Linear Model')
# plot_sim(t(ests_lmmod2[,8,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (1)')
# plot_sim(t(ests_lmmod2[,9,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (2)')
# plot_sim(t(ests_lmmod2[,10,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (3)')
# plot_sim(t(ests_lmmod2[,11,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (4)')
# dev.off()
# 
# 
# pdf('paper-figures/eb-response-splmod.pdf', height = 4, width = 12)
# par(mfrow=c(2,5), oma = c(0,0,0,0), mar = c(4,4,2,2)+0.1)
# plot_sim(t(ests_splmod[,1,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'Unweighted')
# plot_sim(t(ests_splmod[,2,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (1)')
# plot_sim(t(ests_splmod[,3,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (2)')
# plot_sim(t(ests_splmod[,4,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (3)')
# plot_sim(t(ests_splmod[,5,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'Entropy Balancing: (4)')
# 
# plot_sim(t(ests_splmod[,6,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'Linear Model')
# plot_sim(t(ests_splmod[,8,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (1)')
# plot_sim(t(ests_splmod[,9,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (2)')
# plot_sim(t(ests_splmod[,10,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (3)')
# plot_sim(t(ests_splmod[,11,]), atest = A_test, truth = truth, aquants = Aquants, 
#          ylimits = c(0,10), ylines = 2, maintitle = 'CBPS Nonparametric: (4)')
# dev.off()
# 


# Tables -------------------------

baltab1 <- cbind(cortab_mean, kstab_mean)
baltab2 <- cbind(cortab_max, kstab_max)
rownames(baltab2) <- rownames(baltab1) <- c('Unweighted',
                                            'Entropy Balancing (1)',
                                            'Entropy Balancing (2)',
                                            'Entropy Balancing (3)',
                                            'Entropy Balancing (4)',
                                            'Linear Model',
                                            'CBPS - Parametric',
                                            'CBPS - Nonparametric: (1)',
                                            'CBPS - Nonparametric: (2)',
                                            'CBPS - Nonparametric: (3)',
                                            'CBPS - Nonparametric: (4)',
                                            'GBM - Optimal Weights')
xtable::xtable(baltab1, digits = 3)
xtable::xtable(baltab2, digits = 3)

knitr::kable(baltab1, digits = 2)
knitr::kable(baltab2, digits = 2)

# Estimates -----

restab <- cbind(t(ess_tab), bias1, bias2, bias3, bias4, msetab)
rownames(restab) <- c('Unweighted',
                      'Entropy Balancing (1)',
                      'Entropy Balancing (2)',
                      'Entropy Balancing (3)',
                      'Entropy Balancing (4)',
                      'Linear Model',
                      'CBPS - Parametric',
                      'CBPS - Nonparametric: (1)',
                      'CBPS - Nonparametric: (2)',
                      'CBPS - Nonparametric: (3)',
                      'CBPS - Nonparametric: (4)',
                      'GBM - Optimal Weights')
xtable::xtable(restab, digits = 3)
knitr::kable(restab, digits = 2)

# 
# 
# # Investigating
# par(mfrow=c(1,2))
# for(i in 1:10){
#   testsim <- simout[[i]]
# 
#   plot(testsim$dat$A, testsim$dat$Y,
#        pch = 19, col = rgb(0,0,0,testsim$weights$ebwts/max(testsim$weights$ebwts)))
#   plot(testsim$dat$A, testsim$dat$Y,
#        pch = 19, col = rgb(0,0,0,testsim$weights$npcbps2/max(testsim$weights$npcbps2)))
#   message(paste(rep('-',80), collapse = ''))
#   print(summary(lm(Y ~ A, data = testsim$dat, weights = testsim$weights$ebwts))$coef)
#   print(summary(lm(Y ~ A, data = testsim$dat, weights = testsim$weights$npcbps1))$coef)
#   print(summary(lm(Y ~ A, data = testsim$dat, weights = testsim$weights$npcbps2))$coef)
# 
# }
# 

headings <- c('Unweighted',
              'Entropy Balancing (1)',
              'Entropy Balancing (2)',
              'Entropy Balancing (3)',
              'Entropy Balancing (4)',
              'Linear Model',
              'CBPS - Parametric',
              'CBPS - Nonparametric: (1)',
              'CBPS - Nonparametric: (2)',
              'CBPS - Nonparametric: (3)',
              'CBPS - Nonparametric: (4)',
              'GBM')

pdf('paper-figures/ne-sim2-bias-dists.pdf', height = 6, width = 12)
par(mfrow=c(3,4), oma = c(0,0,0,0), mar = c(4,6,2,2)+0.1)

for(c in c(1,6,7,12,2:5,8:11)){
# }
# for(c in 1:length(bias1)){
#   if(c == 7 | c == 1 | c == 6) next 
  print(c)
  plot(0, 
       xlim = range(c(bias_quant1[,c(2:5,8:11)],bias_quant2[,c(2:5,8:11)],bias_quant3[,c(2:5,8:11)],bias_quant4[,c(2:5,8:11)])), 
       ylim = c(0.5, 4.5),
       main = headings[c],
       ylab = '',
       xlab = 'Bias',
       axes = F)
  abline(v = 0, lwd = 3, col = rgb(0,0,0,0.25))
  axis(1)
  axis(2, at = 4:1, c('B-Spline', 'LOESS', 'Quadratic', 'Linear Model'), las = 2)
  points(c(bias1[c], bias2[c], bias3[c], bias4[c]), c(1,2,3,4), pch = 19,)
  lines(bias_quant1[c(1,4),c], rep(1,2), col = rgb(0,0,0,0.5))
  lines(bias_quant2[c(1,4),c], rep(2,2), col = rgb(0,0,0,0.5))
  lines(bias_quant3[c(1,4),c], rep(3,2), col = rgb(0,0,0,0.5))
  lines(bias_quant4[c(1,4),c], rep(4,2), col = rgb(0,0,0,0.5))
  points(bias_quant1[,c], rep(1,4), pch = 4)
  points(bias_quant2[,c], rep(2,4), pch = 4)
  points(bias_quant3[,c], rep(3,4), pch = 4)
  points(bias_quant4[,c], rep(4,4), pch = 4)
  # lines()
}
dev.off()
