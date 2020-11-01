library(parallel)
library(doParallel)

library(dplyr)
library(magrittr)

source('simulation-functions/load_all_funcs_and_libs.R')

source('/Users/bvegetab/Documents/GitHub/ebcontinuous/simulation-functions-params.R')
dat <- read.csv('/Users/bvegetab/Documents/data/imp1.csv')
sessions <- read.csv('/Users/bvegetab/Documents/data/12-17_sescnt_DID.csv')
dat <- merge(dat, sessions, by = 'xpid2')


dat = dat %>%
  mutate(sfs8p_0 = sfs8p_0,
         detox90_0 = case_when(s5a_0 > 0 ~ 1,
                               s5a_0 == 0 ~ 0,
                               TRUE ~ NA_real_),
         cws_0 = cws_0,
         nps_0 = nps_0,
         wphsu = (p6a_al_0 + p6b_to_0 + p6b_to_0 + p6b_mj_0 + p6b_cr_0 + p10b_an_0 +
                    p10c_an_0 + p10e_stm_0 + p10f_to_0 + p10f_mj_0 + p10f_cr_0 + p10h_stm_0 +
                    p10q_an_0 + p10r_an_0 + nsp_al_0 + nsp_stm_0 + nsp_sed_0),
         cis_0 = cis_0,
         dss9_0 = dss9_0,
         mhtrt_0 = case_when(mhu90_0 == 0 & meru90_0 == 0 ~ 0,
                             mhu90_0 > 0 & meru90_0 == 0 ~ 1,
                             meru90_0 > 0 ~ 2,
                             TRUE ~ NA_real_),
         hsts2_0 = hsts_0 - m1c1,
         gcts = gcts,
         csudptx_0 = csudptx_0,
         sati_0 = sati_0,
         sud_sy_0 = sud_sy_0,
         lri7_0 = lri7_0,
         sri7_0 = sri7_0,
         reri_0 = e14a_0 + e14b_0,
         pai_0 = pai_0,
         ada_0 = ada_0,
         recov_0 = recov_0,
         imds_0 = imds_0,
         bcs_0 = bcs_0,
         tss_0 = tss_0,
         eps7p_0 = eps7p_0,
         codis4 = codis4,
         school_0 = school_0,
         homeless_0 = homeless_0)

# Need six month data and only using ACRA data
dset <- dat[dat$has6m == 1 & dat$txtype == '1-ACRA',
            c('sncnt',
              'gender',
              'age',
              'race',
              'victim',
              'sfs8p_0',
              'detox90_0',
              'cws_0',
              'nps_0',
              'wphsu',
              'cis_0',
              'dss9_0',
              'mhtrt_0',
              'hsts2_0',
              'gcts',
              'csudptx_0',
              'sati_0',
              'sud_sy_0',
              'lri7_0',
              'sri7_0',
              'reri_0',
              'pai_0',
              'ada_0',
              'maxce_0',
              'recov_0',
              'imds_0',
              'bcs_0',
              'tss_0',
              'eps7p_0',
              'codis4_n',
              'eps7p_6',          # 31 Emotional Problem Scale at 6 Months
              'recov_6',
              'engage')]          # 31 Adjusted Days Abstinent at 6 Months

tmax <- 50
dset$eps7p_0 <- dset$eps7p_0 / 100
dset$eps7p_6 <- dset$eps7p_6 / 100
dset$recov_0 <- ifelse(dset$recov_0 > 0.5, 1, 0)
dset$recov_6 <- ifelse(dset$recov_6 > 0.5, 1, 0)
dset$mhtrt_0 <- ifelse(dset$mhtrt_0 > 0, 1, 0)
dset <- dset[dset$engage == 1 ,]

# ------------------------------------------------------------------------------
# Figure 1


pdf('paper-figures/sessions-withengagement-appfig1-1.pdf',
    height = 4.5, width = 9)
par(mfrow=c(2,3), oma = c(0,0,0,0), mar = c(5,4,1,3)+0.1)
hist(dset$sncnt, breaks = seq(0,80, by = 5),
     col = rgb(0,0,0,0.25),
     xlab = 'Number of Sessions',
     main = '', axes = F); axis(1); axis(2, las = 2);
hist(dset$eps7p_6, breaks = seq(0,1, by = .1),
     col = rgb(0,0,0,0.25),
     xlab = 'Emotional Problem Scale - 6mo', main = '', axes = F);
axis(1); axis(2, las = 2)

g0 <- table('Baseline' = dset$recov_0, 'Outcome' = dset$recov_6)

plot(0, xlim = c(0, 1), ylim = c(0,1700), pch = 19, col = rgb(0,0,0,0),
     axes = F, ylab = '', xlab = 'Summary of "In Recovery" Variable');
lines(c(0.25,0.25) - 0.05, c(0, sum(dset$recov_0 == 0)))
lines(c(0.25,0.25) + 0.05, c(0, sum(dset$recov_0 == 1)), lty = 3)
lines(c(0.75,0.75) - 0.05, c(0, sum(dset$recov_6 == 0)))
lines(c(0.75,0.75) + 0.05, c(0, sum(dset$recov_6 == 1)), lty = 3)
axis(1, at = c(0.25, 0.75), c('Baseline', '6-Month'))
points(c(0.2,0.3), rowSums(g0), pch = c(19,17), col = rgb(0,0,0,1), cex = 1.5)
points(c(0.7,0.8), colSums(g0), pch = c(19,17), col = rgb(0,0,0,1), cex = 1.5)
text(c(0.2,0.3), rowSums(g0), rowSums(g0), pos = 4)
text(c(0.7,0.8), colSums(g0) - 100, colSums(g0), pos = c(2,4))
axis(2, at = seq(0, 1500, 250), las = 2)
legend('topright',
       c('In Recovery', 'Not In Recovery'),
       pch = c(19,17), bty = 'n')

# hist(dset$recov_6, xlab = 'In Recovery: 6 Months', main = '',
#      col = rgb(0,0,0,0.25), axes = F);
# axis(1, at = c(0.05,0.95), labels = c('No', 'Yes')); axis(2, las = 2);
# text(x = c(0.1,0.9), y = table(dset$recov_6)-25,  labels = table(dset$recov_6), pos = c(4,2))


plot(0, xlab = '', ylab = '', main = '', col = rgb(0,0,0,0), axes = F,
     ylim = c(0,5.25), xlim = c(0,1))
text(x = 0, y = 4.9, labels = paste('Data Summaries:\nNumber of Sessions'), pos = 4)
text(x = 0, y = 4, labels = 'Mean:', pos = 4)
text(x = 0.5, y = 4, labels = round(mean(dset$sncnt),3), pos = 4)
text(x = 0, y = 3.5, labels = 'Median:', pos = 4)
text(x = 0.5, y = 3.5, labels = round(median(dset$sncnt),3), pos = 4)
text(x = 0, y = 3, labels = 'S.D.:', pos = 4)
text(x = 0.5, y = 3, labels = round(sd(dset$sncnt),3), pos = 4)

text(x = 0, y = 2, labels = 'Min:', pos = 4)
text(x = 0.5, y = 2, labels = round(min(dset$sncnt),3), pos = 4)
text(x = 0, y = 1.5, labels = '1st Quantile:', pos = 4)
text(x = 0.5, y = 1.5, labels = round(quantile(dset$sncnt, 0.01),3), pos = 4)
text(x = 0, y = 1, labels = '99th Quantile:', pos = 4)
text(x = 0.5, y = 1, labels = round(quantile(dset$sncnt, 0.99),3), pos = 4)
text(x = 0, y = 0.5, labels = 'Max:', pos = 4)
text(x = 0.5, y = 0.5, labels = round(max(dset$sncnt),3), pos = 4)
plot(dset$sncnt,
     dset$eps7p_6,
     xlim = c(0,50), ylim = c(0,1),
     pch = 19, col = rgb(0,0,0,0.15),
     xlab = 'Number of Sessions',
     ylab = 'Emotional Problem Scale - 6mo', axes = F)
axis(1, las = 2); axis(2, las =2 );

h0 <- hist(dset$sncnt[dset$recov_6==0], breaks = seq(0,90,5), plot = F)
h1 <- hist(dset$sncnt[dset$recov_6==1], breaks = seq(0,90,5), plot = F)

count_recov <- table(dset$recov_6, dset$sncnt)

plot(h0$mids[1:10],
     (h1$counts / (h0$counts + h1$counts))[1:10], type = 'b',
     xlim = c(0,tmax), ylim = c(0,1),
     pch = 19, col = rgb(0,0,0,0.5),
     xlab = 'Number of Sessions',
     ylab = 'Fraction In Recovery: 6 Months', axes = F); axis(1, las = 2); axis(2, las =2 );
abline(v = h1$breaks, lty = 3, col = rgb(0,0,0,0.4))
text(h0$mids[1:10], rep(0.1, length(h0$mids[1:10])), h0$counts[1:10], srt = 90, col = rgb(0,0,0,0.5))
text(h1$mids[1:10], rep(0.9, length(h1$mids[1:10])), h1$counts[1:10], srt = 90, col = rgb(0,0,0,0.5))
axis(4, at = c(0.1, 0.9), labels = c(expression(N[0], N[1])), las = 2)
abline(h = 0.5, lty = 1, lwd = 2, col = rgb(0.5,0,0, 0.25))
dev.off()

# ------------------------------------------------------------------------------
# Analysis

dset <- dset[dset$sncnt < tmax,]

XD <- dset[, 2:29]
XD$victim <- ifelse(XD$victim > 0.5, 1, 0)
XD$recov_0 <- ifelse(XD$recov_0 > 0.5, 1, 0)
XD$csudptx_0 <- ifelse(XD$csudptx_0 > 0.5, 1, 0)
XD$mhtrt_0 <- ifelse(XD$mhtrt_0 > 0, 1, 0)
XD$gender <- ifelse(XD$gender == 'Male', 1, 0)
XD$raceW <- ifelse(XD$race == 'White/Caucasian', 1, 0)
XD$raceB <- ifelse(XD$race == 'Black/African-American', 1, 0)
XD$raceH <- ifelse(XD$race == 'Hispanic', 1, 0)
XD$raceO <- ifelse(XD$race == 'Other', 1, 0)
XD$codis1 <- ifelse(dset$codis4_n == 1, 1, 0)
XD$codis2 <- ifelse(dset$codis4_n == 2, 1, 0)
XD$codis3 <- ifelse(dset$codis4_n == 3, 1, 0)
XD$codis4 <- ifelse(dset$codis4_n == 4, 1, 0)

XD <- XD[, -3]
A <- dset[, 1]
Y1 <- dset[, 31]
Y2 <- dset[,32]

CZ <- makeC2(XD, A, n_moments = 3)

ebz <- entbal_fit(CZ, rep(0,ncol(CZ)),
                  n_moments = 1,
                  verbose = T,
                  bal_tol = 1e-10, max_iters = 10000)

ess <- 1 / sum(ebz$wts^2)
print(c(nrow(XD), ess))

# Covariate Balance

wtd_ecdf <- function (var_data, wts) {
  #-----------------------------------------------------------------------------
  # wtd_ecdf is a modification of the ecdf() function in base R.  It modifies
  # the function to be able to incorporate weights.  This is to visualize
  # balance using the empirical cumulative distribution function for continuous
  # covariates after weighting by the inverse of the propensity score (IPTW)
  #
  # Input variables
  # --- var_data : covariate values - vector of data
  # --- wts      : weights for assessing cov balance by IPTW - vector of data.
  #-----------------------------------------------------------------------------
  ord <- order(var_data)
  var_ordered <- var_data[ord]
  wts_ordered <- wts[ord]
  n <- length(var_data)
  if (n < 1)
    stop("'var_data' must have 1 or more non-missing values")
  vals <- unique(var_ordered)
  matched_vals <- match(var_ordered, vals)
  weight_list <- aggregate(wts_ordered, by=list(matched_vals), sum)
  rval <- approxfun(vals, cumsum(weight_list[,2])/sum(wts_ordered),
                    method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}

ksbal <- function(X, wts, plotit=F) {
  xr <- range(X)
  xs <- seq(xr[1], xr[2], length.out = 500)
  b4dist <- ecdf(X)
  wdist <- wtd_ecdf(X, wts)
  if(plotit){
    plot(xs, b4dist(xs), type = 'l', col = rgb(1,0,0,0.5),
         xlim = xr, ylim = c(0,1),
         ylab = 'ECDF', xlab = 'covariate value',
         axes = F); axis(1); axis(2, las = 2);
    lines(xs, wdist(xs), col = rgb(0,0,1,0.5))
    abline(h=c(0,1), lty = 3, col = rgb(0,0,0,0.5))
  }
  max(abs(b4dist(xs) - wdist(xs)))
}


CB <- cbind(round(apply(XD, 2, mean), 2),
            round(apply(XD, 2, sd), 2),
            round(apply(XD, 2, function(x) cor(x,A)), 2),
            # t(apply(XD, 2, function(x) summary(lm(x ~ scale(A)))$coef[2,c(1,4)])),
            round(apply(XD, 2, function(x) wmean(ebz$wts, x)), 2),
            round(apply(XD, 2, function(x) sqrt(wvar(ebz$wts, x))), 2),
            round(apply(XD, 2, function(x) wcor(ebz$wts, x, A)), 2),
            # t(apply(XD, 2, function(x) lm_ps(x, A, ebz$wts)$ests[2,c(1,4)])),
            apply(XD, 2, function(x) ksbal(x, ebz$wts)))
CB <- rbind('Number of Sessions' = c(mean(A), sd(A), NA,  
                                     wmean(x= A, wts = ebz$wts),
                                     sqrt(wvar(wts = ebz$wts, x = A)), NA, 
                                     ksbal(A, ebz$wts, plotit = T)), CB)
colnames(CB) <- c('Mean', 'SD', 'Cor', 
                  'Mean', 'SD', 'Cor', 'KS')
# colnames(CB) <- c('Before', 'After')
knitr::kable(CB, digits = 3)
xtable::xtable(CB, digits = 2)


# Outcome Modeling -------------------------------------------------------------

# Model 1 ----------------------------------------------------------------------

set.seed(202001)

ss_spans0 <- function(x){
  loess_split_sample(Y1, A, XD, n_moments = 3, model = 'none', spans = seq(0.2, 1, by = 0.05))$best_span
}

ss_spans1 <- function(x){
  loess_split_sample(Y1, A, XD, n_moments = 3, model = 'EB', spans = seq(0.2, 1, by = 0.05))$best_span
}

options(cores = 8)
cl <- makeCluster(8)
registerDoParallel(cl)
system.time(cvspans <- foreach(i = 1:100,
                               .errorhandling = 'remove',
                               .verbose = T)
            %dopar% ss_spans1(i))
stopCluster(cl)
hist(unlist(cvspans))
bspan <- mean(unlist(cvspans))


options(cores = 12)
cl <- makeCluster(12)
registerDoParallel(cl)
system.time(cvspans0 <- foreach(i = 1:100,
                                .errorhandling = 'remove',
                                .verbose = T)
            %dopar% ss_spans0(i))
stopCluster(cl)
hist(unlist(cvspans0))
bspan0 <- mean(unlist(cvspans0))


# Bootstrap ests

bsest_example <- function(j, X, A, Y, testpts = seq(0, tmax, by = 1)) {
  NX <- ncol(X)
  dat <- data.frame(X=X,
                    A=A,
                    Y=Y)
  
  bsamp <- sample(1:nrow(dat), nrow(dat), replace = T)
  datbs <- dat[bsamp, ]
  
  C <- makeC2(datbs[,1:NX], datbs$A, n_moments = 3)
  eb <- entbal_fit(C, rep(0,ncol(C)),
                   n_moments = 1,
                   verbose = F,
                   bal_tol = 1e-8)
  datbs$W <- eb$wts
  #
  bspan <- loess_split_sample(datbs$Y, 
                              datbs$A, 
                              datbs[,1:NX], 
                              n_moments = 3, 
                              model = 'EB', 
                              spans = seq(0.2, 1, by = 0.05))$best_span
  
  # bspan <- 0.25
  
  loessmod <- loess(Y ~ A,
                    weights = W,
                    data = datbs,
                    degree = 1,
                    span = bspan,
                    control = loess.control(surface = 'direct'))
  
  pls <- predict(loessmod, newdata = data.frame(A = testpts))
  pls
}

options(cores = 12)
cl <- makeCluster(12)
registerDoParallel(cl)
system.time(bsests1 <- foreach(i = 1:100,
                               .errorhandling = 'remove',
                               .combine = 'rbind',
                               .verbose = T)
            %dopar% bsest_example(i, XD, A, Y1))
stopCluster(cl)
bspan <- mean(unlist(cvspans))

loessmod <- loess(Y1 ~ A, weights = ebz$wts, span = bspan, degree = 1)
lower <- predict(loessmod, newdata = data.frame(A = seq(0, tmax, by = 1))) - 1.96 * apply(bsests1, 2, sd)
upper <- predict(loessmod, newdata = data.frame(A = seq(0, tmax, by = 1))) + 1.96 * apply(bsests1, 2, sd)

# Other models
loessmod2 <- loess(Y1 ~ A, span = bspan0, degree = 1)
naive_lmmod <- lm(eps7p_6 ~ sncnt, data = dset)
contr_lmmod <- lm(eps7p_6 ~ sncnt + eps7p_0, data = dset)
wtd_lmmod <- lm(eps7p_6 ~ sncnt, data = dset, weights = ebz$wts)

treatvals <- seq(0, tmax, by = 0.5)

pdf('paper-figures/sessions-withengagement-appfig2-1.pdf', height = 4, width = 12)
par(mfrow=c(1,3))
plot(dset$sncnt,
     dset$eps7p_6,
     xlim = c(0,tmax), col = rgb(0,0,0,0.1), pch = 19,
     xlab = 'Number of Sessions',
     ylab = 'Emotional Problem Scale', axes = F)
lines(treatvals,
      predict(loessmod,
              newdata = data.frame(A = treatvals)),
      col = rgb(0,0,0.75,1), lwd = 4)
lines(seq(0, tmax, by = 1),
      lower,
      col = rgb(0,0,0.75,1), lwd = 2, lty = 2)
lines(seq(0, tmax, by = 1),
      upper,
      col = rgb(0,0,0.75,1), lwd = 2, lty = 2)
axis(1); axis(2, las = 2)
abline(h=0, col = rgb(1,0,0,0.5), lwd = 2)
plot(treatvals,
     predict(loessmod,
             newdata = data.frame(A = treatvals)),
     type = 'l',
     xlim = c(0,tmax), ylim = c(min(lower,na.rm=T),max(upper,na.rm=T)),
     col = rgb(0,0,0.75,1), lwd = 3,
     xlab = 'Number of Sessions',
     ylab = 'Emotional Problem Scale', axes = F)
axis(1); axis(2, las = 2)
abline(h = seq(0.15,0.3, 0.05), lty = 3 )
lines(seq(0, tmax, by = 1),
      lower,
      col = rgb(0,0,0.75,1), lwd = 2, lty = 2)
lines(seq(0, tmax, by = 1),
      upper,
      col = rgb(0,0,0.75,1), lwd = 2, lty = 2)
plot(treatvals,
     predict(loessmod,
             newdata = data.frame(A = treatvals)),
     type = 'l',
     xlim = c(0,tmax), ylim = c(min(lower,na.rm=T),max(upper,na.rm=T)),
     col = rgb(0,0,0.75,1), lwd = 3,
     xlab = 'Number of Sessions',
     ylab = 'EPS Difference', axes = F)
axis(1); axis(2, las = 2)
abline(h = seq(0.15,0.3, 0.05), lty = 3 )
lines(treatvals, predict(loessmod2,
                         newdata = data.frame(A = treatvals)),
      lwd = 3, col = rgb(0,0,0,0.5), lty = 2)
lines(treatvals,
      predict(naive_lmmod, newdata = data.frame(sncnt = treatvals)),
      col = rgb(0.75,0,0,0.5), lwd = 3)
lines(treatvals,
      predict(contr_lmmod, newdata = data.frame(sncnt = treatvals,
                                                eps7p_0 = mean(dset$eps7p_0))),
      col = rgb(0.75,0,0,0.5), lwd = 3, lty = 2)
lines(treatvals,
      predict(wtd_lmmod, newdata = data.frame(sncnt = treatvals)),
      col = rgb(0,0,0.75,0.5), lwd = 3, lty=2)
legend('topleft',
       c('EB + LOESS',
         'Naive LOESS',
         'Naive Reg.',
         'Baseline Control',
         'Wtd. Reg. - EB'),
       lwd = 3, lty = c(1,2,1,2,2),
       bg = 'white',
       col = c(rgb(0,0,0.75,0.75),
               rgb(0,0,0,0.5),
               rgb(0.75,0,0,0.5),
               rgb(0.75,0,0,0.5),
               rgb(0,0,0.75,0.5)),
       cex = 1.5)
title('\n\nRelationship between Number of Sessions and Emotional Problem Scale at 6-Months', outer = T)
dev.off()


# Model 2 ----------------------------------------------------------------------

set.seed(202001)

ss_spans2 <- function(x){
  loess_split_sample(Y2, A, XD, n_moments = 3, model = 'none', spans = seq(0.2, 1, by = 0.05))$best_span
}

ss_spans3 <- function(x){
  loess_split_sample(Y2, A, XD, n_moments = 3, model = 'EB', spans = seq(0.2, 1, by = 0.05))$best_span
}

