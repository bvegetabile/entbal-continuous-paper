MX1 <- -0.5
MX2 <- 1
MX3 <- 0.3

amin <- 1.5
amax <- 45

Aquants <- c(1.45, 3.6, 35, 45.5)

A_test <- seq(amin, amax, length.out = 100)
truth <- rep(MX1 + (MX1^2 + 1) + MX2 + (MX2^2 + 1) + MX1 * MX2 + MX3, 100)