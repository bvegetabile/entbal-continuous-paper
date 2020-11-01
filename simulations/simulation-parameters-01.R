MX1 <- -0.5
MX2 <- 1
MX3 <- 0.3

amin <- 1.5
amax <- 45

Aquants <- c(1.45, 3.6, 35, 45.5)

A_test <- seq(amin, amax, length.out = 100)
truth <- - 0.15 * A_test^2 + A_test * (2 + MX1^2 + MX2^2) - 15 #+ 5 * (1 + MX1^2 + 6 * MX1 + 9) + 15 * (1 + MX2^2 + 6 * MX2 + 9) + MX3
truth <- truth / 50
# truth <- - (A_test - 10) * (A_test - 10) / 5 + 5 * A_test * (2 + MX1^2 + MX2^2) / 1 - 15 * (MX1 + MX2) + MX3 + (1 + MX2^2 - 40 * MX2 + 400)