MX1 <- -0.5
MX2 <- 1
MX3 <- 0.3

amin <- 0
amax <- 45

Aquants <- c(0.93, 2.5, 30.5, 40.34)

A_test <- seq(amin, amax, length.out = 100)
truth <- - (A_test + 5) * (A_test - 5) / 300 + A_test * (2 + MX1^2 + MX2^2) / 25 + MX1 + MX2 + MX3