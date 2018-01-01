

rm(list = ls())


#### FSA yield curve

date_valuation <- "2012-08-17"

zcb_fsa <- read.table('data/zcb_yield_curve_danish_fsa.csv', header = T, sep = ';', dec = ',')
zcb_fsa$Dato <- as.Date(zcb_fsa$Dato, format = "%d-%m-%Y")
zcb_fsa <- as.numeric(zcb_fsa[which(zcb_fsa$Dato == date_valuation), -c(1, 2)] / 100)[-1]

time_fsa <- 0:(length(zcb_fsa) - 1)
n <- length(time_fsa)

zcb_fsa_fct <- approxfun(time_fsa, zcb_fsa)
plot(time_fsa, sapply(time_fsa, function(s) zcb_fsa_fct(s)))


#### FSA forward curve

fwd_rate <- function(t, s, r, time) {
  (r[s + 1] * time[s + 1] - r[t + 1] * time[t + 1]) / (time[s + 1] - time[t + 1])
}

fwd_rate_vector <- sapply(1:(n - 1), function(i) fwd_rate(time_fsa[i], time_fsa[i + 1], zcb_fsa, time_fsa))
fwd_rate_fct <- approxfun(c(-1, time_fsa), c(fwd_rate_vector[1], fwd_rate_vector, tail(fwd_rate_vector, 1)))
plot(time_fsa, sapply(time_fsa, function(s) fwd_rate_fct(s)), type = 'l', xlim = c(0, 25))



#### Estimate b1(t)

library(numDeriv)

# fix a and sigma with reasonable values (parameters as in 24.36 in Björk)
a <- 0.02
sigma <- 0.005

# determine theta function based on 24.47 in Björk
g <- function(t, a, sigma) sigma ^ 2 / (2 * a ^ 2) * (1 - exp(-a * t)) ^ 2
g_diff <- function(t, a, sigma) sigma ^ 2 / a * (1 - exp(-a * t)) * exp(-a * t)

b1_vector <- rep(NA, (n - 1))
for (i in 1:(n - 1)) {
  b1_vector[i] <- grad(function(x) fwd_rate_fct(x), time_fsa[i]) + g_diff(time_fsa[i], a, sigma) + a * (fwd_rate_vector[i] + g(time_fsa[i], a, sigma))
}



#### test that b1_vector function mathces zcb curve from FSA

B <- function(t, a) 1 / a * (1 - exp(-a * t))
A <- function(t, a, sigma, theta, time) {
  
  B_fct <- function(s, maturity) 1 / a * (1 - exp(-a * (maturity - s)))
  theta_fct <- approxfun(time, theta)
  integrate(function(s) 0.5 * sigma ^ 2 * B_fct(s, t) ^ 2 - theta_fct(s) * B_fct(s, t), 0, t, subdivisions = 200)$v
  
}

zcb_ext_vasicek <- sapply(time_fsa[1:(n - 1)], function(t) exp(A(t, a, sigma, b1_vector, time_fsa[1:(n - 1)]) - B(t, a) * fwd_rate_vector[1]) )
yield_ext_vasicek <- -log(zcb_ext_vasicek) / time_fsa[1:(n - 1)]

plot(time_fsa[1:(n - 1)], yield_ext_vasicek, type = 'l')
points(time_fsa[1:(n - 1)], zcb_fsa[1:(n - 1)])





#### solve 5.2.3 and 5.2.5

rho <- 0
beta1 <- function(t) c(0.02, 0)
beta2 <- function(t) c(0, 0.02)
sigma1 <- function(t) 0.005
sigma2 <- function(t) 0.15
base_surrender <- function(t) 0.06 - 0.002 * t
b1 <- splinefun(time_fsa, c(b1_vector, tail(b1_vector, 1)))
b2 <- function(t) 0.02

a <- function(t) {
  res <- matrix(0, 2, 2)
  res[1, 1] <- sigma1(t) ^ 2 * (1 - rho ^ 2)
  return(res)
}

alpha1 <- function(t) matrix(0, 2, 2)

alpha2 <- function(t) {
  res <- matrix(0, 2, 2)
  res[1, 1] <- sigma1(t) ^ 2 * rho ^ 2
  res[1, 2] <- res[2, 1] <- sigma1(t) * sigma2(t) * rho
  res[2, 2] <- sigma2(t) ^ 2
  return(res)
}



gamma1 <- function(t) 1
gamma2 <- base_surrender




dt <- 1/12
projectionYears <- 25
time <- seq(0, projectionYears, by = dt)
n <- length(time)

# plot(time, sapply(time, function(i) b1(i)))

result_phi <- rep(NA, n)
result_psi_1 <- rep(NA, n)
result_psi_2 <- rep(NA, n)

result_A_k1 <- rep(NA, n)
result_B_1_k1 <- rep(NA, n)
result_B_2_k1 <- rep(NA, n)

result_A_k2 <- rep(NA, n)
result_B_1_k2 <- rep(NA, n)
result_B_2_k2 <- rep(NA, n)




#### differential functions



diff_phi <- function(s, psi, a, b, c) {
  
  -1 / 2 * psi %*% a(s) %*% psi + b(s) %*% psi + sum(c)
  
}

diff_psi <- function(s, psi, alpha, beta, gamma) {
  
  1 / 2 * psi %*% alpha(s) %*% psi + beta(s) %*% psi - gamma(s)
  
}



diff_A <- function(s, B, psi, a, b) {
  
  psi %*% a(s) %*% B - b(s) %*% B
  
}

diff_B <- function(s, B, psi, alpha, beta) {
  
  psi %*% alpha(s) %*% B + beta(s) %*% B
  
}



f_phi <- function(s, psi) diff_phi(s, psi, a = a, b = function(s) c(b1(s), b2(s)), c = c(0, 0))
f_psi_1 <- function(s, psi) diff_psi(s, psi, alpha = alpha1, beta = beta1, gamma = gamma1)
f_psi_2 <- function(s, psi) diff_psi(s, psi, alpha = alpha2, beta = beta2, gamma = gamma2)

for (k in n:2) {

  # boundaries dependent on u

  u <- time[k]

  result_phi[n] <- 0
  result_psi_1[n] <- 0
  result_psi_2[n] <- 0

  result_A_k1[n] <- c(1, 0) %*% c(0, 0)
  result_B_1_k1[n] <- (c(1, 0) %*% matrix(c(1, 0, 0, base_surrender(u)), 2, 2))[1]
  result_B_2_k1[n] <- (c(1, 0) %*% matrix(c(1, 0, 0, base_surrender(u)), 2, 2))[2]

  result_A_k2[n] <- c(0, 1) %*% c(0, 0)
  result_B_1_k2[n] <- (c(0, 1) %*% matrix(c(1, 0, 0, base_surrender(u)), 2, 2))[1]
  result_B_2_k2[n] <- (c(0, 1) %*% matrix(c(1, 0, 0, base_surrender(u)), 2, 2))[2]


  for (i in n:(n - k + 2)){

    index_time <- i - (n - k)
    s <- time[index_time]

    ## 5.2.3

    phi_temp <- result_phi[i]
    psi_1_temp <- result_psi_1[i]
    psi_2_temp <- result_psi_2[i]
    psi_temp <- c(psi_1_temp, psi_2_temp)

    F1_phi <- f_phi(s, psi_temp) # psi is correct, since phi is not in diff_phi
    F1_psi_1 <- f_psi_1(s, psi_temp)
    F1_psi_2 <- f_psi_2(s, psi_temp)

    result_phi[i - 1] <- phi_temp - dt * F1_phi
    result_psi_1[i - 1] <- psi_1_temp - dt * F1_psi_1
    result_psi_2[i - 1] <- psi_2_temp - dt * F1_psi_2

    ## 5.2.5

    f_A <- function(s, B) diff_A(s, B, psi = psi_temp, a = a, b = function(s) c(b1(s), b2(s)))
    f_B_1 <- function(s, B) diff_B(s, B, psi = psi_temp, alpha = alpha1, beta = beta1)
    f_B_2 <- function(s, B) diff_B(s, B, psi = psi_temp, alpha = alpha2, beta = beta2)

    A_k1_temp <- result_A_k1[i]
    B_1_k1_temp <- result_B_1_k1[i]
    B_2_k1_temp <- result_B_2_k1[i]
    B_k1_temp <- c(B_1_k1_temp, B_2_k1_temp)

    A_k2_temp <- result_A_k2[i]
    B_1_k2_temp <- result_B_1_k2[i]
    B_2_k2_temp <- result_B_2_k2[i]
    B_k2_temp <- c(B_1_k2_temp, B_2_k2_temp)

    F1_A_k1 <- f_A(s, B_k1_temp) # B is correct, since A is not in diff_A
    F1_B_1_k1 <- f_B_1(s, B_k1_temp)
    F1_B_2_k1 <- f_B_2(s, B_k1_temp)

    F1_A_k2 <- f_A(s, B_k2_temp) # B is correct, since A is not in diff_A
    F1_B_1_k2 <- f_B_1(s, B_k2_temp)
    F1_B_2_k2 <- f_B_2(s, B_k2_temp)

    result_A_k1[i - 1] <- A_k1_temp - dt * F1_A_k1
    result_B_1_k1[i - 1] <- B_1_k1_temp - dt * F1_B_1_k1
    result_B_2_k1[i - 1] <- B_2_k1_temp - dt * F1_B_2_k1

    result_A_k2[i - 1] <- A_k2_temp - dt * F1_A_k2
    result_B_1_k2[i - 1] <- B_1_k2_temp - dt * F1_B_1_k2
    result_B_2_k2[i - 1] <- B_2_k2_temp - dt * F1_B_2_k2

  }

  if ( i != n ) {

    index <- i:n

    result_phi[index] <- NA
    result_psi_1[index] <- NA
    result_psi_2[index] <- NA

    result_A_k1[index] <- NA
    result_B_1_k1[index] <- NA
    result_B_2_k1[index] <- NA

    result_A_k2[index] <- NA
    result_B_1_k2[index] <- NA
    result_B_2_k2[index] <- NA

  }


}


result_phi <- rev(result_phi)
result_psi_1 <- rev(result_psi_1)
result_psi_2 <- rev(result_psi_2)
result_A_k1 <- rev(result_A_k1)
result_B_1_k1 <- rev(result_B_1_k1)
result_B_2_k1 <- rev(result_B_2_k1)
result_A_k2 <- rev(result_A_k2)
result_B_1_k2 <- rev(result_B_1_k2)
result_B_2_k2 <- rev(result_B_2_k2)



#### figure 5.3: dependent forward rates

x_0 <- c(fwd_rate_fct(0), 1)

A_r <- approxfun(time, result_A_k1)
B_1_r <- approxfun(time, result_B_1_k1)
B_2_r <- approxfun(time, result_B_2_k1)
fwd_r <- function(s) A_r(s) + x_0 %*% c(B_1_r(s), B_2_r(s))

A_n <- approxfun(time, result_A_k2)
B_1_n <- approxfun(time, result_B_1_k2)
B_2_n <- approxfun(time, result_B_2_k2)
fwd_n <- function(s) A_n(s) + x_0 %*% c(B_1_n(s), B_2_n(s))

plot(time, sapply(time, function(s) fwd_r(s)), type = 'l')
points(time, sapply(time, function(s) fwd_rate_fct(s)))

plot(time, sapply(time, function(s) fwd_n(s)), type = 'l')
points(time, sapply(time, function(s) base_surrender(s)))





