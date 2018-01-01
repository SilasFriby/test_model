

rm(list = ls())




############  DETERMINISTIC VS STOCHASTIC MORTALITY IN SURVIVAL MODEL ############ 




#### initialize

age <- 40
sex <- 1
projection_years <- 50
dt <- 1 / 12
time <- seq(0, projection_years, by = dt)
n_time <- length(time)







#### fsa mortality model

benchmarkFT <- read.table("data/Benchmark_doedelighed2012.csv", header = T, sep = ';', dec = '.')
benchmarkR  <- read.table("data/Benchmark_levetidsforbedringer2012.csv", header = T, sep = ';', dec = '.')

femaleR  <- approxfun(benchmarkR$Alder, benchmarkR$Kvinder)
femaleMu <- approxfun(benchmarkFT$Alder, benchmarkFT$Kvinder)
maleR  <- approxfun(benchmarkR$Alder, benchmarkR$Maend)
maleMu <- approxfun(benchmarkFT$Alder, benchmarkFT$Maend)

mu_fsa <- function(age, time, sex) {
  
  if (sex == 1) {
    
    res <- maleMu(age + time) * (1 - maleR(age + time)) ^ time
    
  } else if (sex == 0) {
    
    res <- femaleMu(age + time) * (1 - femaleR(age + time)) ^ time
    
  }  
  
  return(res)
  
}

mu_ad <- approxfun(age + time, mu_fsa(age, time, sex))



#### stochastic mortality model: dmu(t) = (b(t) - beta * mu) * dt + sigma * dW(t)


beta <- function(t) 0.5
sigma <- function(t) 0.001

a <- function(t) sigma(t) ^ 2
alpha <- function(t) 0
gamma <- function(t) 1
c <- 0




## Calibrate b(t) to fsa mortality

library(numDeriv)

# determine b function based on 24.47 in Björk
g <- function(t, beta, sigma) sigma(t) ^ 2 / (2 * beta(t) ^ 2) * (1 - exp(-beta(t) * t)) ^ 2
g_diff <- function(t, beta, sigma) sigma(t) ^ 2 / beta(t) * (1 - exp(-beta(t) * t)) * exp(-beta(t) * t)

mu_ad <- approxfun(c(age - 1, age + time), c(mu_ad(age), sapply(time, function(i) mu_ad(age + i)))) # for grad function from numDeriv to work
b_vector <- rep(NA, (n_time - 1))
for (i in 1:(n_time - 1)) {
  b_vector[i] <- grad(function(x) mu_ad(x), age + time[i]) + g_diff(time[i], beta, sigma) + beta(time(i)) * (mu_ad(age + time[i]) + g(time[i], beta, sigma))
}

b <- approxfun(time, c(b_vector, tail(b_vector, 1)))



#### solve 5.2.3 and 5.2.5

result_phi <- rep(NA, n_time)
result_psi <- rep(NA, n_time)
result_A <- rep(NA, n_time)
result_B <- rep(NA, n_time)



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


f_phi <- function(s, psi) diff_phi(s, psi, a = a, b = b, c = c)
f_psi <- function(s, psi) diff_psi(s, psi, alpha = alpha, beta = beta, gamma = gamma)

for (k in n_time:2) {
  
  # boundaries dependent on u
  
  u <- time[k]
  
  result_phi[n_time] <- 0
  result_psi[n_time] <- 0
  
  result_A[n_time] <- 0
  result_B[n_time] <- gamma(u)
  
  
  for (i in n_time:(n_time - k + 2)){
    
    index_time <- i - (n_time - k)
    s <- time[index_time]
    
    ## 5.2.3
    
    phi_temp <- result_phi[i]
    psi_temp <- result_psi[i]
    
    F1_phi <- f_phi(s, psi_temp) # psi is correct, since phi is not in diff_phi
    F1_psi <- f_psi(s, psi_temp)
    
    result_phi[i - 1] <- phi_temp - dt * F1_phi
    result_psi[i - 1] <- psi_temp - dt * F1_psi
    
    ## 5.2.5
    
    f_A <- function(s, B) diff_A(s, B, psi = psi_temp, a = a, b = b)
    f_B <- function(s, B) diff_B(s, B, psi = psi_temp, alpha = alpha, beta = beta)
    
    A_temp <- result_A[i]
    B_temp <- result_B[i]
    
    F1_A <- f_A(s, B_temp) # B is correct, since A is not in diff_A
    F1_B <- f_B(s, B_temp)
    
    result_A[i - 1] <- A_temp - dt * F1_A
    result_B[i - 1] <- B_temp - dt * F1_B
    
    
  }
  
  if ( i != n_time ) {
    
    index <- i:n_time
    
    result_phi[index] <- NA
    result_psi[index] <- NA
    result_A[index] <- NA
    result_B[index] <- NA
    
  }
  
  
}


result_phi <- rev(result_phi)
result_psi <- rev(result_psi)
result_A <- rev(result_A)
result_B <- rev(result_B)

mu_0 <- mu_ad(age)

A <- approxfun(time, result_A)
B <- approxfun(time, result_B)
fwd_mort <- function(s) A(s) + mu_0 * B(s) 

plot(time, fwd_mort(time), type = 'l')
points(time, mu_ad(age + time))





#### simulate mortality model


n_sim <- 5
mu_paths <- matrix(NA, n_time, n_sim)
mu_paths[1, ] <- mu_0

for (k in 1:n_sim) {
  for (i in 2:n_time) {
    dmu <- (b(time[i]) - beta(time[k]) * mu_paths[i - 1, k]) * dt + sigma(time[k]) * sqrt(dt) * rnorm(1,0,1)
    mu_paths[i, k] <- mu_paths[i - 1, k] + dmu
  }
  
  
}







#### deterministic mortality vs stochastic mortality

plot(time, mu_ad(age + time), type = 'l')
for (i in 1:n_sim) {
  points(time, mu_paths[, i], type = 'l')
}  



