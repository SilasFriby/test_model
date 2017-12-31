




#### laplace integration med 5 nedstigende differencer ##



# laplace <- function(f, a, b, n = 5) {
#   # a: nedre graense i integralet
#   # b: oevre graense i integralet
#   # n: antal nedstigende differenser
#   # f: funktion
#   
#   laplaceConst <- c(41393, -23719, 22742, -14762, 5449, -863) / 60480
#   
#   sum(sapply(a:b, function(i) f(i) )) +
#     sum(sapply(0:n, function(i) laplaceConst[i + 1] * f(b + i) )) -
#     sum(sapply(0:n, function(i) laplaceConst[i + 1] * f(a + i) ))
#   
# }

## vectorized version of laplace

laplace <- function(f, a, b, n = 5) {
  # a: nedre graense i integralet
  # b: oevre graense i integralet
  # n: antal nedstigende differenser
  # f: funktion
  
  laplaceConst <- c(41393, -23719, 22742, -14762, 5449, -863) / 60480
  sapply(1:length(b), function(k) {
    sum( f(a:b[k]) ) +
      sum( laplaceConst[1:(n + 1)] * f(b[k] + 0:n) ) -
      sum( laplaceConst[1:(n + 1)] * f(a + 0:n) )
  })
  
}


## useful functions

v <- function(i, age) (1 + i) ^ (-age)

l <- function(age, age0, mu_AD) exp( -laplace(mu_AD, age0, age) )

#l1 <- function(age, age0, A, B, C) exp( -A * (age - age0) - 10 ^ (B - 10) / (C * log(10)) * (10 ^ (C * age) - 10 ^ (C * age0)) )

D <- function(i, age, age0, mu_AD) v(i, age) * l(age, age0, mu_AD)

#D1 <- function(i, age, age0, A, B, C) v(i, age) * l1(age, age0, A, B, C)

l_ai <- function(age, age0, mu_AI) exp( -laplace(mu_AI, age0, age) )

l_a <- function(age, age0, mu_AD, mu_AI) l(age, age0, mu_AD) * l_ai(age, age0, mu_AI) 

D_a <- function(i, age, age0, mu_AD, mu_AI) v(i, age) * l_a(age, age0, mu_AD, mu_AI)

N <- function(i, age, age0, mu_AD) laplace(function(x) D(i, x, age0, mu_AD), age, 120)

N_a <- function(i, age, age0, mu_AD, mu_AI) laplace(function(x) D_a(i, x, age0, mu_AD, mu_AI), age, 120)



#### AKTIV ####

aktiv <- function(i, age, age0, mu_AD, mu_AI, ageRetirement) {
  ( N_a(i, age, age0, mu_AD, mu_AI) - N_a(i, ageRetirement, age0, mu_AD, mu_AI) ) / D_a(i, age, age0, mu_AD, mu_AI) *
    (age < ageRetirement)
}

aktiv_uden_prm_fritagelse <- function(i, age, age0, mu_AD, ageRetirement) {
  ( N(i, age, age0, mu_AD) - N(i, ageRetirement, age0, mu_AD) ) / D(i, age, age0, mu_AD) *
    (age < ageRetirement)
}


#### PASSIV 125 ####

passiv_125 <- function(i, age, age0, mu_AD, ageRetirement) D(i, ageRetirement, age0, mu_AD) / D(i, age, age0, mu_AD) * (age < ageRetirement) 



#### PASSIV 211 ####

passiv_211 <- function(i, age, age0, mu_AD, ageRetirement) {
  
  N(i, ageRetirement, age0, mu_AD) /  D(i, age, age0, mu_AD) * (age < ageRetirement) +
    N(i, age, age0, mu_AD) /  D(i, age, age0, mu_AD) * (age >= ageRetirement) 
  
}


#### PASSIV 415 ####

passiv_415 <- function(i, age, age0, mu_AD, mu_AI, ageRetirement) {
  
  ((N(i, age, age0, mu_AD) - N(i, ageRetirement, age0, mu_AD)) / D(i, age, age0, mu_AD) - 
    (N_a(i, age, age0, mu_AD, mu_AI) - N_a(i, ageRetirement, age0, mu_AD, mu_AI)) / D_a(i, age, age0, mu_AD, mu_AI)) * 
    (age < ageRetirement)
  
}





# aktiv(i = 0,
#       age = 30,
#       age0 = 0,
#       mu_AD = intensityListTech$mu13,
#       mu_AI = intensityListTech$mu12,
#       ageRetirement = ageRetirement)


# passiv_125(i = 0,
#            age = 30,
#            n = 35,
#            age0 = 0,
#            mu_AD = intensityListTech$mu13)


# passiv_415(i = 0, 
#            age = 30, 
#            age0 = 0, 
#            mu_AD = intensityListTech$mu13, 
#            mu_AI = intensityListTech$mu12, 
#            ageRetirement = ageRetirement)
# 
# 
# 
# D(i = 0, 
#   age = 30, 
#   age0 = 0, 
#   mu_AD = intensityListTech$mu13)
# 
# f <- function(x) D(i = 0, x, age0 = 0, mu_AD = intensityListTech$mu13)
# g <- function(x) D_a(i = 0, x, age0 = 0, mu_AD = intensityListTech$mu13, mu_AI = intensityListTech$mu12)
# laplace(a = 30, b = 120, n = 5, f = f)
# laplace(a = 30, b = 120, n = 5, f = g)
# 
# 
# intensityListTech$mu13(30)
# l(age = 30, age0 = 0, mu_AD = intensityListTech$mu13)
# D(i = 0, age = 30, age0 = 0, mu_AD = intensityListTech$mu13)
# D_a(i = 0, age = 30, age0 = 0, mu_AD = intensityListTech$mu13, mu_AI = intensityListTech$mu12)
# N(i = 0, age = 30, age0 = 0, mu_AD = intensityListTech$mu13)
# N_a(i = 0, age = 30, age0 = 0, mu_AD = intensityListTech$mu13, mu_AI = intensityListTech$mu12)



