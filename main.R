
# rm(list = ls())

library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
library(rbenchmark)
library(pbapply)
library(testthat)

source('assumptions/assumptionsModel.R')
source('assumptions/mortality_fsa.R')

sourceCpp("rcpp/main.cpp")


#### initialize

m <- 1
age <- rep(40 ,m)
sex <- rep(1, m)
age_retirement <- rep(65, m)
state <- rep("active", m)
time <- seq(0, projectionYears, by = dt)
n <- length(time)
muAD_matrix <- pbsapply(1:m, function(i) approx(age[i] + time, muFT(age[i], time, sex[i]), xout = age[i] + time)$y)


## annual benefits

benefitDisability <- rep(10 ^ 5, m)
benefitRetirement <- rep(10 ^ 5, m)
premium <- rep(46409, m)
rate_tech <- 0.01



#### main function - cash flow



system.time(
  res <- main_cpp(time,
                  age,
                  sex, 
                  age_retirement, 
                  state, 
                  nStates,
                  premium, 
                  benefitDisability, 
                  benefitRetirement,
                  rate_tech, 
                  nStatesTech, 
                  dt, 
                  muAD_matrix, 
                  stateCombi, 
                  states_FP)
)

# n <- 10
# test_time <- benchmark(
#   res <- main_cpp(time,
#                   age,
#                   sex, 
#                   age_retirement, 
#                   state, 
#                   nStates,
#                   premium, 
#                   benefitDisability, 
#                   benefitRetirement,
#                   rate_tech, 
#                   nStatesTech, 
#                   dt, 
#                   muAD_matrix, 
#                   stateCombi, 
#                   states_FP),
#   replications = n
# )
# test_time$elapsed/n/m


#### present value of cash flow

discount_curve <- read.table('test/data/zcb_yield_curve_danish_fsa.csv', header = T, sep = ';', dec = ',')
discount_curve$Dato <- as.Date(discount_curve$Dato, format = "%d-%m-%Y")
discount_curve <- discount_curve[which(discount_curve$Dato == as.Date("2013-09-02")), -c(1, 2)] / 100
discount_curve_fct <- approxfun(seq(0, 135), discount_curve[, -1])
pv_cpp <- sum( (res$cf_matrix[1:720, 1])  / (1 + discount_curve_fct(time[1:720])) ^ time[1:720] )

pv_cpp


expect_that(0, equals(-75092.54 - pv_cpp, tolerance  = 10 ^ -2))
