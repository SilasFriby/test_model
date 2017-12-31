
rm(list =ls())

library(rbenchmark)
library(Rcpp)
library(RcppArmadillo)

sourceCpp("rcpp/test_loop.cpp")

m <- 10^7
n <- 10
benchmark(
  f1(rep(1, m)),
  f2(rep(1, m)),
  f3(rep(1, m)),
  replications = n
  
)

