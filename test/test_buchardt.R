
# rm(list = ls())

library(testthat)
library(timeDate) # provides the function timeFirstDayInMonth

source('test/assumptionsModel_buchardt.R')
source('test/assumptionsIntensities_buchardt.R')
source('functions/thiele.R')
source('functions/rkSolve.R')
source("functions/kolmogorov.R")
source("functions/probabilities_with_FP.R")
source('functions/helpFunctions.R')


#### initialize

test_tolerance <- 10 ^ (-2) # accepted relative difference between result found here and buchardt's result

dateValuation <- as.Date('2013-09-02')
ageRetirement <- 65
probInitial <- c(1, rep(0, nStates - 1)) # active
sex <- 1
age <- 40 
dateBirth <- firstDateOfNextPeriod(addToYearsDate(dateValuation, -age), period = dtString)
dateRetirement <- addToYearsDate(dateBirth, ageRetirement)


## annual benefits

benefitDisability <- 10 ^ 5
benefitRetirement <- 10 ^ 5
premium <- 46409
pv_buchardt <- -72641

rate_tech <- 0.01


#### model dates 

dateFirst <- firstDateOfNextPeriod(dateValuation, period = dtString)
dateLast <- addToYearsDate(dateFirst, projectionYears) 
datesModel <- seq.Date(dateFirst, dateLast, by = dtString )
nDates <- length(datesModel)
time <- seq(0, (nDates - 1) * dt, by = dt)



#### fair premium

## 1. order intensities

intensityListTech <- intensityTechAssumptionFct(age = age, 
                                                sex = sex, 
                                                ageRetirement = ageRetirement,
                                                r = rate_tech)


## payment functions

premiumFct <- function(s,premium){
  # premium: the premium paid per year while in state 1
  prm    <- matrix(0, length(s), nStatesTech)
  prm[, 1] <- premium * (s < ageRetirement - age)
  return(prm)
}

benefitContFct <- function(s,bCont) {
  # bCont: vector with entry i indicating the 'continuous' benefits paid in state i
  b      <- matrix(0, length(s), nStatesTech)
  b[, 1] <- bCont[2] * (s >= ageRetirement - age) 
  b[, 2] <- bCont[2] * (s >= ageRetirement - age) + bCont[1] * (s < ageRetirement - age)
  return(b)
}

benefitDisFct <- function(s,bDis) {
  # bDis: matrix indicating the benefits paid once when jumping from i to j 
  b      <- matrix(0, nStatesTech * length(s), nStatesTech)
  b[seq(1, nStatesTech * length(s), by = nStatesTech), nStatesTech] <- bDis[1, nStatesTech] * (s < ageRetirement - age)
  b[seq(2, nStatesTech * length(s), by = nStatesTech), nStatesTech]  <- bDis[2, nStatesTech] * (s < ageRetirement - age)
  return(b)
}

# product specifications
bCont <- c(benefitDisability, benefitRetirement) 
bDis <- matrix(0, nStatesTech, nStatesTech) 




#### 1. order reserve at different times given active

fair_premium <- uniroot(function(p){
  rksolve(c(0, 0, 0),
          function(t, x) thiele(t, x, p, bCont, bDis, intensityListTech$mu, rate_tech),
          rev(time))$xHist[721, 1]
  },
  c(-1, 10^5)
)$r

reserve_tech <- rksolve(c(0, 0, 0),
                        function(t, x) thiele(t, x, premium, bCont, bDis, intensityListTech$mu, rate_tech),
                        rev(time))$xHist[, 1]
reserve_tech <- rev(reserve_tech)

reserve_tech_plus <- rksolve(c(0, 0, 0),
                             function(t, x) thiele(t, x, 0, bCont, bDis, intensityListTech$mu, rate_tech),
                             rev(time))$xHist[, 1]
reserve_tech_plus <- rev(reserve_tech_plus)

reserve_first_order <- rbind(reserve_tech, reserve_tech_plus)

#### test: find reserve_tech_plus given reserve_tech

pMatrix_tech <-probabilityProjectionFct_FP(time = time, 
                                           probInitial = c(1, 0, 0), 
                                           stateCombi = c("12", "13", "21", "23"), 
                                           nStates = 3,
                                           intensityListTech$mu, 
                                           w1Vector = rep(0, n), 
                                           states_FP = states_FP)$pMatrix


reserve_tech_plus_test <- sapply(1:n, function(i) {
  
  if (age + time[i] < age_retirement) {
    index_retirement <- which(age + time == age_retirement) - 1
    reserve_tech[i] + sum(1 / (1 + rate_tech) ^ (time[i:index_retirement] - time[i]) * pMatrix_tech[i:index_retirement] * premium) * dt
  } else {
    reserve_tech[i]
  }
  
})

plot(time, reserve_tech_plus_test, type = "l")
points(time, reserve_tech)
points(time, reserve_tech_plus)

#### free policy factors - rho = 1. order reserve / 1. order benefit reserve

w1Vector <- ifelse(reserve_first_order[2, ] != 0, reserve_first_order[1, ] / reserve_first_order[2, ], 0)
w1Vector[w1Vector[] < 0] <- 0
w1Vector[nDates] <- 1
plot(time, w1Vector)


#### market state probabilities and FP probabilities

## intensities

intensityList <- intensityAssumptionFct(age = age, 
                                        sex = sex, 
                                        ageRetirement = ageRetirement, 
                                        timeVector = time, 
                                        muFT = muFT)

## probabilities

pMatrix <- probabilityProjectionFct_FP(time, probInitial, stateCombi, nStates,
                                       intensityList$mu, w1Vector, states_FP)




#### cash flow

muFct <- intensityList$mu
ageVec <- age + time
cf <- pMatrix$pMatrix[, 1] * (benefitRetirement * (ageVec > ageRetirement) - premium * (ageVec <= ageRetirement) + sapply(time, function(i) muFct(i, 1, 7)) * reserve_first_order[1, ]) +
  pMatrix$pMatrix[, 2] * (benefitDisability * (ageVec <= ageRetirement) + benefitRetirement * (ageVec > ageRetirement)) +
  pMatrix$pMatrix_FP[, 1] * (benefitRetirement * (ageVec > ageRetirement) + sapply(time, function(i) muFct(i, 4, 8)) * reserve_first_order[2, ]) +
  pMatrix$pMatrix_FP[, 2] * (benefitDisability * (ageVec <= ageRetirement) + benefitRetirement * (ageVec > ageRetirement))

cfMonthly <- cf * dt



plot(time, cfMonthly, type = 'l')


#### present value of cash flow

discount_curve <- read.table('test/data/zcb_yield_curve_danish_fsa.csv', header = T, sep = ';', dec = ',')
discount_curve$Dato <- as.Date(discount_curve$Dato, format = "%d-%m-%Y")
discount_curve <- discount_curve[which(discount_curve$Dato == dateValuation), -c(1, 2)] / 100
discount_curve_fct <- approxfun(seq(0, 135), discount_curve[, -1])
pv <- sum( (cfMonthly)  / (1 + discount_curve_fct(time)) ^ time )

relative_diff <- abs( (pv - pv_buchardt) / pv_buchardt )

expect_that( 0, equals(relative_diff, tolerance  = test_tolerance) )




