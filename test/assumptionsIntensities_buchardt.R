

#########################################################################
############             INTENSITY ASSUMPTIONS               ############
#########################################################################


## FTs levetidsmodel

benchmarkFT <- read.table("test/data/Benchmark_doedelighed2012.csv", header = T, sep = ';', dec = ',')
benchmarkR  <- read.table("test/data/Benchmark_levetidsforbedringer2012.csv", header = T, sep = ';', dec = ',')
# head(benchmarkFT)
# head(benchmarkR)

femaleR  <- benchmarkR$Kvinder
femaleMu <- benchmarkFT$Kvinder
maleR    <- benchmarkR$Maend 
maleMu   <- benchmarkFT$Maend #* (1 - maleR) ^ (1 + 2/12)

# # FTs levetidsmodel
# muFT <- function(x,t,k) {
#   # x: alder
#   # t: ?r ud i fremtiden
#   # k: k?n
#   if (k == 'K') (1 - femaleR[x + 1 + t]) ^ t * femaleMu[x + 1 + t]
#   else if (k == 'M') (1 - maleR[x + 1 + t]) ^ t * maleMu[x + 1 + t]
# }

# FTs levetidsmodel
muFT <- function(x, t, k) {
  # x: alder
  # t: ?r ud i fremtiden
  # k: k?n
  
  lower <- floor(x + t)
  upper <- ceiling(x + t)
  
  if (k == 'K') {
    femaleMu_temp <- ((x + t) - lower) * femaleMu[upper + 1] + (upper - (x + t)) * femaleMu[lower + 1]
    femaleR_temp <- ((x + t) - lower) * femaleR[upper + 1] + (upper - (x + t)) * femaleR[lower + 1]
    res <- (1 - femaleR_temp) ^ t * femaleMu_temp
  }
  else if (k == 1) {
    maleMu_temp <- ((x + t) - lower) * maleMu[upper + 1] + (upper - (x + t)) * maleMu[lower + 1]
    maleR_temp <- ((x + t) - lower) * maleR[upper + 1] + (upper - (x + t)) * maleR[lower + 1]
    res <- (1 - maleR_temp) ^ t * maleMu_temp
    
  }
  
  res <- ifelse(res == 0,
                (1 - maleR[x + 1 + t]) ^ t * maleMu[x + 1 + t],
                res)
  
         
  return(res)
  
}



#### intensities ####


intensityAssumptionFct <- function(age, sex, ageRetirement, timeVector, muFT) {
  
  
  # disability
  muAI <- function(age, sex, ageRetirement) {

    (10 ^ (5.662015 + 0.033462 * age - 10)) * (age <= ageRetirement) 
  }
  
  # mortality - FT mortality model
  
  muAD <- approxfun(age + timeVector, muFT(age, timeVector, sex) ) # approxfun(age + 0:projectionYears, muFT(age, 0:projectionYears, sex) ) 

  muID <- function(x) {

    # if ( x <= ageRetirement ) {
    #   0.010339 + 10 ^ (5.070927 + 0.05049 * x - 10)
    # } else {
    #   approxfun(age + timeVector, muFT(age, timeVector, sex))(x)
    # }

    ifelse(x <= ageRetirement, 
           0.010339 + 10 ^ (5.070927 + 0.05049 * x - 10),
           muAD(x))
    
  }
  
  
  # reactivation
  muIA   <- function(age, ageRetirement) 4.0116 * exp(-0.117 * age) * (age <= ageRetirement)
  
  # free policy 
  muAF <- function(age, ageRetirement) 0.05 * (age <= ageRetirement)
  
  # surrender
  
  #surRates$rate[surRates$age >= ageRetirement] <- 0
  muAS <- function(age, ageRetirement) (0.06 - 0.002 * max(age - 40, 0)) * (age <= ageRetirement)
  
  
  # mu - intensity function of time and states
  mu12 <- function(age) muAI(age, sex, ageRetirement)
  mu13 <- function(age) muAD(age)
  mu21 <- function(age) muIA(age, ageRetirement)
  mu23 <- function(x) muID(x)
  mu14 <- function(age) muAF(age, ageRetirement)
  mu45 <- function(age) muAI(age, sex, ageRetirement)
  mu46 <- function(age) muAD(age)
  mu54 <- function(age) muIA(age, ageRetirement)
  mu56 <- function(age) muID(age)
  mu17 <- function(age) muAS(age, ageRetirement)
  mu48 <- function(age) muAS(age, ageRetirement)

  mu <- function(s, from, to) {
    ifelse(from == 1 & to == 2,
           mu12(age + s),
           ifelse(from == 1 & to == 3,
                  mu13(age + s),
                  ifelse(from == 2 & to == 1,
                         mu21(age + s),
                         ifelse(from == 2 & to == 3,
                                mu23(age + s),
                                ifelse(from == 1 & to == 4,
                                       mu14(age + s),
                                       ifelse(from == 4 & to == 5,
                                              mu45(age + s),
                                              ifelse(from == 4 & to == 6,
                                                     mu46(age + s),
                                                     ifelse(from == 5 & to == 4,
                                                            mu54(age + s),
                                                            ifelse(from == 5 & to == 6,
                                                                   mu56(age + s), 
                                                                   ifelse(from == 1 & to == 7,
                                                                          mu17(age + s),
                                                                          ifelse(from == 4 & to == 8,
                                                                                 mu48(age + s),
                                                                                 0)))))))))))
  }
  
  
  return(
    list(
      mu12 = mu12,
      mu13 = mu13,
      mu21 = mu21,
      mu23 = mu23,
      mu14 = mu14,
      mu45 = mu45,
      mu46 = mu46,
      mu54 = mu54,
      mu56 = mu56,
      mu17 = mu17,
      mu48 = mu48,
      mu   = mu
    )
  )
  
  
}


## technical basis

intensityTechAssumptionFct <- function(age, sex, ageRetirement, r) {

  # disability
  muAI <- function(age, ageRetirement, sex) {
   
     (0.0004 + 10^(4.54 + 0.06 * age - 10)) * (age <= ageRetirement)
    
  }
  
  # mortality
  muAD <- function(age, sex) {
    
    0.0005 + 10^(5.88 + 0.038 * age - 10)
    
  }
  
  muID <- function(age, ageRetirement, sex) {
    
    muAD(age, sex) * (1 + (age <= ageRetirement)) 
    
  }
  
  # reactivation
  muIA <- function(age, ageRetirement) {
    
    (2.0058 * exp(-0.117 * age)) * (age <= ageRetirement)
    
  }
  
  
  # mu - intensity function of age
  mu12 <- function(age) muAI(age, ageRetirement, sex)
  mu13 <- function(age) muAD(age, sex)
  mu21 <- function(age) muIA(age, ageRetirement)
  mu23 <- function(age) muID(age, ageRetirement, sex)
  
  mu <- function(s, from, to) {
    if (from == 1 && to == 2) {
      return( mu12(age + s) )
    } else if (from == 1 && to == 3) {
      return( mu13(age + s) )
    } else if (from == 2 && to == 1) {
      return( mu21(age + s) )
    } else if (from == 2 && to == 3) {
      return( mu23(age + s) )
    } else {
      return( 0 )
    }
  }
    
    return(
      list(
        mu12 = mu12,
        mu13 = mu13,
        mu21 = mu21,
        mu23 = mu23,
        mu = mu
      )
    )
    
    
}




## market basis no PHB

intensityNoPHBAssumptionFct <- function(age, sex, ageRetirement, timeVector, muFT) {
  
  # disability
  muAI <- function(age, sex, ageRetirement) {
    
    (10 ^ (5.662015 + 0.033462 * age - 10)) * (age <= ageRetirement) 
  }
  
  # mortality - FT mortality model
  
  muAD <- approxfun(age + timeVector, muFT(age, timeVector, sex) )
  
  muID <- function(x) {
    
    if ( x <= ageRetirement ) {
      0.010339 + 10 ^ (5.070927 + 0.05049 * x - 10)
    } else {
      approxfun(age + timeVector, muFT(age, timeVector, sex))(x)
    }
    
  }
  
  
  # reactivation
  muIA   <- function(age, ageRetirement) 4.0116 * exp(-0.117 * age) * (age <= ageRetirement)
  
  
  # mu - intensity function of time and states
  mu12 <- function(age) muAI(age, ageRetirement, sex)
  mu13 <- function(age) muAD(age)
  mu21 <- function(age) muIA(age, ageRetirement)
  mu23 <- function(x) muID(x)
  
  mu <- function(s, from, to) {
    if (from == 1 && to == 2) {
      return( mu12(age + s) )
    } else if (from == 1 && to == 3) {
      return( mu13(age + s) )
    } else if (from == 2 && to == 1) {
      return( mu21(age + s) )
    } else if (from == 2 && to == 3) {
      return( mu23(age + s) )
    } else {
      return( 0 )
    }
  }
  
  return(
    list(
      mu12 = mu12,
      mu13 = mu13,
      mu21 = mu21,
      mu23 = mu23,
      mu   = mu
    )
  )
  
  
}



## test basis

intensityTestAssumptionFct <- function(age, sex) {
  
  
  # mortality
  muAD <- function(age, sex) {
    
    0.0005 + 10^(5.88 + 0.038 * age - 10)
    
  }

  
  
  # mu - intensity function of time and states
  mu12 <- function(age) muAD(age, sex)

  mu <- function(s, from, to) {
    if (from == 1 && to == 2) {
      return( mu12(age + s) )
    } else {
      return( 0 )
    }
  }
  
  return(
    list(
      mu12 = mu12,
      mu   = mu
    )
  )
  
  
}






