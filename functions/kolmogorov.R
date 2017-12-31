

## Kolmogorov forward differential equation

klm <- function(s, x, mu){
  
  states <- length(x)   
  xDeriv <- rep(0, states)
  for (j in 1:states) {
    for (k in 1:states) {
      if (k != j) {
        xDeriv[j]    <- xDeriv[j] - x[j] * mu(s, j, k) + x[k] * mu(s, k, j)
      }
    }
  }
  
  return(xDeriv)
}


klm_FP <- function(s, x, mu, states_FP, rhoFct, probActiveFct){
  
  xDeriv <- rep(0, length(states_FP))
  for (j in states_FP) {
    index_j <- which(states_FP == j)
    xDeriv[index_j] <- xDeriv[index_j] + (j == states_FP[1]) * probActiveFct(s) * mu(s, 1, states_FP[1]) * rhoFct(s)
    for (k in states_FP) {
      if (k != j) {
        index_k <- which(states_FP == k)
        xDeriv[index_j] <- xDeriv[index_j] - x[index_j] * mu(s, j, k) + x[index_k] * mu(s, k, j)
      }
    }
    if (j == states_FP[1]) { # include the intensity from 4 to 8 - not dynamically written
      xDeriv[index_j] <- xDeriv[index_j] - x[index_j] * mu(s, j, 8)
    }
  }
  
  return(xDeriv)
}







## Aggregated probability of jumping from i to j (i !=j )  

P_ij_fct <- function(s, x, mu, pMatrix, timeVec, nStates, stateCombi){
  # pFct: fct of time and state - gives the prob of being in some state at some time
  
  states <- nStates 
  xDeriv <- matrix(NA, states, states)
  
  for (i in 1:states) {
    for (j in 1:states) {
      if (j != i & paste(i, j, sep = '') %in% stateCombi ) {
        xDeriv[i, j] <- approxfun(timeVec, pMatrix[,i], rule = 2)(s) * mu(s, i, j)
      }
    }
  }
  
  xDeriv <- as.vector(t(xDeriv))
  xDeriv <- xDeriv[!is.na(xDeriv)]
  
  return(xDeriv)
  
}

## Aggregated probability of leaving state i

P_i._fct <- function(s, x, mu, p, nStates){
  # pFct: fct of time and state - gives the prob of being in some state at some time
  
  states <- nStates
  xDeriv <- rep(NA, states)
  
  for (i in 1:states) {
        mu_i <- sum( sapply(setdiff(1:states, i), function(k) mu(s, i, k)) ) # sum of intensities away from i
        xDeriv[i] <- splinefun(p$time, p$xHist[,i])(s) * mu_i 
  }
  
  return(xDeriv)
  
}




