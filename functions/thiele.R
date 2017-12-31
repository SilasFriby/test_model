
## Thiele differential equation

thiele <- function(s, x, premium, bCont, bDis, mu, r){
  # premium: the premium paid per year while in state 1
  # bCont: vector with entry i indicating the 'continuous' benefits paid in state i
  # bDis: matrix indicating the benefits paid once when jumping from i to j
  
  states <- length(x)
  
  # interest rate and cont payments
  xderiv <- r * x - benefitContFct(s, bCont) + 
    premiumFct(s, premium)
  
  # Discrete payments
  for (i in 1:states) {
    for (j in 1:states) {
      xderiv[i] <- xderiv[i] - mu(s, i, j) *
        (x[j] - x[i] + benefitDisFct(s, bDis)[i, j])
    }
  }
  
  return(xderiv)
}


# thiele_2state <- function(s, x, premium, bCont, mu){
#   
#   xderiv <- r * x - bCont * (s >= (retirementAge-40)) + 
#     premium * (s < (retirementAge-40)) + mu(s, 1, 2) * x
#   
#   return(xderiv)
# }
