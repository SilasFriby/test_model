

## Runge-Kutta solver

rksolve <- function(x0,f,timeVector) {
  # x0 - starting values
  # f - differential equation function: function(time, x)
  
  n <- length(timeVector)
  a <- timeVector[1]
  b <- timeVector[n]
  
  h <- (b-a)/(n-1) # step size (if negative we go backwards in time)
  
  xVec     <- matrix(NA, nrow=n, ncol=length(x0)) # results saved in matrix [time, dimension]
  xVec[1,] <- x0
  #pb <- txtProgressBar(min = 1, max = n, style = 3)
  for (i in 1:(n-1)){
    
    xFoo <- xVec[i, ]
    s    <- a + h * (i-1) # i'th time point
    
    F1 <- f(s, xFoo)
    F2 <- f(s + 1/2*h, xFoo + 1/2*h*F1)
    F3 <- f(s + 1/2*h, xFoo + 1/2*h*F2)
    F4 <- f(s + h, xFoo + h*F3)
    
    xVec[i+1, ] <- xFoo + h * 1/6 * (F1 + 2*F2 + 2*F3 + F4)
    
    #setTxtProgressBar(pb,i)
  }
  #close(pb)
  return(list("time"=a+h*(0:(n-1)), "xHist"=xVec, "x" = xVec[n,]))
}


## Euler

euler <- function(x0,f,a,b,n) {
  # x0 - starting values
  # f - differential equation function: function(time, x)
  # a - start time
  # b - end time
  # n - number of time points 
  
  h <- (b-a)/(n-1) # step size (if negative we go backwards in time)
  
  xVec     <- matrix(NA, nrow=n, ncol=length(x0)) # results saved in matrix [time, dimension]
  xVec[1,] <- x0
  #pb <- txtProgressBar(min = 1, max = n, style = 3)
  for (i in 1:(n-1)){
    
    xFoo <- xVec[i, ]
    s    <- a + h * (i-1) # i'th time point
    xVec[i+1, ] <- xFoo + h * f(s, xFoo)
    
    #setTxtProgressBar(pb,i)
  }
  #close(pb)
  return(list("time"=a+h*(0:(n-1)), "xHist"=xVec, "x" = xVec[n,]))
}






