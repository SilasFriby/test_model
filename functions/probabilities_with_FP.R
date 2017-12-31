



probabilityProjectionFct_FP <- function(timeProbability, probInitial, stateCombi, nStates,
                                        mu, w1Vector, states_FP, ... ) {
  
  
  #### step size (if negative we go backwards in time)
  
  nProbability <- length(timeProbability)
  h            <- (timeProbability[nProbability] -  timeProbability[1])/(nProbability - 1)
  
  #### prepare for runge-kutta algorithm
  
  ## (1) p_i: prob of being in state i at time s given state 1 at time zero - solve kolmogorov
  
  # result matrix
  p_i_initial     <- probInitial
  p_i_matrix      <- matrix(NA, nrow = nProbability, ncol = length(p_i_initial)) # results saved in matrix [time, dimension]
  p_i_matrix[1, ] <- p_i_initial
  
  p_i_initial_FP <- rep(0, nStatesFP)
  p_i_matrix_FP <- matrix(NA, nrow = nProbability, ncol = length(p_i_initial_FP))
  p_i_matrix_FP[1, ] <- p_i_initial_FP
  
  # differential function
  f_p_i <- function(s, x) klm(s, x, mu)
  rhoFct =  approxfun(timeProbability, w1Vector)
  
  #### runge-kutta algorithm
  
  #pb <- txtProgressBar(min = 1, max = nProbability, style = 3)
  for (step in 1:(nProbability - 1)) {
    #foreach(step = 1:(nProbability - 1), .combine = rbind) %dopar% {
    
    s <- timeProbability[step] # i'th time point
    
    ## (1) p_i
    
    foo_p_i <- p_i_matrix[step, ]
    
    F1_p_i <- f_p_i(s, foo_p_i)
    F2_p_i <- f_p_i(s + 1 / 2 * h, foo_p_i + 1 / 2 * h * F1_p_i)
    F3_p_i <- f_p_i(s + 1 / 2 * h, foo_p_i + 1 / 2 * h * F2_p_i)
    F4_p_i <- f_p_i(s + h, foo_p_i + h * F3_p_i)
    
    p_i_matrix[step + 1, ] <- foo_p_i + h * 1 / 6 * (F1_p_i + 2 * F2_p_i + 2 * F3_p_i + F4_p_i)
    
    
    ## (2) p_i_FP
    
    probActiveFct <- approxfun(timeProbability[1:(step + 1)], p_i_matrix[1:(step + 1), 1])
    f_p_i_FP  <- function(s, x) klm_FP(s, x, mu, states_FP, rhoFct, probActiveFct)
    
    foo_p_i_FP <- p_i_matrix_FP[step, ]
    
    F1_p_i_FP <- f_p_i_FP(s, foo_p_i_FP)
    F2_p_i_FP <- f_p_i_FP(s + 1 / 2 * h, foo_p_i_FP + 1 / 2 * h * F1_p_i_FP)
    F3_p_i_FP <- f_p_i_FP(s + 1 / 2 * h, foo_p_i_FP + 1 / 2 * h * F2_p_i_FP)
    F4_p_i_FP <- f_p_i_FP(s + h - 10^(-10), foo_p_i_FP + h * F3_p_i_FP)
    
    p_i_matrix_FP[step + 1, ] <- foo_p_i_FP + h * 1 / 6 * (F1_p_i_FP + 2 * F2_p_i_FP + 2 * F3_p_i_FP + F4_p_i_FP)
    
    
    #setTxtProgressBar(pb,step)
  }
  #close(pb)
  
  
  
  return(
    list(
      "pMatrix" = p_i_matrix,
      "pMatrix_FP" = p_i_matrix_FP
    )
  )
  
}
