

#### FTs levetidsmodel

## data

# benchmarkFT <- read.table("data/Benchmark_doedelighed2015.csv", header = T, sep = ';', dec = ',')
# benchmarkR  <- read.table("data/Benchmark-levetidsforbedringer2015.csv", header = T, sep = ';', dec = ',')
# head(benchmarkFT)
# head(benchmarkR)
# 
# initial_mortality <- benchmarkFT[, 2]
# R <- benchmarkR[, 2]
# age <- benchmarkFT[, 1]
# betaVector <- c(0, 0, 0)

#' Danish FSA mortality model
#'
#' @param t - years into the future (double)
#' @param age - vector of ages (numeric vector)
#' @param initial_mortality - vector of mortalities (numeric vector)
#' @param R - vector of death decline factors (numeric vector)
#' @param beta_vector - vector of beta's from statistical procedure (numeric vector)
#'
#' @return updated mortality benchmark
#'
#' @examples
ft_benchmark_fct <- function(t, age, initial_mortality, R, beta_vector) {
  
  
  
  ## regressors
  
  regressor_fct <- function(i, x) {
   
    xi <- 20 * (2 + 0:3) 
    
    ifelse( x <= xi[i - 1 + 1], 
            1, 
            ifelse( x > xi[i - 1 + 1] & x < xi[i + 1], 
                    (xi[i + 1] - x) / (xi[i + 1] - xi[i - 1 + 1]), 
                    0) 
            )
    
  }
  
  
  
  ## central mortality
  
  if ( all(beta_vector == 0) ) {
    central_mortality <- head(initial_mortality, -1)
  } else {
    central_mortality <- ( tail(initial_mortality, -1) + head(initial_mortality, -1) ) / 2
  }
  
  
  ## population adjustments - using beta's obtained from statistical analysis
  
  model_adjustment <- exp( beta_vector[1] * regressor_fct(1, age) +
                              beta_vector[2] * regressor_fct(2, age) + 
                              beta_vector[3] * regressor_fct(3, age) ) 
    
  ## model mortality
  
  if ( all(beta_vector == 0) ) {
    
    model_mortality <- head(initial_mortality, -1)
    
  } else {
    
    model_mortality <-  ( tail(model_adjustment, -1) * tail(central_mortality, -1) + 
                            head(model_adjustment, -1) * head(central_mortality, -1) ) / 2
    
    model_mortality <- c(initial_mortality[1] * model_adjustment[1], model_mortality)
    
  }
  
  
                            
  
  
  ## project t years forward
  
  final_mortality <- model_mortality * (1 - R) ^ t
  
  return(model_mortality)

  
}


