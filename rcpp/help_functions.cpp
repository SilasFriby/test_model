#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

namespace help_functions {


vec arma_reverse(vec x) {
  
  int n = x.size();
  vec res = zeros(n);
  
  for (int i = 0; i < n; i++) {
    int index_start = i;
    int index_end = (n - 1) - i;
    res[index_start] = x[index_end];
  }
  
  return res;
  
}

arma::colvec arma_sub_cond(arma::colvec x, arma::colvec y, double val, double replace) {
  
  arma::uvec ids = find(y >= val); // Find indices
  
  x.elem(ids).fill(replace);       // Assign value to condition
  
  return x;
  
}

arma::colvec arma_sub_cond_less(arma::colvec x, arma::colvec y, double val, double replace) {
  
  arma::uvec ids = find(y < val); // Find indices
  
  x.elem(ids).fill(replace);       // Assign value to condition
  
  return x;
  
}


}