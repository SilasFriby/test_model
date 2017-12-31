#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

  

// [[Rcpp::export]]
vec f1(vec& y) {
  
  int n = y.size();
  vec res = zeros(n);
  for (int i = 0; i < n; i++) {
    res[i] = pow(y[i], 2);
  }
  
  return res;
  
}


// [[Rcpp::export]]
vec f2(vec& y) {
  
  uword n = y.size();
  vec res = zeros(n);
  for (uword i = 0; i < n; i++) {
    res[i] = pow(y[i], 2);
  }
  
  return res;
  
}


// // [[Rcpp::export]]
// vec f2(vec& y) {
//   
//   int n = y.size();
//   vec res = zeros(n);
//   
//   for(arma::vec::iterator it = y.begin();
//       it != y.end(); ++it)
//   {
//     int index = it - y.begin();
//     res[index] = pow(*it, 2);
//   }
//   
//   return res;
//   
// }



