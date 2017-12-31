
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;



namespace assumptions_market {

// disability

double muAI_C(double age, int sex, double age_retirement) {
  
  double res = 0;
  
  res = pow(10, 5.662015 + 0.033462 * age - 10) * (age <= age_retirement);
  
  return res;
  
}


// mortality - Danish fsa model


double muAD_C(double time, int sex, vec& muAD_vector, double dt){

  double index;
  double linear_combination;
  int lower;
  int upper;
  
  index = time / dt;
  lower = floor(index);
  upper = ceil(index);
  
  if ((index == lower) & (index == upper)) {
    linear_combination = muAD_vector[index];
  } else {
    linear_combination = (index - lower) * muAD_vector[upper] + (upper - index)* muAD_vector[lower];
  }
  
  return linear_combination;
  
}


double muID_C(double age, int sex, double time, double age_retirement, vec& muAD_vector, double dt){
  
  double res = 0;
  
  if (age <= age_retirement) {
    res = 0.010339 + pow(10, 5.070927 + 0.05049 * age - 10);
  } else {
    res = muAD_C(time, sex, muAD_vector, dt);
  }
  
  return res;
  
}


// reativation

double muIA_C(double age, int sex, double age_retirement) {
  
  double res = 0;
  
  res = 4.0116 * exp(-0.117 * age) * (age <= age_retirement);
  
  return res;
  
}



// free policy

double muAF_C(double age, double age_retirement) {
  
  double res = 0;
  
  res = 0.05 * (age <= age_retirement);
  
  return res;
  
}



// surrender

double muAS_C(double age, double age_retirement) {
  
  double res = 0;
  vec max_input = zeros(2);
  
  max_input(0) = age - 40;
  max_input(1) = 0;
  
  res = (0.06 - 0.002 * max(max_input)) * (age <= age_retirement);
  
  return res;
}


}


