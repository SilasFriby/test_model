
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;



namespace assumptions_technical {


// disability

double muAI_C(double age, int sex, double age_retirement) {
  
  double res = 0;
  
  res = (0.0004 + pow(10, 4.54 + 0.06 * age - 10)) * (age <= age_retirement);
  
  return res;
  
}


// mortality

double muAD_C(double age, int sex) {
  
  double res = 0;
  
  res = 0.0005 + pow(10, 5.88 + 0.038 * age - 10);
  
  return res;
  
}

double muID_C(double age, int sex, double age_retirement) {
  
  double res = 0;
  
  res = (0.0005 + pow(10, 5.88 + 0.038 * age - 10)) * (1 + (age <= age_retirement));
  
  return res;
  
}


// reativation

double muIA_C(double age, int sex, double age_retirement) {
  
  double res = 0;
  
  res = (2.0058 * exp(-0.117 * age)) * (age <= age_retirement);
  
  return res;
  
}



double mu_cpp(double s, int from, int to, double age, int sex, double age_retirement){
  
  
  double res = 0;
  
  
  if (from == 1 && to == 2) {
    res = muAI_C(age + s, sex, age_retirement);
  } else if (from == 1 && to == 3) {
    res = muAD_C(age + s, sex);
  } else if (from == 2 && to == 1) {
    res = muIA_C(age + s, sex, age_retirement);
  } else if (from == 2 && to == 3) {
    res = muID_C(age + s, sex, age_retirement);
  } else {
    res = 0;
  }
  
  
  return res;
  
}



}