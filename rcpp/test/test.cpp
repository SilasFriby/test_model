#include <RcppArmadillo.h>
#include <assumptions_market.cpp>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace assumptions_market;



double mu_cpp(double s, int from, int to, double age, int sex, double age_retirement, vec& muAD_vector, double dt){
  
  
  double res = 0;
  
  
  if (from == 1 && to == 2) {
    res = assumptions_market::muAI_C(age + s, sex, age_retirement);
  } else if (from == 1 && to == 3) {
    res = assumptions_market::muAD_C(s, sex, muAD_vector, dt);
  } else if (from == 2 && to == 1) {
    res = assumptions_market::muIA_C(age + s, sex, age_retirement);
  } else if (from == 2 && to == 3) {
    res = assumptions_market::muID_C(age + s, sex, s, age_retirement, muAD_vector, dt);
  } else if (from == 1 && to == 4) {
    res = assumptions_market::muAF_C(age + s, age_retirement);
  } else if (from == 4 && to == 5) {
    res = assumptions_market::muAI_C(age + s, sex, age_retirement);
  } else if (from == 4 && to == 6) {
    res = assumptions_market::muAD_C(s, sex, muAD_vector, dt);
  } else if (from == 5 && to == 4) {
    res = assumptions_market::muIA_C(age + s, sex, age_retirement);
  } else if (from == 5 && to == 6) {
    res = assumptions_market::muID_C(age + s, sex, s, age_retirement, muAD_vector, dt);
  } else if (from == 1 && to == 7) {
    res = assumptions_market::muAS_C(age + s, age_retirement);
  } else if (from == 4 && to == 8) {
    res = assumptions_market::muAS_C(age + s, age_retirement);
  } else {
    res = 0;
  }
  
  
  return res;
  
}






vec kolmogorov_cpp(double s, const vec x, double age, int sex, double age_retirement, vec& muAD_vector, 
                   double dt){
  
  int n = x.size();
  vec xDeriv = zeros(n);
  
  for (int j = 0; j < n; ++j) {
    for (int k = 0; k < n; ++k) {
      if ( k != j ){
        xDeriv[j] = xDeriv[j] - x[j] * mu_cpp(s, j + 1, k + 1, age, sex, age_retirement, muAD_vector, dt) +
          x[k] * mu_cpp(s, k + 1, j + 1, age, sex, age_retirement, muAD_vector, dt);
      }
      
    }
  }
  
  return xDeriv;
}



vec kolmogorov_FP_cpp(double s, vec x, double age, int sex, double age_retirement, vec& muAD_vector,
                      double dt, vec states_FP, double rho, double p_active){
  
  int n = x.size();
  vec xDeriv = zeros(n);
  
  for (int j = 0; j < n; ++j) {
    xDeriv[j] = xDeriv[j] + ((j + 4) == states_FP[0]) * p_active * 
      mu_cpp(s, 1, states_FP[0], age, sex, age_retirement, muAD_vector, dt) * rho;
    for (int k = 0; k < n; ++k) {
      if (k != j) {
        xDeriv[j] = xDeriv[j] - x[j] * mu_cpp(s, j + 4, k + 4, age, sex, age_retirement, muAD_vector, dt)
        + x[k] * mu_cpp(s, k + 4, j + 4, age, sex, age_retirement, muAD_vector, dt);
      }
    }
    if (j + 4 == states_FP[0]) { // include the intensity from 4 to 8 - not dynamically written
      xDeriv[j] = xDeriv[j] - x[j] * mu_cpp(s, j + 4, 8, age, sex, age_retirement, muAD_vector, dt);
    }
  }
  
  return xDeriv;
  
}


// [[Rcpp::export]]
mat probabilities_FP_cpp(vec& timeProbability, vec& probInitial, CharacterVector& stateCombi,
                         int& nStates, double& age, int& sex, double& age_retirement, vec& muAD_vector, 
                         double& dt, vec& states_FP, vec& rho) {
  
  ////// step size (if negative we go backwards in time)
  
  int nProbability = timeProbability.size();
  double h = (timeProbability[nProbability - 1] -  timeProbability[0])/(nProbability - 1);
  
  ////// prepare for runge-kutta algorithm
  
  //// p_i: prob of being in state i at time s given state 1 at time zero - solve kolmogorov
  
  // result matrix
  vec p_i_initial = probInitial;
  mat p_i_matrix = zeros(nProbability, nStates); // results saved in matrix [time, dimension]
  p_i_matrix.row(0) = p_i_initial.t();
  
  int nStates_FP = states_FP.size();  
  mat p_i_matrix_FP = zeros(nProbability, nStates_FP); // results saved in matrix [time, dimension]
  
  
  ////// runge-kutta algorithm
  
  
  vec foo, foo_p_i, F1_p_i, F2_p_i, F3_p_i, F4_p_i;
  vec foo_FP, foo_p_i_FP, F1_p_i_FP, F2_p_i_FP, F3_p_i_FP, F4_p_i_FP;
  double s;
  
  for (int step = 0; step < (nProbability - 1); ++step) {
    
    s = timeProbability[step];
    
    //// (1) p_i
    
    foo_p_i = p_i_matrix.row(step).t();
    F1_p_i = kolmogorov_cpp(s, foo_p_i, age, sex, age_retirement, muAD_vector, dt);
    F2_p_i = kolmogorov_cpp(s + 1.0 / 2.0 * h, foo_p_i + 1.0 / 2.0 * h * F1_p_i, age, sex, age_retirement, muAD_vector, dt);
    F3_p_i = kolmogorov_cpp(s + 1.0 / 2.0 * h, foo_p_i + 1.0 / 2.0 * h * F2_p_i, age, sex, age_retirement, muAD_vector, dt);
    F4_p_i = kolmogorov_cpp(s + h, foo_p_i + h * F3_p_i, age, sex, age_retirement, muAD_vector, dt);
    
    foo = foo_p_i + h * 1.0 / 6.0 * (F1_p_i + 2.0 * F2_p_i + 2.0 * F3_p_i + F4_p_i);
    
    p_i_matrix.row(step + 1) = foo.t();
    
    
    // (2) p_i_FP
    
    double index1 = round(s / dt);
    
    foo_p_i_FP = p_i_matrix_FP.row(step).t();
    
    double foo_rho1 = rho[index1];
    double p_active1 = p_i_matrix(step, 0);
    F1_p_i_FP = kolmogorov_FP_cpp(s, foo_p_i_FP, age, sex, age_retirement, muAD_vector, dt, states_FP, foo_rho1, p_active1);
    
    double s2 = s + 1.0 / 2.0 * h;
    double index2 = s2 / dt;
    int lower2 = floor(index2);
    int upper2 = ceil(index2);
    double foo_rho2 = (index2 - lower2) * rho[upper2] + (upper2 - index2)* rho[lower2];
    double p_active2 = (index2 - lower2) * p_i_matrix(step + 1, 0) + (upper2 - index2) * p_i_matrix(step, 0);
    F2_p_i_FP = kolmogorov_FP_cpp(s2, foo_p_i_FP + 1.0 / 2.0 * h * F1_p_i_FP, age, sex, age_retirement, muAD_vector, dt, states_FP, foo_rho2, p_active2);
    
    double s3 = s + 1.0 / 2.0 * h;
    double index3 = s3 / dt;
    int lower3 = floor(index3);
    int upper3 = ceil(index3);
    double foo_rho3 = (index3 - lower3) * rho[upper3] + (upper3 - index3)* rho[lower3];
    double p_active3 = (index3 - lower3) * p_i_matrix(step + 1, 0) + (upper3 - index3) * p_i_matrix(step, 0);
    F3_p_i_FP = kolmogorov_FP_cpp(s3, foo_p_i_FP + 1.0 / 2.0 * h * F2_p_i_FP, age, sex, age_retirement, muAD_vector, dt, states_FP, foo_rho3, p_active3);
    
    double s4 = s + h;
    double index4 = round(s4 / dt);
    double foo_rho4 = rho[index4];
    double p_active4 = p_i_matrix(step + 1, 0);
    F4_p_i_FP = kolmogorov_FP_cpp(s4, foo_p_i_FP + h * F3_p_i_FP, age, sex, age_retirement, muAD_vector, dt, states_FP, foo_rho4, p_active4);
    
    foo_FP = foo_p_i_FP + h * 1.0 / 6.0 * (F1_p_i_FP + 2.0 * F2_p_i_FP + 2.0 * F3_p_i_FP + F4_p_i_FP);
    
    p_i_matrix_FP.row(step + 1) = foo_FP.t();
    
    
  }
  
  int m = states_FP.size(); 
  int index1;
  int index2;
  for (int i = states_FP[0]; i <= states_FP[m - 1]; i++) {
    index1 = i - 1;
    index2 = i - 1 - m;
    p_i_matrix.col(index1) = p_i_matrix_FP.col(index2);
  } 
  
  return p_i_matrix;
  
  // return List::create(
  //   _["p_matrix"] = p_i_matrix,
  //   _["p_matrix_FP"] = p_i_matrix_FP
  // );
  
}

