
#include <RcppArmadillo.h>
#include <assumptions_technical.cpp>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace assumptions_technical;




// Thiele differential equation

namespace thiele_technical {

vec thiele_diff_cpp(double s, vec x, mat& premium, mat& benefit_continuous, mat& benefit_discrete, double r, double dt,
                    double age, int sex, double age_retirement){

  int n_states = x.size();
  vec res = zeros(n_states);
  int index;
  index = floor(s / dt);

  // interest rate and continuous payments 
  for (int i = 0; i < n_states; i++) {
    res[i] = res[i] + r * x[i] - benefit_continuous(index, i) + premium(index, i);
  }
    
  // discrete payments
  for (int i = 0; i < n_states; i++) {
    for (int j = 0; j < n_states; j++) {
      mat benefit_discrete_temp = benefit_discrete.rows(index * n_states, index * n_states + n_states - 1);
      res[i] = res[i] - assumptions_technical::mu_cpp(s, i + 1, j + 1, age, sex, age_retirement) *
        (x[j] - x[i] + benefit_discrete_temp(i, j));

    }
  }
    
    return res;
}



mat thiele_cpp(vec& time, vec& x, mat& premium, mat& benefit_continuous, mat& benefit_discrete,
               double& r, double& dt, double& age, int& sex, double& age_retirement) {

  // result matrix
  int n = time.size();
  double h = (time[n - 1] -  time[0])/(n - 1);
  mat res = zeros(n, x.size()); // results saved in matrix [time, dimension]
  res.row(0) = x.t();

  vec foo, foo_res, F1, F2, F3, F4;
  double s;
  
  for (int step = 0; step < (n - 1); ++step) {
    
    s = time[step];
    
    //// (1) p_i
    
    foo = res.row(step).t();
    F1 = thiele_diff_cpp(s, foo, premium, benefit_continuous, benefit_discrete, r, dt, age, sex, age_retirement);
    F2 = thiele_diff_cpp(s + 1.0 / 2.0 * h, foo + 1.0 / 2.0 * h * F1, premium, benefit_continuous, benefit_discrete, r, dt, age, sex, age_retirement);
    F3 = thiele_diff_cpp(s + 1.0 / 2.0 * h, foo + 1.0 / 2.0 * h * F2, premium, benefit_continuous, benefit_discrete, r, dt, age, sex, age_retirement);
    F4 = thiele_diff_cpp(s + h, foo + h * F3, premium, benefit_continuous, benefit_discrete, r, dt, age, sex, age_retirement);
    
    foo_res = foo + h * 1.0 / 6.0 * (F1 + 2 * F2 + 2 * F3 + F4);
    
    res.row(step + 1) = foo_res.t();
  
  }
  
  return res;

}

}



