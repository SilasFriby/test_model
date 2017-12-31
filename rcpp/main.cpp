
#include <RcppArmadillo.h>
#include <thiele_technical.cpp>
#include <probabilities_8_state_market.cpp>
#include <help_functions.cpp>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace thiele_technical;
using namespace probabilities_market;
using namespace help_functions;




// [[Rcpp::export]]
List main_cpp(colvec& time, vec& age_vector, vec& sex_vector, vec& age_retirement_vector, StringVector& state_vector,
              int& nStates, vec& premium_vector, vec& benefitDisability_vector, vec& benefitRetirement_vector,
              double& rate_technical, int& nStatesTech, double& dt, vec& muAD_vector, CharacterVector& stateCombi,
              vec& states_FP, bool display_progress = true){
  
  int n = time.size();
  int m = age_vector.size();
  
  mat benefit_continuous_matrix = zeros(n, nStatesTech);
  vec reserve_technical;
  vec reserve_technical_plus;
  vec rho;
  mat p_matrix;
  mat premium_matrix = zeros(n, nStatesTech);
  mat cf_matrix = zeros(n, m);

  
  Progress p(m, display_progress);
  for (int i = 0; i < m; i++) {
    
    if (Progress::check_abort() )
      return -1.0;
    
    double age = age_vector[i];
    int sex = sex_vector[i];
    double age_retirement = age_retirement_vector[i];
    vec probInitial = zeros(nStates);
    if (state_vector[i] == "active") probInitial(0) = 1;
    
    
    // annual benefits
    
    double premium = premium_vector[i];
    double benefitDisability = benefitDisability_vector[i];
    double benefitRetirement = benefitRetirement_vector[i];
    
    
    
    
    //// 1. order reserve at different times given active
    
    // mat premium_matrix = zeros(n, nStatesTech);
    // mat benefit_continuous_matrix = zeros(n, nStatesTech);
    mat benefit_discrete_matrix = zeros(nStatesTech * n, nStatesTech); // not implemented yet
    
    // premiums
    premium_matrix.col(0) = help_functions::arma_sub_cond_less(premium_matrix.col(0), 
                                  time, age_retirement - age, 
                                  premium);
    
    // disability benefits
    benefit_continuous_matrix.col(1) = help_functions::arma_sub_cond_less(benefit_continuous_matrix.col(1), 
                                  time, age_retirement - age, 
                                  benefitDisability);
    
    // retirement benefits
    benefit_continuous_matrix.col(0) = help_functions::arma_sub_cond(benefit_continuous_matrix.col(0), 
                                  time, age_retirement - age, 
                                  benefitRetirement);
    benefit_continuous_matrix.col(1) = help_functions::arma_sub_cond(benefit_continuous_matrix.col(1), 
                                  time, age_retirement - age, 
                                  benefitRetirement);
    
    
    // 1. order reserve 
    
    vec time_thiele = help_functions::arma_reverse(time);
    vec boundary_thiele = zeros(nStatesTech);
    reserve_technical = thiele_technical::thiele_cpp(time_thiele,
                                                        boundary_thiele,
                                                        premium_matrix,
                                                        benefit_continuous_matrix,
                                                        benefit_discrete_matrix,
                                                        rate_technical,
                                                        dt,
                                                        age,
                                                        sex,
                                                        age_retirement).col(0);
    reserve_technical = help_functions::arma_reverse(reserve_technical);

    // first order reserve given zero premium
    
    mat premium_matrix_plus = zeros(n, nStatesTech);
    reserve_technical_plus = thiele_technical::thiele_cpp(time_thiele,
                                                              boundary_thiele,
                                                              premium_matrix_plus,
                                                              benefit_continuous_matrix,
                                                              benefit_discrete_matrix,
                                                              rate_technical,
                                                              dt,
                                                              age,
                                                              sex,
                                                              age_retirement).col(0);
    reserve_technical_plus = help_functions::arma_reverse(reserve_technical_plus);

    // free policy factor
    
    rho = reserve_technical / reserve_technical_plus;


    //// state probabilities

    p_matrix = probabilities_market::probabilities_FP_cpp(time,
                                                        probInitial,
                                                        stateCombi,
                                                        nStates,
                                                        age,
                                                        sex,
                                                        age_retirement,
                                                        muAD_vector,
                                                        dt,
                                                        states_FP,
                                                        rho
    );
    

    
    //// cash flow

    // annual cash flow
    vec cf_annual =  zeros(n);
    
    double b1, b2, b3;
    for (int k = 0; k < n; k++) {

      if ( (time[k] + age) > age_retirement) {
        b1 = benefitRetirement;
        b2 = benefitRetirement;
        b3 = benefitRetirement;
      } else {
        b1 = -premium;
        b2 = benefitDisability;
        b3 = 0;
      }

      cf_annual[k] =
        p_matrix(k, 0) * (
            b1 + probabilities_market::mu_cpp(time[k], 1, 7, age, sex, age_retirement, muAD_vector, dt) *
              reserve_technical[k]
        ) +
          p_matrix(k, 1) * (
            b2
          ) +
            p_matrix(k, 3) * (
                b3 + probabilities_market::mu_cpp(time[k], 4, 8, age, sex, age_retirement, muAD_vector, dt) *
                  reserve_technical_plus[k]
            ) +
              p_matrix(k, 4) * (
                b2
              );


    }
    
    // save dt cash flow

    cf_matrix.col(i) = cf_annual * dt;
    
    p.increment();
      
  }
  
  return List::create(
    _["p_matrix"] = p_matrix,
    _["cf_matrix"] = cf_matrix,
    _["rho"] = rho,
    _["reserve_technical_plus"] = reserve_technical_plus,
    _["reserve_technical"] = reserve_technical,
    _["premium_matrix"] = premium_matrix
  );
  
  
}