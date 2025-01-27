#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::interfaces(r, cpp)]]

//' Evalute matrics
//' @name eva
//' @param num_Sim A numeric vector of simulated values.
//' @param num_Obs A numeric vector of observed values. NA values are removed along with corresponding values in num_Sim.
//' @return A double representing the Evalute matrics.
//' @export
// [[Rcpp::export]]
double eva_NSE(NumericVector num_Sim, NumericVector num_Obs) {
  // Ensure both vectors have the same length
  if (num_Sim.size() != num_Obs.size()) {
    stop("Simulated and observed vectors must have the same length.");
  }

  // Filter out NA values
  NumericVector filtered_Sim, filtered_Obs;
  for (int i = 0; i < num_Obs.size(); i++) {
    if (!NumericVector::is_na(num_Obs[i])) {
      filtered_Obs.push_back(num_Obs[i]);
      filtered_Sim.push_back(num_Sim[i]);
    }
  }

  // Ensure there are remaining values to compute NSE
  if (filtered_Obs.size() == 0) {
    stop("All observed values are NA; cannot calculate NSE.");
  }

  // Calculate the mean of observed values
  double obs_mean = mean(filtered_Obs);

  // Compute numerator and denominator for NSE
  double numerator = 0.0, denominator = 0.0;
  for (int i = 0; i < filtered_Sim.size(); i++) {
    numerator += pow(filtered_Obs[i] - filtered_Sim[i], 2);
    denominator += pow(filtered_Obs[i] - obs_mean, 2);
  }

  // Return NSE
  return 1.0 - (numerator / denominator);
}

//' @rdname eva
//' @param factor_r,factor_alpha,factor_beta A double specifying the weight for the correlation term (r - 1), (alpha - 1) and (beta - 1). Default is 1.0.
//' @export
// [[Rcpp::export]]
double eva_KGE(NumericVector num_Sim, NumericVector num_Obs,
               double factor_r = 1.0, double factor_alpha = 1.0, double factor_beta = 1.0) {
  // Ensure both vectors have the same length
  if (num_Sim.size() != num_Obs.size()) {
    stop("Simulated and observed vectors must have the same length.");
  }

  // Filter out NA values
  NumericVector filtered_Sim, filtered_Obs;
  for (int i = 0; i < num_Obs.size(); i++) {
    if (!NumericVector::is_na(num_Obs[i])) {
      filtered_Obs.push_back(num_Obs[i]);
      filtered_Sim.push_back(num_Sim[i]);
    }
  }

  // Ensure there are remaining values to compute KGE
  if (filtered_Obs.size() == 0) {
    stop("All observed values are NA; cannot calculate KGE.");
  }

  // Calculate means of observed and simulated values
  double obs_mean = mean(filtered_Obs);
  double sim_mean = mean(filtered_Sim);

  // Calculate standard deviations
  double obs_sd = sd(filtered_Obs);
  double sim_sd = sd(filtered_Sim);

  // Calculate correlation coefficient
  double numerator = 0.0;
  for (int i = 0; i < filtered_Sim.size(); i++) {
    numerator += (filtered_Obs[i] - obs_mean) * (filtered_Sim[i] - sim_mean);
  }
  double correlation = numerator / ((filtered_Sim.size() - 1) * obs_sd * sim_sd);

  // Calculate KGE components
  double alpha = sim_sd / obs_sd;
  double beta = sim_mean / obs_mean;

  // Apply factors to each component
  double sum_Factor = factor_r + factor_alpha + factor_beta;
  factor_r = factor_r / sum_Factor;
  factor_alpha = factor_alpha / sum_Factor;
  factor_beta = factor_beta / sum_Factor;


  double term_r = pow(factor_r * (correlation - 1.0), 2);
  double term_alpha = pow(factor_alpha * (alpha - 1.0), 2);
  double term_beta = pow(factor_beta * (beta - 1.0), 2);

  // Calculate and return KGE
  return 1.0 - sqrt(term_r + term_alpha + term_beta);
}
