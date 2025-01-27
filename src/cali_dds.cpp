#include <Rcpp.h>
#include <random>
#include <algorithm>
#include <cmath>
#include <vector>
// [[Rcpp::interfaces(r, cpp)]]

using namespace Rcpp;
// Function to generate a random number from a normal distribution
double generate_random(double mean, double stddev) {
  // Create a random number generator
  static std::random_device rd;  // Seed
  static std::mt19937 gen(rd()); // Mersenne Twister engine
  std::normal_distribution<> d(mean, stddev); // Normal distribution

  return d(gen); // Generate a number
}

// Usage in your DDS function
NumericVector generate_random_vector(int n, double mean, double stddev) {
  NumericVector result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = generate_random(mean, stddev);
  }
  return result;
}
// Helper function to generate a binary vector with probability P_i
std::vector<bool> generate_binary_vector(int n_x, double P_i, std::default_random_engine &generator) {
  std::bernoulli_distribution distribution(P_i);
  std::vector<bool> binary_vector(n_x);
  for (int i = 0; i < n_x; i++) {
    binary_vector[i] = distribution(generator);
  }
  return binary_vector;
}


//' This function implements a calibration algorithm based on the Direct Search method.
//'
//' It attempts to find the optimal parameter values that minimize the given objective (fitness) function.
//'
//' @param fitness A function handle to the objective function that returns a numeric value
//' @param lst_OtherData A list of other data that may be passed to the fitness function
//' @param x_Min A NumericVector of lower bounds for each parameter
//' @param x_Max A NumericVector of upper bounds for each parameter
//' @param x_Init An optional initial solution (NumericVector), if NULL, the midpoint of `x_Min` and `x_Max` is used
//' @param max_iter The maximum number of iterations to perform (default is 100)
//' @param r The perturbation factor (default is 0.2)
//'
//' @return A List containing:
//'   - `x_Best`: The best parameter set found during the search
//'   - `y_Best`: The corresponding objective value of the best parameter set
//'
//' @export
// [[Rcpp::export]]
List cali_DDS(Function fitness,
              List lst_OtherData,
              NumericVector x_Min, NumericVector x_Max,
              Nullable<NumericVector> x_Init = R_NilValue, int max_iter = 100,
              double r = 0.2) {
  // Step 1: Set initial parameters
  NumericVector x_Best = x_Init.isNull() ? (x_Min + x_Max) / 2 : as<NumericVector>(x_Init);
  int n_x = x_Min.size();

  // Evaluate the initial objective function
  double y_Best = as<double>(fitness(x_Best, lst_OtherData));

  // Precompute P_i for perturbation probabilities
  NumericVector P_i(max_iter);
  for (int i = 0; i < max_iter; i++) {
    P_i[i] = 1.0 - std::log(i + 1.0) / std::log(max_iter);
  }

  // Random number generator
  std::random_device rd;
  std::default_random_engine generator(rd());

  // Create lst_Cali_x: perturbation indices for each iteration
  std::vector<std::vector<int>> lst_Cali_x(max_iter);
  for (int i = 0; i < max_iter; i++) {
    std::vector<bool> randm_Para = generate_binary_vector(n_x, P_i[i], generator);
    std::vector<int> selected_indices;
    for (int j = 0; j < n_x; j++) {
      if (randm_Para[j]) {
        selected_indices.push_back(j);
      }
    }
    // Ensure at least one index is perturbed
    if (selected_indices.empty()) {
      selected_indices.push_back(std::uniform_int_distribution<int>(0, n_x - 1)(generator));
    }
    lst_Cali_x[i] = selected_indices;
  }

  NumericVector sigma_ = x_Max - x_Min;

  // Progress bar simulation
  Rcout << "Calibration in progress...\n";

  for (int i = 1; i < max_iter; i++) {
    // Generate new candidate solution
    NumericVector x_New = clone(x_Best);
    NumericVector N_01 = generate_random_vector(n_x, 0.0, 1.0);

    NumericVector x_New0 = x_Best + r * N_01 * sigma_;
    NumericVector x_New1 = pmax(2 * x_Min - x_New0, x_Min);
    NumericVector x_New2 = pmin(2 * x_Max - x_New0, x_Max);
    for (int j = 0; j < n_x; j++) {
      if (x_New0[j] < x_Min[j]) x_New0[j] = x_New1[j];
      if (x_New0[j] > x_Max[j]) x_New0[j] = x_New2[j];
    }

    // Apply perturbation to selected indices
    for (int idx : lst_Cali_x[i]) {
      x_New[idx] = x_New0[idx];
    }

    // Evaluate the objective function
    double y_New = as<double>(fitness(x_New, lst_OtherData));
    if (y_New < y_Best) {
      x_Best = x_New;
      y_Best = y_New;
    }

    // Print progress
    if (i % (max_iter / 10) == 0) {
      Rcout << "Iteration " << i << " / " << max_iter << " completed.\n";
    }
  }

  // Return the results
  return List::create(Named("x_Best") = x_Best,
                      Named("y_Best") = y_Best);
}




