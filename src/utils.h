#ifndef UTILS
#define UTILS

#include <RcppArmadillo.h>
#include <string>
#include <set>
#include <algorithm>

arma::mat load_mat(const std::string& path);
arma::vec load_vec(const std::string& path);

#endif // UTILS
