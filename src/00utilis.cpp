#include "00utilis.h"

// [[Rcpp::interfaces(r, cpp)]]

//' **utilize functions**
//' @name utilis
NumericVector vecpow(NumericVector base, NumericVector exp) {
  NumericVector out(base.size());
  std::transform(base.begin(), base.end(),
                 exp.begin(), out.begin(), ::pow);
  return (out);
}


NumericVector vec_const_pow(NumericVector base, double exp) {
  NumericVector out(base.size());
  std::transform(base.begin(), base.end(), out.begin(), 
                 [exp](double b) { return std::pow(b, exp); });
  return out;
}


NumericVector vecpow10(NumericVector exp) {

  NumericVector base10 (exp.size(), 10.0);

  return (vecpow(base10, exp));
}

double sum_product(NumericVector lhs,
                   NumericVector rhs)
{
  int n = lhs.size();
  double s = 0.0;
  for (int i= 0; i < n; i++) {
    s += lhs(i) * rhs(i);
  }
  return (s);
}

void resetVector(Rcpp::NumericVector& x) {
  // Fill the vector with zeros
  std::fill(x.begin(), x.end(), 0.0);
}

template <typename T>
T subset_get_(T data, IntegerVector int_Index) {
  int n = int_Index.size();
  T result(n);  // Create a result vector of the same type as data

  for (int i = 0; i < n; ++i) {
    result[i] = data[int_Index[i] - 1];  // Adjust for R's 1-based indexing
  }

  return result;
}

//' @rdname utilis
//' @param vec_Data (vector of num / int) data
//' @param int_Index (vector of int) index of subset
//' @export
// [[Rcpp::export]]
NumericVector subset_get(NumericVector vec_Data, IntegerVector int_Index) {
  return subset_get_(vec_Data, int_Index);
}

//' @rdname utilis
//' @export
// [[Rcpp::export]]
LogicalVector subset_get_logical(LogicalVector vec_Data, IntegerVector int_Index) {
  return subset_get_(vec_Data, int_Index);
}


//' @rdname utilis
//' @param vec_DataPut (vector of num / int) data, to refresh
//' @export
// [[Rcpp::export]]
void subset_put(NumericVector &vec_Data, IntegerVector int_Index, NumericVector vec_DataPut) {
  int n = int_Index.size();

  for (int i = 0; i < n; ++i) {
    vec_Data[int_Index[i] - 1] = vec_DataPut[i];
  }

}

//' @rdname utilis
//' @param vec_DataAdd (vector of num / int) data, to refresh
//' @export
// [[Rcpp::export]]
void subset_add(NumericVector &vec_Data, IntegerVector int_Index, NumericVector vec_DataAdd) {
 int n = int_Index.size();

 for (int i = 0; i < n; ++i) {
   vec_Data[int_Index[i] - 1] += vec_DataAdd[i];
 }

}




// Template function for concatenating two vectors
template <typename T>
T concat_vectors(const T& vec_A, const T& vec_B) {
  int n_A = vec_A.size();
  int n_B = vec_B.size();
  T result(n_A + n_B);

  // Use std::copy to copy the contents of vec_A and vec_B into result
  std::copy(vec_A.begin(), vec_A.end(), result.begin());
  std::copy(vec_B.begin(), vec_B.end(), result.begin() + n_A);

  return result;
}

IntegerVector c_int(IntegerVector vec_A, IntegerVector vec_B) {
  return concat_vectors(vec_A, vec_B);
}

NumericVector c_num(NumericVector vec_A, NumericVector vec_B) {
  return concat_vectors(vec_A, vec_B);
}



IntegerVector find_locations(IntegerVector x, int x0) {
  std::vector<int> locations;

  for (int i = 0; i < x.size(); ++i) {
    if (x[i] == x0) {
      locations.push_back(i + 1); // Use 1-based indexing for R compatibility
    }
  }

  return wrap(locations);
}


//' @rdname utilis
//' @param mat_X,mat_Y (matrix of num / int) data
//' @param num_X0 (vector of num) data
//' @export
// [[Rcpp::export]]
NumericVector linear_interpolate_vec(const NumericMatrix& mat_X, const NumericMatrix& mat_Y, const NumericVector& num_X0) {
  int n_Vec = mat_X.ncol();  // Now vectors are along columns
  int n_Interpolat = mat_X.nrow(); // Interpolation points along rows
  NumericVector num_Y0(n_Vec);
  
  for (int i = 0; i < n_Vec; ++i) {
    double x0 = num_X0[i];
    std::vector<double> X(n_Interpolat);
    std::vector<double> Y(n_Interpolat);
    
    for (int j = 0; j < n_Interpolat; ++j) {
      X[j] = mat_X(j, i);  // Access column-wise
      Y[j] = mat_Y(j, i);  // Access column-wise
    }
    
    int idx = std::lower_bound(X.begin(), X.end(), x0) - X.begin();
    if (idx == 0) {
      num_Y0[i] = Y[0];
    } else if (idx == n_Interpolat) {
      num_Y0[i] = Y[n_Interpolat - 1];
    } else {
      double x1 = X[idx - 1], x2 = X[idx];
      double y1 = Y[idx - 1], y2 = Y[idx];
      num_Y0[i] = y1 + (x0 - x1) * (y2 - y1) / (x2 - x1);
    }
  }
  
  return num_Y0;
}
