#include <Rcpp.h>
using namespace Rcpp;
#include <netcdf.h>
// [[Rcpp::interfaces(r, cpp)]]

//' Evalute matrics
//' @name nc
//' @param fileName NC-File name.
//' @param varName Variable Name.
//' @return 2D-Matrix.
//' @export
// [[Rcpp::export]]
NumericMatrix read_nc_2d(const std::string& fileName, const std::string& varName) {
  // Open the NetCDF file
  int ncid;
  if (nc_open(fileName.c_str(), NC_NOWRITE, &ncid)) {
    throw std::runtime_error("Error opening NetCDF file");
  }
  
  // Get the ID of the variable
  int varid;
  if (nc_inq_varid(ncid, varName.c_str(), &varid)) {
    throw std::runtime_error("Error finding variable in NetCDF file");
  }
  
  // Get the dimensions of the variable
  int ndim;
  if (nc_inq_varndims(ncid, varid, &ndim)) {
    throw std::runtime_error("Error getting number of dimensions for the variable");
  }
  
  if (ndim != 2) {
    throw std::runtime_error("Variable is not 2D.");
  }
  
  // Get the dimensions of the 2D variable
  size_t dims[2];
  int dimids[2];
  if (nc_inq_vardimid(ncid, varid, dimids)) {
    throw std::runtime_error("Error getting dimensions for the variable");
  }
  
  if (nc_inq_dimlen(ncid, dimids[0], &dims[0]) || nc_inq_dimlen(ncid, dimids[1], &dims[1])) {
    throw std::runtime_error("Error getting dimension lengths");
  }
  
  // Allocate a buffer for the data
  std::vector<float> data(dims[0] * dims[1]);
  
  // Read the data into the buffer
  if (nc_get_var_float(ncid, varid, data.data())) {
    throw std::runtime_error("Error reading variable data");
  }
  
  // // Create the NumericMatrix preserving the original NetCDF dimension ordering
  // NumericMatrix matrix(dims[0], dims[1]);
  // 
  // // Copy data from the buffer into the NumericMatrix
  // for (size_t i = 0; i < dims[0]; ++i) {
  //   for (size_t j = 0; j < dims[1]; ++j) {
  //     matrix(i, j) = data[i * dims[1] + j];
  //   }
  // }
  
  
  // Create the NumericMatrix with R's expected ordering (transposed)
  NumericMatrix matrix(dims[1], dims[0]); // Swap dimensions

  // Copy data from buffer with transposition
  for (size_t i = 0; i < dims[0]; ++i) {
    for (size_t j = 0; j < dims[1]; ++j) {
      matrix(j, i) = data[i * dims[1] + j]; // Swap i and j
    }
  }
  
  // Close the NetCDF file
  nc_close(ncid);
  
  return matrix;
}



