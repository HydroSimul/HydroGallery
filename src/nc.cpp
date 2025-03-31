#include <Rcpp.h>
#include <netcdf.h>
#include <string>
#include <vector>
#include <stdexcept>

using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

//' NetCDF Interface Functions
 //' 
 //' @name nc
 //' @description Functions for reading and writing 2D data from/to NetCDF files
 //' 
 //' @param fileName Path to the NetCDF file
 //' @param varName Name of the variable to read or write
 //' @param mat_Data NumericMatrix containing the data to write
 //' @param dimName1 Name for the first dimension
 //' @param dimName2 Name for the second dimension
 //' @param varAttrs Named list of variable attributes (optional)
 //' @param globalAttrs Named list of global attributes (optional)
 //' 
 //' @details
 //' These functions provide a simple interface for NetCDF operations in R,
 //' allowing reading and writing of 2D matrices to NetCDF files.
 //' 
 //' 
 //' @export
 // [[Rcpp::export]]
 NumericMatrix read_nc_2d(const std::string& fileName, const std::string& varName) {
   int ncid, status;
   
   // Open the NetCDF file
   status = nc_open(fileName.c_str(), NC_NOWRITE, &ncid);
   if (status != NC_NOERR) {
     throw std::runtime_error("Failed to open NetCDF file: " + 
                              std::string(nc_strerror(status)));
   }
   
   // Use RAII to ensure file gets closed even if exceptions occur
   struct NCFileCloser {
     int ncid;
     NCFileCloser(int id) : ncid(id) {}
     ~NCFileCloser() { nc_close(ncid); }
   } fileCloser(ncid);
   
   // Get the variable ID
   int varid;
   status = nc_inq_varid(ncid, varName.c_str(), &varid);
   if (status != NC_NOERR) {
     throw std::runtime_error("Failed to find variable '" + varName + 
                              "': " + std::string(nc_strerror(status)));
   }
   
   // Check dimensions
   int ndims;
   status = nc_inq_varndims(ncid, varid, &ndims);
   if (status != NC_NOERR) {
     throw std::runtime_error("Failed to get number of dimensions: " + 
                              std::string(nc_strerror(status)));
   }
   
   if (ndims != 2) {
     throw std::runtime_error("Variable '" + varName + 
                              "' is not 2D. It has " + std::to_string(ndims) + " dimensions.");
   }
   
   // Get dimension information
   int dimids[2];
   status = nc_inq_vardimid(ncid, varid, dimids);
   if (status != NC_NOERR) {
     throw std::runtime_error("Failed to get dimension IDs: " + 
                              std::string(nc_strerror(status)));
   }
   
   size_t dims[2];
   for (int i = 0; i < 2; i++) {
     status = nc_inq_dimlen(ncid, dimids[i], &dims[i]);
     if (status != NC_NOERR) {
       throw std::runtime_error("Failed to get dimension length: " + 
                                std::string(nc_strerror(status)));
     }
   }
   
   // Create output matrix with correct dimensions (R uses column-major order)
   NumericMatrix result(dims[1], dims[0]);
   
   // Allocate buffer for reading data in native layout
   std::vector<double> buffer(dims[0] * dims[1]);
   
   // Read the data
   status = nc_get_var_double(ncid, varid, buffer.data());
   if (status != NC_NOERR) {
     throw std::runtime_error("Failed to read variable data: " + 
                              std::string(nc_strerror(status)));
   }
   
   // Transfer data to R matrix with correct orientation
   for (size_t i = 0; i < dims[0]; ++i) {
     for (size_t j = 0; j < dims[1]; ++j) {
       result(j, i) = buffer[i * dims[1] + j];
     }
   }
   
   // File will be automatically closed by NCFileCloser destructor
   return result;
 }
 
 //' @rdname nc
 //' @export
 // [[Rcpp::export]]
 void write_nc_2d(NumericMatrix mat_Data, 
                  const std::string& fileName, 
                  const std::string& varName = "data", 
                  const std::string& dimName1 = "rows", 
                  const std::string& dimName2 = "cols",
                  Rcpp::Nullable<List> varAttrs = R_NilValue,
                  Rcpp::Nullable<List> globalAttrs = R_NilValue) {
   
   int ncid, status;
   int varid, dimids[2];
   int nrows = mat_Data.nrow();
   int ncols = mat_Data.ncol();
   
   // Create NetCDF file (NC_CLOBBER will overwrite existing file)
   status = nc_create(fileName.c_str(), NC_CLOBBER, &ncid);
   if (status != NC_NOERR) {
     throw std::runtime_error("Failed to create NetCDF file '" + fileName + 
                              "': " + std::string(nc_strerror(status)));
   }
   
   // Use RAII to ensure file gets closed even if exceptions occur
   struct NCFileCloser {
     int ncid;
     NCFileCloser(int id) : ncid(id) {}
     ~NCFileCloser() { nc_close(ncid); }
   } fileCloser(ncid);
   
   // Add global attributes if provided
   if (globalAttrs.isNotNull()) {
     List attrs(globalAttrs);
     CharacterVector attrNames = attrs.names();
     
     for (int i = 0; i < attrs.size(); i++) {
       std::string attrName = Rcpp::as<std::string>(attrNames[i]);
       
       // Handle different attribute types
       SEXP attrValue = attrs[i];
       switch (TYPEOF(attrValue)) {
       case REALSXP: {
         // Handle both single values and vectors
         if (Rf_length(attrValue) == 1) {
         double value = Rcpp::as<double>(attrValue);
         status = nc_put_att_double(ncid, NC_GLOBAL, attrName.c_str(), NC_DOUBLE, 1, &value);
       } else {
         Rcpp::NumericVector values = Rcpp::as<Rcpp::NumericVector>(attrValue);
         status = nc_put_att_double(ncid, NC_GLOBAL, attrName.c_str(), NC_DOUBLE, 
                                    values.size(), values.begin());
       }
       break;
       }
       case INTSXP: {
         // Handle both single values and vectors
         if (Rf_length(attrValue) == 1) {
         int value = Rcpp::as<int>(attrValue);
         status = nc_put_att_int(ncid, NC_GLOBAL, attrName.c_str(), NC_INT, 1, &value);
       } else {
         Rcpp::IntegerVector values = Rcpp::as<Rcpp::IntegerVector>(attrValue);
         status = nc_put_att_int(ncid, NC_GLOBAL, attrName.c_str(), NC_INT, 
                                 values.size(), values.begin());
       }
       break;
       }
       case STRSXP: {
         std::string value = Rcpp::as<std::string>(attrValue);
         status = nc_put_att_text(ncid, NC_GLOBAL, attrName.c_str(), value.length(), value.c_str());
         break;
       }
       default:
         Rcpp::warning("Skipping global attribute '%s' with unsupported type", attrName.c_str());
         continue;
       }
       
       if (status != NC_NOERR) {
         Rcpp::warning("Failed to add global attribute '%s': %s", 
                       attrName.c_str(), nc_strerror(status));
       }
     }
   }
   
   // Define dimensions (column dimension first, then row dimension)
   status = nc_def_dim(ncid, dimName2.c_str(), ncols, &dimids[0]);
   if (status != NC_NOERR) {
     throw std::runtime_error("Failed to define column dimension: " + 
                              std::string(nc_strerror(status)));
   }
   
   status = nc_def_dim(ncid, dimName1.c_str(), nrows, &dimids[1]);
   if (status != NC_NOERR) {
     throw std::runtime_error("Failed to define row dimension: " + 
                              std::string(nc_strerror(status)));
   }
   
   // Define the variable
   status = nc_def_var(ncid, varName.c_str(), NC_DOUBLE, 2, dimids, &varid);
   if (status != NC_NOERR) {
     throw std::runtime_error("Failed to define variable: " + 
                              std::string(nc_strerror(status)));
   }
   
   // Add variable attributes if provided
   if (varAttrs.isNotNull()) {
     List attrs(varAttrs);
     CharacterVector attrNames = attrs.names();
     
     for (int i = 0; i < attrs.size(); i++) {
       std::string attrName = Rcpp::as<std::string>(attrNames[i]);
       
       // Handle different attribute types
       SEXP attrValue = attrs[i];
       switch (TYPEOF(attrValue)) {
       case REALSXP: {
         // Handle both single values and vectors
         if (Rf_length(attrValue) == 1) {
         double value = Rcpp::as<double>(attrValue);
         status = nc_put_att_double(ncid, varid, attrName.c_str(), NC_DOUBLE, 1, &value);
       } else {
         Rcpp::NumericVector values = Rcpp::as<Rcpp::NumericVector>(attrValue);
         status = nc_put_att_double(ncid, varid, attrName.c_str(), NC_DOUBLE, 
                                    values.size(), values.begin());
       }
       break;
       }
       case INTSXP: {
         // Handle both single values and vectors
         if (Rf_length(attrValue) == 1) {
         int value = Rcpp::as<int>(attrValue);
         status = nc_put_att_int(ncid, varid, attrName.c_str(), NC_INT, 1, &value);
       } else {
         Rcpp::IntegerVector values = Rcpp::as<Rcpp::IntegerVector>(attrValue);
         status = nc_put_att_int(ncid, varid, attrName.c_str(), NC_INT, 
                                 values.size(), values.begin());
       }
       break;
       }
       case STRSXP: {
         std::string value = Rcpp::as<std::string>(attrValue);
         status = nc_put_att_text(ncid, varid, attrName.c_str(), value.length(), value.c_str());
         break;
       }
       default:
         Rcpp::warning("Skipping variable attribute '%s' with unsupported type", attrName.c_str());
         continue;
       }
       
       if (status != NC_NOERR) {
         Rcpp::warning("Failed to add variable attribute '%s': %s", 
                       attrName.c_str(), nc_strerror(status));
       }
     }
   }
   
   // End definition mode
   status = nc_enddef(ncid);
   if (status != NC_NOERR) {
     throw std::runtime_error("Failed to end definition mode: " + 
                              std::string(nc_strerror(status)));
   }
   
   // Copy data to buffer with proper layout transformation
   std::vector<double> buffer(nrows * ncols);
   for (int j = 0; j < ncols; ++j) {
     for (int i = 0; i < nrows; ++i) {
       buffer[j * nrows + i] = mat_Data(i, j);
     }
   }
   
   // Write the data
   status = nc_put_var_double(ncid, varid, buffer.data());
   if (status != NC_NOERR) {
     throw std::runtime_error("Failed to write data: " + 
                              std::string(nc_strerror(status)));
   }
   
   // File will be automatically closed by NCFileCloser destructor
   return;
 }