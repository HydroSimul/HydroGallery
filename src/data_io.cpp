#include <Rcpp.h>
#include <fstream>
#include <vector>

using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

// Matrix type enumeration
enum DataType { TYPE_NUMERIC = 1, TYPE_INTEGER = 2, TYPE_LOGICAL = 3 };

// Helper function to determine matrix type
DataType get_data_type(SEXP mat) {
switch(TYPEOF(mat)) {
case REALSXP:  return TYPE_NUMERIC;
case INTSXP:   return (Rf_isLogical(mat)) ? TYPE_LOGICAL : TYPE_INTEGER;
case LGLSXP:   return TYPE_LOGICAL;
default:       stop("Unsupported matrix type");
}
}

//' Save a matrix (Numeric, Integer, or Logical) to a binary file
//'
//' This function saves a matrix to a binary file with a header indicating
//' the matrix type (1=numeric, 2=integer, 3=logical), dimensions, and data.
//' All types are stored using 4 bytes per element.
//' @name dataio
//' @param matrix The matrix to be saved
//' @param filename The name of the binary file
//' @export
// [[Rcpp::export]]
void save_matbin(SEXP matrix, const std::string& filename) {
 DataType mtype = get_data_type(matrix);
 int rows, cols;
 std::ofstream file(filename, std::ios::binary);
 
 if (!file.is_open()) {
   stop("Failed to open file for writing: " + filename);
 }
 
 // Write matrix type code (1, 2, or 3)
 file.write(reinterpret_cast<const char*>(&mtype), sizeof(int));
 
 if (mtype == TYPE_NUMERIC) {
   NumericMatrix mat(matrix);
   rows = mat.nrow();
   cols = mat.ncol();
   
   // Write dimensions
   file.write(reinterpret_cast<const char*>(&rows), sizeof(int));
   file.write(reinterpret_cast<const char*>(&cols), sizeof(int));
   
   // Write data as floats (4 bytes per element)
   std::vector<float> buffer(rows * cols);
   for (int i = 0; i < rows; ++i) {
     for (int j = 0; j < cols; ++j) {
       buffer[j + i * cols] = static_cast<float>(mat(i, j));
     }
   }
   file.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(float));
 } 
 else if (mtype == TYPE_INTEGER) {
   IntegerMatrix mat(matrix);
   rows = mat.nrow();
   cols = mat.ncol();
   
   // Write dimensions
   file.write(reinterpret_cast<const char*>(&rows), sizeof(int));
   file.write(reinterpret_cast<const char*>(&cols), sizeof(int));
   
   // Write data as int (4 bytes per element)
   std::vector<int> buffer(rows * cols);
   for (int i = 0; i < rows; ++i) {
     for (int j = 0; j < cols; ++j) {
       buffer[j + i * cols] = static_cast<int>(mat(i, j));
     }
   }
   file.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(int));
 }
 else { // TYPE_LOGICAL
   LogicalMatrix mat(matrix);
   rows = mat.nrow();
   cols = mat.ncol();
   
   // Write dimensions
   file.write(reinterpret_cast<const char*>(&rows), sizeof(int));
   file.write(reinterpret_cast<const char*>(&cols), sizeof(int));
   
   // Write data as int (4 bytes per element, 0=FALSE, 1=TRUE)
   std::vector<int> buffer(rows * cols);
   for (int i = 0; i < rows; ++i) {
     for (int j = 0; j < cols; ++j) {
       buffer[j + i * cols] = static_cast<int>(mat(i, j));
     }
   }
   file.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(int));
 }
 
 file.close();
}

//' Load a matrix from a binary file
//'
//' This function loads a matrix from a binary file that was saved
//' using save_matbin(). It automatically detects the matrix type.
//' @rdname dataio
//' @param filename The binary file to load
//' @return Either a NumericMatrix, IntegerMatrix, or LogicalMatrix
//' @export
// [[Rcpp::export]]
SEXP load_matbin(const std::string& filename) {
 int type_code, rows, cols;
 std::ifstream file(filename, std::ios::binary);
 
 if (!file.is_open()) {
   stop("Failed to open file for reading: " + filename);
 }
 
 // Read matrix type code
 file.read(reinterpret_cast<char*>(&type_code), sizeof(int));
 
 // Read dimensions
 file.read(reinterpret_cast<char*>(&rows), sizeof(int));
 file.read(reinterpret_cast<char*>(&cols), sizeof(int));
 
 if (type_code == TYPE_NUMERIC) {
   std::vector<float> buffer(rows * cols);
   file.read(reinterpret_cast<char*>(buffer.data()), buffer.size() * sizeof(float));
   file.close();
   
   NumericMatrix mat(rows, cols);
   for (int i = 0; i < rows; ++i) {
     for (int j = 0; j < cols; ++j) {
       mat(i, j) = static_cast<double>(buffer[j + i * cols]);
     }
   }
   return mat;
 } 
 else if (type_code == TYPE_INTEGER) {
   std::vector<int> buffer(rows * cols);
   file.read(reinterpret_cast<char*>(buffer.data()), buffer.size() * sizeof(int));
   file.close();
   
   IntegerMatrix mat(rows, cols);
   for (int i = 0; i < rows; ++i) {
     for (int j = 0; j < cols; ++j) {
       mat(i, j) = static_cast<int>(buffer[j + i * cols]);
     }
   }
   return mat;
 }
 else if (type_code == TYPE_LOGICAL) {
   std::vector<int> buffer(rows * cols);
   file.read(reinterpret_cast<char*>(buffer.data()), buffer.size() * sizeof(int));
   file.close();
   
   LogicalMatrix mat(rows, cols);
   for (int i = 0; i < rows; ++i) {
     for (int j = 0; j < cols; ++j) {
       mat(i, j) = static_cast<bool>(buffer[j + i * cols]);
     }
   }
   return mat;
 }
 else {
   file.close();
   stop("Unknown matrix type code in file: " + std::to_string(type_code));
 }
}

//' Bind multiple binary matrix files into a single file
//'
//' This function combines several matrix files into one. All matrices
//' must be of the same type and have matching column counts.
//' @rdname dataio
//' @param input_files Vector of input file paths
//' @param output_file Path for the combined output file
//' @export
// [[Rcpp::export]]
void bind_matbin(const StringVector& input_files, const std::string& output_file) {
 if (input_files.size() == 0) {
   stop("No input files provided");
 }
 
 int total_rows = 0;
 int total_cols = 0;
 int expected_cols = 0;
 int expected_type = 0;
 bool first_file = true;
 
 // First pass: validate all files
 for (int f = 0; f < input_files.size(); f++) {
   std::string filename = as<std::string>(input_files[f]);
   std::ifstream file(filename, std::ios::binary);
   
   if (!file.is_open()) {
     stop("Failed to open file: " + filename);
   }
   
   int type_code, rows, cols;
   file.read(reinterpret_cast<char*>(&type_code), sizeof(int));
   file.read(reinterpret_cast<char*>(&rows), sizeof(int));
   file.read(reinterpret_cast<char*>(&cols), sizeof(int));
   file.close();
   
   if (first_file) {
     expected_type = type_code;
     expected_cols = cols;
     first_file = false;
   } else {
     if (type_code != expected_type) {
       stop("Type mismatch in file " + filename + 
         ": expected type " + std::to_string(expected_type) +
         ", got " + std::to_string(type_code));
     }
     if (cols != expected_cols) {
       stop("Column count mismatch in file " + filename + 
         ": expected " + std::to_string(expected_cols) +
         " columns, got " + std::to_string(cols));
     }
   }
   
   total_rows += rows;
   total_cols = expected_cols;
 }
 
 // Create output file
 std::ofstream outfile(output_file, std::ios::binary);
 if (!outfile.is_open()) {
   stop("Failed to create output file: " + output_file);
 }
 
 // Write header
 outfile.write(reinterpret_cast<const char*>(&expected_type), sizeof(int));
 outfile.write(reinterpret_cast<const char*>(&total_rows), sizeof(int));
 outfile.write(reinterpret_cast<const char*>(&total_cols), sizeof(int));
 
 // Concatenate data from all files
 for (int f = 0; f < input_files.size(); f++) {
   std::string filename = as<std::string>(input_files[f]);
   std::ifstream file(filename, std::ios::binary);
   
   int type_code, rows, cols;
   file.read(reinterpret_cast<char*>(&type_code), sizeof(int));
   file.read(reinterpret_cast<char*>(&rows), sizeof(int));
   file.read(reinterpret_cast<char*>(&cols), sizeof(int));
   
   // Skip header and copy data
   if (type_code == TYPE_NUMERIC) {
     std::vector<float> buffer(rows * cols);
     file.read(reinterpret_cast<char*>(buffer.data()), buffer.size() * sizeof(float));
     outfile.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(float));
   } 
   else { // TYPE_INTEGER or TYPE_LOGICAL both use int
     std::vector<int> buffer(rows * cols);
     file.read(reinterpret_cast<char*>(buffer.data()), buffer.size() * sizeof(int));
     outfile.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(int));
   }
   
   file.close();
 }
 
 outfile.close();
}


//' Save a vector (Numeric, Integer, or Logical) to a binary file
//'
//' This function saves a vector to a binary file with a header indicating
//' the vector type (1=numeric, 2=integer, 3=logical) and length.
//' All types are stored using 4 bytes per element.
//' @name dataio
//' @param vector The vector to be saved
//' @param filename The name of the binary file
//' @export
// [[Rcpp::export]]
void save_vecbin(SEXP vector, const std::string& filename) {
 DataType vtype = get_data_type(vector);
 int length;
 std::ofstream file(filename, std::ios::binary);
 
 if (!file.is_open()) {
   stop("Failed to open file for writing: " + filename);
 }
 
 // Write vector type code (1, 2, or 3)
 file.write(reinterpret_cast<const char*>(&vtype), sizeof(int));
 
 if (vtype == TYPE_NUMERIC) {
   NumericVector vec(vector);
   length = vec.size();
   
   // Write length
   file.write(reinterpret_cast<const char*>(&length), sizeof(int));
   
   // Write data as floats (4 bytes per element)
   std::vector<float> buffer(length);
   for (int i = 0; i < length; ++i) {
     buffer[i] = static_cast<float>(vec[i]);
   }
   file.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(float));
 } 
 else if (vtype == TYPE_INTEGER) {
   IntegerVector vec(vector);
   length = vec.size();
   
   // Write length
   file.write(reinterpret_cast<const char*>(&length), sizeof(int));
   
   // Write data as int (4 bytes per element)
   std::vector<int> buffer(length);
   for (int i = 0; i < length; ++i) {
     buffer[i] = static_cast<int>(vec[i]);
   }
   file.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(int));
 }
 else { // TYPE_LOGICAL
   LogicalVector vec(vector);
   length = vec.size();
   
   // Write length
   file.write(reinterpret_cast<const char*>(&length), sizeof(int));
   
   // Write data as int (4 bytes per element, 0=FALSE, 1=TRUE)
   std::vector<int> buffer(length);
   for (int i = 0; i < length; ++i) {
     buffer[i] = static_cast<int>(vec[i]);
   }
   file.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(int));
 }
 
 file.close();
}

//' Load a vector from a binary file
//'
//' This function loads a vector from a binary file that was saved
//' using save_vecbin(). It automatically detects the vector type.
//' @rdname dataio
//' @param filename The binary file to load
//' @return Either a NumericVector, IntegerVector, or LogicalVector
//' @export
// [[Rcpp::export]]
SEXP load_vecbin(const std::string& filename) {
 int type_code, length;
 std::ifstream file(filename, std::ios::binary);
 
 if (!file.is_open()) {
   stop("Failed to open file for reading: " + filename);
 }
 
 // Read vector type code
 file.read(reinterpret_cast<char*>(&type_code), sizeof(int));
 
 // Read length
 file.read(reinterpret_cast<char*>(&length), sizeof(int));
 
 if (type_code == TYPE_NUMERIC) {
   std::vector<float> buffer(length);
   file.read(reinterpret_cast<char*>(buffer.data()), buffer.size() * sizeof(float));
   file.close();
   
   NumericVector vec(length);
   for (int i = 0; i < length; ++i) {
     vec[i] = static_cast<double>(buffer[i]);
   }
   return vec;
 } 
 else if (type_code == TYPE_INTEGER) {
   std::vector<int> buffer(length);
   file.read(reinterpret_cast<char*>(buffer.data()), buffer.size() * sizeof(int));
   file.close();
   
   IntegerVector vec(length);
   for (int i = 0; i < length; ++i) {
     vec[i] = static_cast<int>(buffer[i]);
   }
   return vec;
 }
 else if (type_code == TYPE_LOGICAL) {
   std::vector<int> buffer(length);
   file.read(reinterpret_cast<char*>(buffer.data()), buffer.size() * sizeof(int));
   file.close();
   
   LogicalVector vec(length);
   for (int i = 0; i < length; ++i) {
     vec[i] = static_cast<bool>(buffer[i]);
   }
   return vec;
 }
 else {
   file.close();
   stop("Unknown vector type code in file: " + std::to_string(type_code));
 }
}

