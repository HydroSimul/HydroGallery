#include <RcppArmadillo.h>
#include <unordered_map>
#include <set>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

//' @rdname routingtopology
//' @name get_inflow_cells
//' @title Get Inflow Cells
//' @description This function calculates inflow cells based on the outflow vector.
//' @param int_Outflow (vector of int) The cell number of the next cell. The cell number must range from 1 to the length of the cells.
//' If the cell has no outflow, the number should be set to itself.
//' @return A list where each element is an IntegerVector containing the inflow cells for each cell.
//' @export
// [[Rcpp::export]]
std::vector<arma::ivec> get_inflow_cells(const arma::ivec& int_Outflow) {
  int n = int_Outflow.n_elem;
  std::vector<std::vector<int>> temp(n);
  
  for (int i = 0; i < n; ++i) {
    int origin = i + 1;
    int next = int_Outflow(i);
    temp[i].push_back(origin);
    
    while (next != origin) {
      temp[next - 1].push_back(origin);
      origin = next;
      next = int_Outflow(origin - 1);
    }
  }
  
  std::vector<arma::ivec> result(n);
  for (int i = 0; i < n; ++i)
    result[i] = arma::ivec(temp[i]);
  
  return result;
}

//' @rdname routingtopology
//' @name get_inflow_lastcell
//' @title Get Inflow Last Cell Matrix
//' @description This function creates a matrix of inflow cells for each cell based on the outflow vector.
//' @return A NumericMatrix where each row corresponds to a cell, and each column represents an inflow cell.
//' @export
// [[Rcpp::export]]
arma::mat get_inflow_lastcell(const arma::ivec& int_Outflow) {
  int n = int_Outflow.n_elem;
  std::vector<std::vector<double>> temp(n);
  std::size_t max_size = 0;
  
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (int_Outflow(j) == i + 1) {
        temp[i].push_back(j + 1);
      }
    }
    if (temp[i].size() > max_size)
      max_size = temp[i].size();
  }
  
  arma::mat out(n, max_size, arma::fill::value(NA_REAL));
  for (int i = 0; i < n; ++i)
    for (std::size_t j = 0; j < temp[i].size(); ++j)
      out(i, j) = temp[i][j];
  
  return out;
}

std::vector<arma::ivec> get_step_cells(const std::vector<arma::ivec>& inflow_cells) {
  int n = inflow_cells.size();
  arma::ivec lengths(n);
  for (int i = 0; i < n; ++i)
    lengths(i) = inflow_cells[i].n_elem;
  
  std::set<int> unique_lengths(lengths.begin(), lengths.end());
  std::vector<int> sorted_lengths(unique_lengths.begin(), unique_lengths.end());
  
  std::unordered_map<int, int> length_to_step;
  for (std::size_t i = 0; i < sorted_lengths.size(); ++i)
    length_to_step[sorted_lengths[i]] = i + 1;
  
  std::vector<std::vector<int>> step_groups(sorted_lengths.size());
  for (int i = 0; i < n; ++i)
    step_groups[length_to_step[lengths(i)] - 1].push_back(i + 1);
  
  std::vector<arma::ivec> result(sorted_lengths.size());
  for (std::size_t i = 0; i < step_groups.size(); ++i)
    result[i] = arma::ivec(step_groups[i]);
  
  return result;
}

std::vector<arma::mat> get_step_lastcell(const std::vector<arma::ivec>& step_cells,
                                              const arma::mat& inflow_lastcell) {
  std::vector<arma::mat> result(step_cells.size());
  result[0].reset();  // First step is empty
  
  for (std::size_t i = 1; i < step_cells.size(); ++i) {
    arma::uvec idx = arma::conv_to<arma::uvec>::from(step_cells[i]) - 1;
    result[i] = inflow_lastcell.rows(idx);
  }
  
  return result;
}

//' @rdname routingtopology
//' @param lst_Inflow_Cell A vector of arma::ivec, where each arma::ivec contains the cells that flow into the respective cell.
//' @param int_OutLet An integer representing the outlet cell (1-based index).
//' @param int_TestCell An integer vector, cells to test.
//' @return An integer vector of cells in the intersection of the station cells and the basin.
//' @export
// [[Rcpp::export]]
arma::ivec get_cell_in_basin(const std::vector<arma::ivec>& lst_Inflow_Cell, int int_OutLet, const arma::ivec& int_TestCell) {
  
  // Extract the Big Basin (1-based indexing)
  arma::ivec int_BigBasin = lst_Inflow_Cell[int_OutLet - 1];

  // Create a set from int_BigBasin for fast lookup
  std::set<int> big_basin_set(int_BigBasin.begin(), int_BigBasin.end());

  // Remove int_OutLet from int_TestCell if present
  std::vector<int> int_TestCell_no_outlet;
  for (size_t i = 0; i < int_TestCell.n_elem; ++i) {
    if (int_TestCell[i] != int_OutLet) {
      int_TestCell_no_outlet.push_back(int_TestCell[i]);
    }
  }

  // Convert to std::set for fast lookup in intersection
  std::set<int> station_cells(int_TestCell_no_outlet.begin(), int_TestCell_no_outlet.end());

  // Find the intersection
  std::vector<int> intersection;
  std::set_intersection(
    big_basin_set.begin(), big_basin_set.end(),
    station_cells.begin(), station_cells.end(),
    std::back_inserter(intersection)
  );

  // Return the intersection as an arma::ivec
  return arma::ivec(intersection);
}

//' @rdname routingtopology
//' @param int_UpstreamCell An integer vector containing the upstream cells to find the upstream basin.
//' @return An integer vector representing the new upstream basin, which includes the upstream cells and the set difference of the basin cells.
//' This function identifies the upstream basin of a given outlet cell by first finding the intersection of the upstream cells
//' with the cells that flow into the outlet. It then computes the set difference between the upstream basin and the outlet basin.
//' @export
// [[Rcpp::export]]
arma::ivec get_inter_basin(const arma::ivec& int_Cell, const arma::ivec& int_Outflow) {
  int n = int_Cell.n_elem;
  arma::ivec out(n);
  
  for (int i = 0; i < n; ++i) {
    int id = int_Cell[i] - 1;
    int next = int_Outflow[id];
    out[i] = (std::find(int_Cell.begin(), int_Cell.end(), next) == int_Cell.end()) ? next : NA_INTEGER;
  }
  
  return out;
}

//' @rdname routingtopology
//' @param int_Outflow_Ori An integer vector representing the original outflow indices (1-based).
//' @param int_CellNew An integer vector representing the cells within the new basin.
//' @return An integer vector of the new outflow indices adjusted for the sub-basin.
//' @export
// [[Rcpp::export]]
arma::ivec get_new_outflow(const arma::ivec& int_Cell, const arma::ivec& int_Outflow) {
  int n = int_Cell.n_elem;
  arma::ivec result(n, arma::fill::value(NA_INTEGER));
  std::unordered_map<int, int> old_to_new;
  
  for (int i = 0; i < n; ++i)
    old_to_new[int_Cell[i]] = i + 1;
  
  for (int i = 0; i < n; ++i) {
    int id = int_Cell[i];
    int next = int_Outflow[id - 1];
    result(i) = old_to_new.count(next) ? old_to_new[next] : id;
  }
  
  return result;
}

//' @rdname routingtopology
//' @param int_CaliCell An integer vector of calibration cells.
//' @return A list of integer vectors (`lst_Step_Cali`), where each element represents calibration cells at a specific step.
//' @export
// [[Rcpp::export]]
arma::ivec get_cali_step(const std::vector<arma::ivec>& step_cells,
                              const arma::ivec& int_Cali) {
  std::vector<int> steps;
  for (std::size_t i = 0; i < step_cells.size(); ++i) {
    for (int id : step_cells[i]) {
      if (arma::any(int_Cali == id)) {
        steps.push_back(i + 1);
      }
    }
  }
  
  return arma::ivec(steps);
}

//' @rdname routingtopology
//' @title Get Step Parameters
//' @description This function returns a list of step cells and the corresponding last cell matrices.
//' @return A list containing the step cells and last cell matrices.
//' @export
// [[Rcpp::export]]
Rcpp::List get_routing_info(const arma::ivec& int_Outflow) {
  std::vector<arma::ivec> inflow_cells = get_inflow_cells(int_Outflow);
  std::vector<arma::ivec> step_cells = get_step_cells(inflow_cells);
  arma::mat inflow_lastcell = get_inflow_lastcell(int_Outflow);
  std::vector<arma::mat> step_lastcell = get_step_lastcell(step_cells, inflow_lastcell);
  
  return Rcpp::List::create(
    Rcpp::Named("int_Cell") = step_cells,
    Rcpp::Named("mat_LastCell") = step_lastcell
  );
}

//' @rdname routingtopology
//' @return A list of integer vectors (`lst_Step_Cali`), where each element represents calibration cells at a specific step.
//' @export
// [[Rcpp::export]]
std::vector<arma::ivec> get_upstream_cali_cell(const std::vector<arma::ivec>& lst_Inflow_Cell, const arma::ivec& int_CaliCell) {
  int n_CaliCells = int_CaliCell.n_elem;

  // Step 1: Get upstream cells for each calibration cell
  std::vector<arma::ivec> lst_Cali_Upstream(n_CaliCells);
  for (int i = 0; i < n_CaliCells; ++i) {
    // Get upstream cells for the current calibration cell
    lst_Cali_Upstream[i] = get_cell_in_basin(lst_Inflow_Cell, int_CaliCell(i), int_CaliCell);
  }

  // Step 2: Determine the last calibration cell for each calibration cell
  std::vector<arma::ivec> lst_LastCaliCell(n_CaliCells);
  for (int i = 0; i < n_CaliCells; ++i) {
    arma::ivec upstream_cells = lst_Cali_Upstream[i];

    // Map upstream cells to their indices in int_CaliCell
    std::unordered_set<int> upstream_indices;
    for (int j = 0; j < upstream_cells.n_elem; ++j) {
      auto it = std::find(int_CaliCell.begin(), int_CaliCell.end(), upstream_cells(j));
      if (it != int_CaliCell.end()) {
        upstream_indices.insert(it - int_CaliCell.begin());
      }
    }

    // Collect last calibration cells by removing the common upstream cells
    std::unordered_set<int> temp_set;
    for (int index : upstream_indices) {
      arma::ivec temp = lst_Cali_Upstream[index];
      temp_set.insert(temp.begin(), temp.end());
    }

    std::vector<int> result;
    for (int cell : upstream_cells) {
      if (temp_set.find(cell) == temp_set.end()) {
        result.push_back(cell);
      }
    }

    lst_LastCaliCell[i] = arma::conv_to<arma::ivec>::from(result);
  }

  return lst_LastCaliCell;
}



//' @rdname routingtopology
//' @export
// [[Rcpp::export]]
void write_int_vector_list(const std::vector<arma::ivec>& vec_list, const std::string& file_path) {
  std::ofstream fout(file_path, std::ios::binary);
  if (!fout) Rcpp::stop("Cannot open file for writing");
  
  int32_t num_vectors = vec_list.size();
  fout.write(reinterpret_cast<char*>(&num_vectors), sizeof(int32_t));
  
  for (const auto& vec : vec_list) {
    int32_t len = vec.n_elem;
    fout.write(reinterpret_cast<char*>(&len), sizeof(int32_t));
    fout.write(reinterpret_cast<const char*>(vec.memptr()), len * sizeof(int32_t));
  }
  
  fout.close();
}

//' @rdname routingtopology
//' @export
// [[Rcpp::export]]
std::vector<arma::ivec> read_int_vector_list(const std::string& file_path) {
  std::ifstream fin(file_path, std::ios::binary);
  if (!fin) Rcpp::stop("Cannot open file for reading");
  
  int32_t num_vectors;
  fin.read(reinterpret_cast<char*>(&num_vectors), sizeof(int32_t));
  if (fin.eof() || num_vectors < 0) Rcpp::stop("Invalid or corrupted file");
  
  std::vector<arma::ivec> result;
  result.reserve(num_vectors);
  
  for (int i = 0; i < num_vectors; ++i) {
    int32_t len;
    fin.read(reinterpret_cast<char*>(&len), sizeof(int32_t));
    if (fin.eof() || len < 0) Rcpp::stop("Invalid vector length");
    
    arma::ivec vec(len);
    fin.read(reinterpret_cast<char*>(vec.memptr()), len * sizeof(int32_t));
    if (fin.eof()) Rcpp::stop("Unexpected end of file");
    
    result.push_back(std::move(vec));
  }
  
  fin.close();
  return result;
}

//' @rdname routingtopology
//' @export
// [[Rcpp::export]]
void write_int_matrix_list(const std::vector<arma::imat>& mat_list, const std::string& file_path) {
  std::ofstream fout(file_path, std::ios::binary);
  if (!fout) Rcpp::stop("Cannot open file for writing");
  
  int32_t num_matrices = mat_list.size();
  fout.write(reinterpret_cast<char*>(&num_matrices), sizeof(int32_t));
  
  for (const auto& mat : mat_list) {
    if (mat.n_cols != 9) Rcpp::stop("Matrix must have 9 columns");
    
    int32_t nrow = mat.n_rows;
    fout.write(reinterpret_cast<char*>(&nrow), sizeof(int32_t));
    fout.write(reinterpret_cast<const char*>(mat.memptr()), nrow * 9 * sizeof(int32_t));
  }
  
  fout.close();
}

//' @rdname routingtopology
//' @export
// [[Rcpp::export]]
std::vector<arma::imat> read_int_matrix_list(const std::string& file_path) {
  std::ifstream fin(file_path, std::ios::binary);
  if (!fin) Rcpp::stop("Cannot open file for reading");
  
  int32_t num_matrices;
  fin.read(reinterpret_cast<char*>(&num_matrices), sizeof(int32_t));
  if (fin.eof() || num_matrices < 0) Rcpp::stop("Corrupted or empty file");
  
  std::vector<arma::imat> result;
  result.reserve(num_matrices);
  
  for (int i = 0; i < num_matrices; ++i) {
    int32_t nrow;
    fin.read(reinterpret_cast<char*>(&nrow), sizeof(int32_t));
    if (fin.eof() || nrow < 0) Rcpp::stop("Corrupted file: invalid row count");
    
    arma::imat mat(nrow, 9);
    fin.read(reinterpret_cast<char*>(mat.memptr()), nrow * 9 * sizeof(int32_t));
    if (fin.eof()) Rcpp::stop("Unexpected end of file");
    
    result.push_back(std::move(mat));
  }
  
  fin.close();
  return result;
}
