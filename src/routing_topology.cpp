#include <Rcpp.h>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <set>       // for std::set
using namespace Rcpp;
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
List get_inflow_cells(IntegerVector int_Outflow) {
  int n_Cell = int_Outflow.size();
  List lst_Inflow(n_Cell);

  // Initialize each element of lst_Inflow with an empty IntegerVector
  for (int i = 0; i < n_Cell; i++) {
    lst_Inflow[i] = IntegerVector();
  }

  for (int i = 0; i < n_Cell; i++) {
    IntegerVector inflow_i = as<IntegerVector>(lst_Inflow[i]);
    inflow_i.push_back(i + 1); // Start with the current cell (1-based index)
    lst_Inflow[i] = inflow_i;

    int num_Ori = i + 1;  // Start with the current cell (1-based index)
    int num_Next = int_Outflow[i]; // Get the next cell from the flow direction vector

    while (num_Next != num_Ori) {
      // Append the current cell to the inflow list of num_Next
      IntegerVector inflow_next = as<IntegerVector>(lst_Inflow[num_Next - 1]);
      inflow_next.push_back(i + 1);  // Adjust for 1-based indexing
      lst_Inflow[num_Next - 1] = inflow_next;

      // Update num_Ori and num_Next for the next iteration
      num_Ori = num_Next;
      num_Next = int_Outflow[num_Next - 1]; // Move to the next cell
    }
  }

  return lst_Inflow;
}

//' @rdname routingtopology
//' @name get_inflow_lastcell
//' @title Get Inflow Last Cell Matrix
//' @description This function creates a matrix of inflow cells for each cell based on the outflow vector.
//' @return A NumericMatrix where each row corresponds to a cell, and each column represents an inflow cell.
//' @export
// [[Rcpp::export]]
NumericMatrix get_inflow_lastcell(IntegerVector int_Outflow) {
  int n_Cell = int_Outflow.size();

  // Create a list to store inflow cells
  std::vector<std::vector<int>> lst_Inflow_LastCell(n_Cell);

  // Populate lst_Inflow_LastCell with indices
  for (int i = 0; i < n_Cell; i++) {
    for (int j = 0; j < n_Cell; j++) {
      if (int_Outflow[j] == i + 1) {  // Adjust for 1-based indexing
        lst_Inflow_LastCell[i].push_back(j + 1);  // Store 1-based index
      }
    }
  }

  // Find the maximum number of inflow cells for any cell
  int max_n_LastCell = 0;
  for (int i = 0; i < n_Cell; i++) {
    if (lst_Inflow_LastCell[i].size() > max_n_LastCell) {
      max_n_LastCell = lst_Inflow_LastCell[i].size();
    }
  }

  // Initialize the matrix with NA values
  NumericMatrix mat_Inflow_LastCell(n_Cell, max_n_LastCell);
  std::fill(mat_Inflow_LastCell.begin(), mat_Inflow_LastCell.end(), NA_REAL);

  // Fill the matrix with the inflow cells
  for (int i = 0; i < n_Cell; i++) {
    for (int j = 0; j < lst_Inflow_LastCell[i].size(); j++) {
      mat_Inflow_LastCell(i, j) = lst_Inflow_LastCell[i][j];
    }
  }

  return mat_Inflow_LastCell;
}


List get_step_cells(List lst_Flow_Cells) {
  int n_Cell = lst_Flow_Cells.size();

  // Step 1: Calculate the lengths of each element in lst_Flow_Cells
  IntegerVector length_Inflow(n_Cell);
  for (int i = 0; i < n_Cell; i++) {
    length_Inflow[i] = Rf_length(lst_Flow_Cells[i]);
  }

  // Step 2: Create a table to count occurrences of each length
  std::unordered_map<int, int> table_Count;
  for (int i = 0; i < n_Cell; i++) {
    table_Count[length_Inflow[i]]++;
  }

  // Step 3: Sort the lengths and assign step numbers
  std::vector<int> sorted_lengths(n_Cell);
  std::copy(length_Inflow.begin(), length_Inflow.end(), sorted_lengths.begin());
  std::sort(sorted_lengths.begin(), sorted_lengths.end());

  IntegerVector i_Step(n_Cell);
  int step = 1;
  for (int i = 0; i < n_Cell; i++) {
    if (i == 0 || sorted_lengths[i] != sorted_lengths[i - 1]) {
      i_Step[i] = step++;
    } else {
      i_Step[i] = i_Step[i - 1];
    }
  }

  // Step 4: Match length_Inflow to sorted_lengths and get corresponding step numbers
  IntegerVector idx_Step_Cell(n_Cell);
  for (int i = 0; i < n_Cell; i++) {
    auto it = std::find(sorted_lengths.begin(), sorted_lengths.end(), length_Inflow[i]);
    idx_Step_Cell[i] = i_Step[std::distance(sorted_lengths.begin(), it)];
  }

  // Step 5: Group cells by their step numbers
  List lst_Step_Cell(step - 1);
  for (int i = 0; i < n_Cell; i++) {
    int idx = idx_Step_Cell[i] - 1;  // Adjust for 0-based indexing
    if (lst_Step_Cell[idx] == R_NilValue) {
      lst_Step_Cell[idx] = IntegerVector::create(i + 1);  // Initialize with 1-based index
    } else {
      IntegerVector temp = as<IntegerVector>(lst_Step_Cell[idx]);
      temp.push_back(i + 1);  // Add 1-based index
      lst_Step_Cell[idx] = temp;
    }
  }

  return lst_Step_Cell;
}


List get_step_lastcell(List lst_Step_Cell, NumericMatrix mat_Inflow_LastCell) {
  int n_Step = lst_Step_Cell.size();
  List lst_Step_InflowLastCell(n_Step);

  // Set the first element to NA (equivalent in Rcpp is R_NilValue)
  lst_Step_InflowLastCell[0] = R_NilValue;

  for (int i = 1; i < n_Step; i++) {
    IntegerVector indices = lst_Step_Cell[i];
    int n_Indices = indices.size();
    NumericMatrix submatrix(n_Indices, mat_Inflow_LastCell.ncol());

    for (int j = 0; j < n_Indices; j++) {
      int idx = indices[j] - 1;  // Adjust for 0-based indexing in C++
      submatrix(j, _) = mat_Inflow_LastCell(idx, _);
    }

    lst_Step_InflowLastCell[i] = submatrix;
  }

  return lst_Step_InflowLastCell;
}


//' @rdname routingtopology
//' @title Get Step Parameters
//' @description This function returns a list of step cells and the corresponding last cell matrices.
//' @return A list containing the step cells and last cell matrices.
//' @export
// [[Rcpp::export]]
List get_routing_info(IntegerVector int_Outflow) {
  // Step 1: Get inflow cells
  List lst_Flow_Cells = get_inflow_cells(int_Outflow);

  // Step 2: Get step cells based on inflow cells
  List lst_Step_Cell = get_step_cells(lst_Flow_Cells);

  // Step 3: Get inflow last cell matrix based on outflow data
  NumericMatrix mat_Inflow_LastCell = get_inflow_lastcell(int_Outflow);

  // Step 4: Get step last cell list based on step cells and inflow last cell matrix
  List lst_Step_LastCell = get_step_lastcell(lst_Step_Cell, mat_Inflow_LastCell);

  // Return the results as a list
  return List::create(
    Named("int_Cell") = lst_Step_Cell,
    Named("mat_LastCell") = lst_Step_LastCell
  );
}

//' @rdname routingtopology
//' @param lst_Inflow_Cell A list of integer vectors, where each vector contains the cells that flow into the respective cell.
//' @param int_OutLet An integer representing the outlet cell (1-based index).
//' @param int_TestCell An integer vector, cells to test.
//' @return An integer vector of cells in the intersection of the station cells and the basin.
//' @export
// [[Rcpp::export]]
IntegerVector get_cell_in_basin(List lst_Inflow_Cell, int int_OutLet, IntegerVector int_TestCell) {
  // Extract the Big Basin (1-based indexing)
  IntegerVector int_BigBasin = lst_Inflow_Cell[int_OutLet - 1];

  // Remove int_OutLet from int_TestCell if present
  IntegerVector int_TestCell_no_outlet = setdiff(int_TestCell, IntegerVector::create(int_OutLet));

  // Convert vectors to sets for efficient operations
  std::set<int> big_basin(int_BigBasin.begin(), int_BigBasin.end());
  std::set<int> station_cells(int_TestCell_no_outlet.begin(), int_TestCell_no_outlet.end());

  // Find the intersection
  std::vector<int> intersection;
  std::set_intersection(
    big_basin.begin(), big_basin.end(),
    station_cells.begin(), station_cells.end(),
    std::back_inserter(intersection)
  );

  // Return the intersection as an IntegerVector
  return wrap(intersection);
}

//' @rdname routingtopology
//' @param int_UpstreamCell An integer vector containing the upstream cells to find the upstream basin.
//' @return An integer vector representing the new upstream basin, which includes the upstream cells and the set difference of the basin cells.
//' This function identifies the upstream basin of a given outlet cell by first finding the intersection of the upstream cells
//' with the cells that flow into the outlet. It then computes the set difference between the upstream basin and the outlet basin.
//' @export
// [[Rcpp::export]]
IntegerVector get_inter_basin(List lst_Inflow_Cell, int int_OutLet, IntegerVector int_UpstreamCell) {
  // Extract the Big Basin
  IntegerVector int_BigBasin = lst_Inflow_Cell[int_OutLet - 1];

  // Collect inflow cells for upstream cells
  std::set<int> upstream_basin;
  for (int cell : int_UpstreamCell) {
    IntegerVector inflow_cells = lst_Inflow_Cell[cell - 1];
    upstream_basin.insert(inflow_cells.begin(), inflow_cells.end());
  }

  // Convert Big Basin to a set
  std::set<int> big_basin(int_BigBasin.begin(), int_BigBasin.end());

  // Compute the set difference: BigBasin - UpstreamBasin
  std::vector<int> remaining_basin;
  std::set_difference(
    big_basin.begin(), big_basin.end(),
    upstream_basin.begin(), upstream_basin.end(),
    std::back_inserter(remaining_basin)
  );

  // Combine upstream cells and remaining basin
  std::vector<int> inter_basin = as<std::vector<int>>(int_UpstreamCell);
  inter_basin.insert(inter_basin.end(), remaining_basin.begin(), remaining_basin.end());

  // Return the result as an IntegerVector
  return wrap(inter_basin);
}




//' @rdname routingtopology
//' @param int_Outflow_Ori An integer vector representing the original outflow indices (1-based).
//' @param int_CellNew An integer vector representing the cells within the new basin.
//' @return An integer vector of the new outflow indices adjusted for the sub-basin.
//' @export
// [[Rcpp::export]]
IntegerVector get_new_outflow(IntegerVector int_Outflow_Ori, IntegerVector int_CellNew) {
  int n_Cell_Ori = int_Outflow_Ori.size();
  int n_interBasin = int_CellNew.size();

  // Create a mapping for the inter-basin
  IntegerVector match_interBasin(n_Cell_Ori, NA_INTEGER);
  for (int i = 0; i < n_interBasin; ++i) {
    if (int_CellNew[i] - 1 < n_Cell_Ori) {
      match_interBasin[int_CellNew[i] - 1] = i + 1; // 1-based index
    }
  }

  // Compute the new outflow indices for the sub-basin
  IntegerVector int_Outflow_Subbasin(n_interBasin, NA_INTEGER);
  for (int i = 0; i < n_interBasin; ++i) {
    int outflow_idx = int_Outflow_Ori[int_CellNew[i] - 1] - 1; // Convert to 0-based index
    if (outflow_idx >= 0 && outflow_idx < n_Cell_Ori) {
      int_Outflow_Subbasin[i] = match_interBasin[outflow_idx];
    }
  }

  // Replace the first NA value with its position (1-based index) and exit the loop
  for (int i = 0; i < n_interBasin; ++i) {
    if (IntegerVector::is_na(int_Outflow_Subbasin[i])) {
      int_Outflow_Subbasin[i] = i + 1; // 1-based index
      break;
    }
  }

  return int_Outflow_Subbasin;
}



//' @rdname routingtopology
//' @param int_CaliCell An integer vector of calibration cells.
//' @return A list of integer vectors (`lst_Step_Cali`), where each element represents calibration cells at a specific step.
//' @export
// [[Rcpp::export]]
List get_cali_step(List lst_Inflow_Cell, IntegerVector int_CaliCell) {
  // Flatten the inflow cells for the calibration cells
  std::vector<int> int_Cell_All;
  for (int i = 0; i < int_CaliCell.size(); ++i) {
    IntegerVector inflow_cells = lst_Inflow_Cell[int_CaliCell[i] - 1]; // 1-based indexing
    int_Cell_All.insert(int_Cell_All.end(), inflow_cells.begin(), inflow_cells.end());
  }

  // Count occurrences of cells within int_CaliCell
  std::unordered_map<int, int> cell_count_map;
  for (int cell : int_Cell_All) {
    if (std::find(int_CaliCell.begin(), int_CaliCell.end(), cell) != int_CaliCell.end()) {
      cell_count_map[cell]++;
    }
  }

  // Get counts for each calibration cell
  std::vector<int> counts(int_CaliCell.size());
  for (int i = 0; i < int_CaliCell.size(); ++i) {
    counts[i] = cell_count_map[int_CaliCell[i]];
  }

  // Find max count
  int max_count = *std::max_element(counts.begin(), counts.end());

  // Directly calculate steps using vectorized arithmetic
  std::vector<int> cell_steps(int_CaliCell.size());
  for (int i = 0; i < counts.size(); ++i) {
    cell_steps[i] = max_count - counts[i] + 1;
  }

  // Group calibration cells by their step
  std::vector<std::vector<int>> step_groups(max_count);
  for (int i = 0; i < int_CaliCell.size(); ++i) {
    step_groups[cell_steps[i] - 1].push_back(int_CaliCell[i]); // Adjust to 0-based index
  }

  // Convert step groups to Rcpp List
  List lst_Step_Cali(max_count);
  for (int i = 0; i < max_count; ++i) {
    lst_Step_Cali[i] = wrap(step_groups[i]);
  }

  return lst_Step_Cali;
}



//' @rdname routingtopology
//' @return A list of integer vectors (`lst_Step_Cali`), where each element represents calibration cells at a specific step.
//' @export
// [[Rcpp::export]]
List get_upstream_cali_cell(List lst_Inflow_Cell, IntegerVector int_CaliCell) {
  int n_CaliCells = int_CaliCell.size();

  // Step 1: Get upstream cells for each calibration cell
  List lst_Cali_Upstream(n_CaliCells);
  for (int i = 0; i < n_CaliCells; ++i) {
    // Get upstream cells for the current calibration cell
    lst_Cali_Upstream[i] = get_cell_in_basin(lst_Inflow_Cell, int_CaliCell[i], int_CaliCell);
  }

  // Step 2: Determine the last calibration cell for each calibration cell
  List lst_LastCaliCell(n_CaliCells);
  for (int i = 0; i < n_CaliCells; ++i) {
    IntegerVector upstream_cells = lst_Cali_Upstream[i];

    // Map upstream cells to their indices in int_CaliCell
    std::unordered_set<int> upstream_indices;
    for (int j = 0; j < upstream_cells.size(); ++j) {
      auto it = std::find(int_CaliCell.begin(), int_CaliCell.end(), upstream_cells[j]);
      if (it != int_CaliCell.end()) {
        upstream_indices.insert(it - int_CaliCell.begin());
      }
    }

    // Collect last calibration cells by removing the common upstream cells
    std::unordered_set<int> temp_set;
    for (int index : upstream_indices) {
      IntegerVector temp = lst_Cali_Upstream[index];
      temp_set.insert(temp.begin(), temp.end());
    }

    std::vector<int> result;
    for (int cell : upstream_cells) {
      if (temp_set.find(cell) == temp_set.end()) {
        result.push_back(cell);
      }
    }

    lst_LastCaliCell[i] = wrap(result);
  }

  return lst_LastCaliCell;
}




//' @rdname routingtopology
//' @export
// [[Rcpp::export]]
void write_int_vector_list(List lst, std::string file_path) {
  std::ofstream fout(file_path, std::ios::binary);
  if (!fout) stop("Cannot open file for writing");
  
  int32_t num_vectors = lst.size();
  fout.write(reinterpret_cast<char*>(&num_vectors), sizeof(int32_t));
  
  for (int i = 0; i < num_vectors; ++i) {
    IntegerVector vec = lst[i];
    int32_t len = vec.size();
    fout.write(reinterpret_cast<char*>(&len), sizeof(int32_t));
    fout.write(reinterpret_cast<char*>(vec.begin()), len * sizeof(int32_t));
  }
  
  fout.close();
}

//' @rdname routingtopology
//' @export
// [[Rcpp::export]]
List read_int_vector_list(std::string file_path) {
  std::ifstream fin(file_path, std::ios::binary);
  if (!fin) stop("Cannot open file for reading");
  
  int32_t num_vectors;
  fin.read(reinterpret_cast<char*>(&num_vectors), sizeof(int32_t));
  if (fin.eof() || num_vectors < 0) stop("Invalid or corrupted file");
  
  List result(num_vectors);
  
  for (int i = 0; i < num_vectors; ++i) {
    int32_t len;
    fin.read(reinterpret_cast<char*>(&len), sizeof(int32_t));
    if (fin.eof() || len < 0) stop("Invalid vector length");
    
    IntegerVector vec(len);
    fin.read(reinterpret_cast<char*>(vec.begin()), len * sizeof(int32_t));
    if (fin.eof()) stop("Unexpected end of file");
    
    result[i] = vec;
  }
  
  fin.close();
  return result;
}


//' @rdname routingtopology
//' @export
// [[Rcpp::export]]
void write_int_matrix_list(List mat_list, std::string file_path) {
  std::ofstream fout(file_path, std::ios::binary);
  if (!fout) stop("Cannot open file for writing");
  
  int32_t num_matrices = mat_list.size();
  fout.write(reinterpret_cast<char*>(&num_matrices), sizeof(int32_t));
  
  for (int i = 0; i < num_matrices; ++i) {
    IntegerMatrix mat = mat_list[i];
    
    int32_t nrow = mat.nrow();
    fout.write(reinterpret_cast<char*>(&nrow), sizeof(int32_t));
    
    fout.write(reinterpret_cast<char*>(mat.begin()), nrow * 9 * sizeof(int32_t));
  }
  
  fout.close();
}

//' @rdname routingtopology
//' @export
// [[Rcpp::export]]
List read_int_matrix_list(std::string file_path) {
  std::ifstream fin(file_path, std::ios::binary);
  if (!fin) stop("Cannot open file for reading");
  
  int32_t num_matrices;
  fin.read(reinterpret_cast<char*>(&num_matrices), sizeof(int32_t));
  if (fin.eof() || num_matrices < 0) stop("Corrupted or empty file");
  
  List result(num_matrices);
  
  for (int i = 0; i < num_matrices; ++i) {
    int32_t nrow;
    fin.read(reinterpret_cast<char*>(&nrow), sizeof(int32_t));
    if (fin.eof() || nrow < 0) stop("Corrupted file: invalid row count");
    
    IntegerMatrix mat(nrow, 9);
    fin.read(reinterpret_cast<char*>(mat.begin()), nrow * 9 * sizeof(int32_t));
    if (fin.eof()) stop("Unexpected end of file");
    
    result[i] = mat;
  }
  
  fin.close();
  return result;
}


