#include <RcppArmadillo.h>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

//'  routingtopology
//' @name routingtopology
//' @title Get Inflow Cells
//' @description This function calculates inflow cells based on the outflow vector.
//' @param int_Outflow A vector of integers representing the cell number of the next cell. 1-based indexing.
//' @return A field of uvecs containing the inflow cells for each cell.
//' @export
// [[Rcpp::export]]
arma::field<arma::uvec> get_inflow_cells(const arma::uvec& int_Outflow) {
  int n = int_Outflow.n_elem;
  std::vector<std::vector<arma::uword>> temp(n);

  for (arma::uword i = 0; i < n; ++i) {
    arma::uword origin = i + 1;
    arma::uword next = int_Outflow(i);
    temp[i].push_back(origin);

    while (next != origin) {
      temp[next - 1].push_back(origin);
      origin = next;
      next = int_Outflow(origin - 1);
    }
  }

  arma::field<arma::uvec> result(n);
  for (arma::uword i = 0; i < n; ++i)
    result(i) = arma::uvec(temp[i]);

  return result;
}

//' @rdname routingtopology
//' @title Get Inflow Last Cell Matrix
//' @description Creates a matrix of inflow cells for each cell.
//' @export
// [[Rcpp::export]]
arma::umat get_inflow_lastcell(const arma::uvec& int_Outflow) {
  const arma::uword n = int_Outflow.n_elem;
  std::vector<std::vector<arma::uword>> lst_Inflow_LastCell(n);
  arma::uword max_size = 0;

  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = 0; j < n; ++j) {
      if (int_Outflow(j) == i + 1) {
        lst_Inflow_LastCell[i].push_back(j + 1);
      }
    }
    if (lst_Inflow_LastCell[i].size() > max_size)
      max_size = lst_Inflow_LastCell[i].size();
  }

  arma::umat mat_Inflow_LastCell(n, max_size, arma::fill::zeros);
  for (arma::uword i = 0; i < n; ++i)
    for (arma::uword j = 0; j < lst_Inflow_LastCell[i].size(); ++j)
      mat_Inflow_LastCell(i, j) = lst_Inflow_LastCell[i][j];

  return mat_Inflow_LastCell;
}

//' @rdname routingtopology
//' @description Creates a list of inflow cells for each step.
//' @export
// [[Rcpp::export]]
arma::field<arma::uvec> get_step_cells(const arma::field<arma::uvec>& inflow_cells) {
  int n = inflow_cells.n_elem;
  arma::uvec lengths(n);
  for (int i = 0; i < n; ++i)
    lengths(i) = inflow_cells(i).n_elem;

  std::set<arma::uword> unique_lengths(lengths.begin(), lengths.end());
  std::vector<arma::uword> sorted_lengths(unique_lengths.begin(), unique_lengths.end());

  std::unordered_map<arma::uword, arma::uword> length_to_step;
  for (std::size_t i = 0; i < sorted_lengths.size(); ++i)
    length_to_step[sorted_lengths[i]] = i + 1;

  std::vector<std::vector<arma::uword>> step_groups(sorted_lengths.size());
  for (int i = 0; i < n; ++i)
    step_groups[length_to_step[lengths(i)] - 1].push_back(i + 1);

  arma::field<arma::uvec> result(sorted_lengths.size());
  for (std::size_t i = 0; i < step_groups.size(); ++i)
    result(i) = arma::uvec(step_groups[i]);

  return result;
}


//' @rdname routingtopology
//' @description Creates a field of last cells for each step.
//' @export
// [[Rcpp::export]]
arma::field<arma::umat> get_step_lastcell(const arma::field<arma::uvec>& step_cells,
                                         const arma::umat& inflow_lastcell) {
  arma::field<arma::umat> result(step_cells.n_elem);
  result(0).reset();  // First step is empty

  for (arma::uword i = 1; i < step_cells.n_elem; ++i) {
    arma::uvec idx = step_cells(i) - 1;
    result(i) = inflow_lastcell.rows(idx);
  }

  return result;
}



//' @rdname routingtopology
//' @param int_Outflow An integer vector of outflow cell indices.
//' @param filepath_step_cells Path to save the step_cells field.
//' @param filepath_step_lastcell Path to save the step_lastcell field.
//'
//' @return NULL (invisible)
//' @export
// [[Rcpp::export]]
void generate_step_cell(const arma::uvec& int_Outflow,
                                     const std::string& filepath_step_cells,
                                     const std::string& filepath_step_lastcell) {
  // Generate topology components
  arma::field<arma::uvec> inflow_cells = get_inflow_cells(int_Outflow);
  arma::umat inflow_lastcell = get_inflow_lastcell(int_Outflow);
  arma::field<arma::uvec> step_cells = get_step_cells(inflow_cells);
  arma::field<arma::umat> step_lastcell = get_step_lastcell(step_cells, inflow_lastcell);

  // Save to files
  step_cells.save(filepath_step_cells, arma::arma_binary);
  step_lastcell.save(filepath_step_lastcell, arma::arma_binary);
}



//' @rdname routingtopology
//' @param lst_Inflow_Cell A list of integer vectors, where each vector contains the cells that flow into the respective cell.
//' @param int_OutLet An integer representing the outlet cell (1-based index).
//' @param int_TestCell An integer vector, cells to test.
//' @return An integer vector of cells in the intersection of the station cells and the basin.
//' @export
// [[Rcpp::export]]
arma::uvec get_cell_in_basin(const arma::field<arma::uvec>& lst_Inflow_Cell,
                              int int_OutLet, const arma::uvec& int_TestCell) {
  arma::uvec int_BigBasin = lst_Inflow_Cell(int_OutLet - 1);
  std::set<arma::uword> big_basin_set(int_BigBasin.begin(), int_BigBasin.end());

  std::vector<arma::uword> int_TestCell_no_outlet;
  for (arma::uword val : int_TestCell) {
    if (val != static_cast<arma::uword>(int_OutLet)) {
      int_TestCell_no_outlet.push_back(val);
    }
  }

  std::set<arma::uword> station_cells(int_TestCell_no_outlet.begin(), int_TestCell_no_outlet.end());

  std::vector<arma::uword> intersection;
  std::set_intersection(
    big_basin_set.begin(), big_basin_set.end(),
    station_cells.begin(), station_cells.end(),
    std::back_inserter(intersection)
  );

  return arma::uvec(intersection);
}

//' @rdname routingtopology
//' @param int_UpstreamCell An integer vector containing the upstream cells to find the upstream basin.
//' @return An integer vector representing the new upstream basin, which includes the upstream cells and the set difference of the basin cells.
//' This function identifies the upstream basin of a given outlet cell by first finding the intersection of the upstream cells
//' with the cells that flow into the outlet. It then computes the set difference between the upstream basin and the outlet basin.
//' @export
// [[Rcpp::export]]
arma::uvec get_inter_basin(const arma::uvec& int_Cell, const arma::uvec& int_Outflow) {
  int n = int_Cell.n_elem;
  arma::ivec out(n, arma::fill::value(NA_INTEGER));

  for (int i = 0; i < n; ++i) {
    int id = int_Cell[i] - 1;
    int next = int_Outflow[id];
    out[i] = (std::find(int_Cell.begin(), int_Cell.end(), next) == int_Cell.end()) ? next : NA_INTEGER;
  }

  return arma::conv_to<arma::uvec>::from(out);
}

//' @rdname routingtopology
//' @param int_Outflow_Ori An integer vector representing the original outflow indices (1-based).
//' @param int_CellNew An integer vector representing the cells within the new basin.
//' @return An integer vector of the new outflow indices adjusted for the sub-basin.
//' @export
// [[Rcpp::export]]
arma::uvec get_new_outflow(const arma::uvec& int_Cell, const arma::uvec& int_Outflow) {
  int n = int_Cell.n_elem;
  arma::ivec result(n, arma::fill::value(NA_INTEGER));
  std::unordered_map<arma::uword, arma::uword> old_to_new;

  for (arma::uword i = 0; i < n; ++i)
    old_to_new[int_Cell[i]] = i + 1;

  for (arma::uword i = 0; i < n; ++i) {
    arma::uword id = int_Cell[i];
    arma::uword next = int_Outflow[id - 1];
    result[i] = old_to_new.count(next) ? old_to_new[next] : id;
  }

  return arma::conv_to<arma::uvec>::from(result);
}

//' @rdname routingtopology
//' @param int_CaliCell An integer vector of calibration cells.
//' @return A list of integer vectors (`lst_Step_Cali`), where each element represents calibration cells at a specific step.
//' @export
// [[Rcpp::export]]
arma::uvec get_cali_step(const arma::field<arma::uvec>& step_cells,
                         const arma::uvec& int_Cali) {
  std::vector<arma::uword> steps;
  for (arma::uword i = 0; i < step_cells.n_elem; ++i) {
    for (arma::uword id : step_cells(i)) {
      if (arma::any(int_Cali == id)) {
        steps.push_back(i + 1);
      }
    }
  }

  return arma::uvec(steps);
}

//' @rdname routingtopology
//' @return A list of integer vectors (`lst_Step_Cali`), where each element represents calibration cells at a specific step.
//' @export
// [[Rcpp::export]]
arma::field<arma::uvec> get_upstream_cali_cell(const arma::field<arma::uvec>& lst_Inflow_Cell,
                                               const arma::uvec& int_CaliCell) {
  int n_CaliCells = int_CaliCell.n_elem;

  arma::field<arma::uvec> lst_Cali_Upstream(n_CaliCells);
  for (int i = 0; i < n_CaliCells; ++i) {
    lst_Cali_Upstream(i) = get_cell_in_basin(lst_Inflow_Cell, int_CaliCell(i), int_CaliCell);
  }

  arma::field<arma::uvec> lst_LastCaliCell(n_CaliCells);
  for (int i = 0; i < n_CaliCells; ++i) {
    arma::uvec upstream_cells = lst_Cali_Upstream(i);
    std::unordered_set<arma::uword> upstream_indices;

    for (arma::uword cell : upstream_cells) {
      auto it = std::find(int_CaliCell.begin(), int_CaliCell.end(), cell);
      if (it != int_CaliCell.end()) {
        upstream_indices.insert(it - int_CaliCell.begin());
      }
    }

    std::unordered_set<arma::uword> temp_set;
    for (arma::uword index : upstream_indices) {
      arma::uvec temp = lst_Cali_Upstream(index);
      temp_set.insert(temp.begin(), temp.end());
    }

    std::vector<arma::uword> result;
    for (arma::uword cell : upstream_cells) {
      if (temp_set.find(cell) == temp_set.end()) {
        result.push_back(cell);
      }
    }

    lst_LastCaliCell(i) = arma::uvec(result);
  }

  return lst_LastCaliCell;
}
