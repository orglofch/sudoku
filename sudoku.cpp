/*
* Copyright (c) 2015 Owen Glofcheski
*
* This software is provided 'as-is', without any express or implied
* warranty. In no event will the authors be held liable for any damages
* arising from the use of this software.
*
* Permission is granted to anyone to use this software for any purpose,
* including commercial applications, and to alter it and redistribute it
* freely, subject to the following restrictions:
*
*    1. The origin of this software must not be misrepresented; you must not
*    claim that you wrote the original software. If you use this software
*    in a product, an acknowledgment in the product documentation would be
*    appreciated but is not required.
*
*    2. Altered source versions must be plainly marked as such, and must not
*    be misrepresented as being the original software.
*
*    3. This notice may not be removed or altered from any source
*    distribution.
*/

#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <stack>
#include <string>
#include <unordered_set>
#include <vector>

#define FAILURE_THRESHOLD 10000

using namespace std;

// Types
struct Cell {
  Cell() : r(-1), c(-1) {
  }

  Cell(int r, int c) : r(r), c(c) {
  }

  int r, c;
};

typedef int Puzzle[9][9];
typedef unordered_set<int> PossibilitySet;
typedef PossibilitySet PossibilityMap[9][9];
typedef vector<Cell> CellList;

// Output operators
ostream &operator<< (ostream &os, const Cell &cell) {
  os << "[" << cell.r << ", " << cell.c << "]";
  return os;
}

ostream &operator<< (ostream &os, const Puzzle &puzzle) {
  for (int r = 0; r < 9; ++r) {
    for (int c = 0; c < 9; ++c) {
      os << puzzle[r][c] << " ";
    }
    if (r != 8) {
      os << endl;
    }
  }
  return os;
}

ostream &operator<< (ostream &os, const PossibilitySet &possibility_set) {
  os << "{";
  for (PossibilitySet::const_iterator it = possibility_set.begin();
      it != possibility_set.end(); ++it) {
    if (it != possibility_set.begin()) {
      os << ", ";
    }
    os << *it;
  }
  os << "}";
  return os;
}

ostream &operator<< (ostream &os, const CellList &cell_list) {
  os << "{";
  for (CellList::const_iterator it = cell_list.begin(); 
      it != cell_list.end(); ++it) {
    if (it != cell_list.begin()) {
      os << ", ";
    }
    os << *it;
  }
  os << "}";
  return os;
}

/*
 * Print usage to stdout
 */
void Usage() {
  cout << "Usage: sudoku [data-file]" << endl;
}

/*
 * Load a Sudoku puzzle in the file <filename> into a 9x9 array
 * Returns true if the load succeeded 
 *         false otherwise
 */
bool LoadPuzzle(const string &filename, Puzzle &puzzle) {
  ifstream ifs(filename);
  if (!ifs.is_open()) {
    return false;
  }

  for (int r = 0; r < 9; ++r) {
    for (int c = 0; c < 9; ++c) {
      ifs >> puzzle[r][c];
    }
  }
  return true;
}

/*
 * Returns true if the given row is legal in the puzzle
 *         false otherwise
 */
bool IsLegalRow(const Puzzle &puzzle, int r) {
  assert(r < 9);
  bool checker[9] = {false};
  for (int c = 0; c < 9; ++c) {
    if (puzzle[r][c] != 0 && checker[puzzle[r][c] - 1] == true) {
      cout << "Row " << r << " is not legal" << endl;
      return false;
    } else {
      checker[puzzle[r][c] - 1] = true;
    }
  }
  return true;
}

/*
 * Returns true if the given column is legal in the puzzle
 *         false otherwise
 */
bool IsLegalColumn(const Puzzle &puzzle, int c) {
  assert(c < 9);
  bool checker[9] = {false};
  for (int r = 0; r < 9; ++r) {
    if (puzzle[r][c] != 0 && checker[puzzle[r][c] - 1] == true) {
      cout << "Column " << c << " is not legal" << endl;
      return false;
    } else {
      checker[puzzle[r][c] - 1] = true;
    }
  }
  return true;
}

/*
 * Returns true of the given submatrix is legal in the puzzle
 *         false otherwise
 */
bool IsLegalSubmatrix(const Puzzle &puzzle, int sub_r, int sub_c) {
  assert(sub_r <= 2 && sub_c <= 2);
  bool checker[9] = {false};
  for (int r = sub_r * 3; r < sub_r * 3 + 3; ++r) {
    for (int c = sub_c * 3; c < sub_c * 3 + 3; ++c) {
      if (puzzle[r][c] != 0 && checker[puzzle[r][c] - 1] == true) {
        cout << "Submatrix " << sub_r << " " << sub_c << " is not legal" << endl;
        return false;
      } else {
        checker[puzzle[r][c] - 1] = true;
      }
    }
  }
  return true;
}

/*
 * Returns true if a Sudoku Puzzle is legal
 *         false otherwise
 */
bool IsLegal(const Puzzle &puzzle) {
  // Check rows and columns
  for (int i = 0; i < 9; ++i) {
    if (!IsLegalRow(puzzle, i) || !IsLegalColumn(puzzle, i)) {
      return false;
    }
  }

  // Check submatrices
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (!IsLegalSubmatrix(puzzle, i, j)) {
        return false;
      }
    }
  }
  return true;
}

/*
 * Adds a constraint to the possibility_map based on an action.
 * Updates the action to store the list of cells it constrained
 * Returns true if the possibility map passes forward checking
 *         false otherwise
 */
bool AddConstraints(PossibilityMap &possibility_map, 
                    const int row,
                    const int col,
                    const int val,
                    CellList &constrained_cells) {
  // Add constraints to rows
  for (int c = 0; c < 9; ++c) {
    if (c != col) {
      PossibilitySet &possibility_set = possibility_map[row][c];
      PossibilitySet::iterator it = possibility_set.find(val);
      if (it != possibility_set.end()) {
        possibility_set.erase(it);
        constrained_cells.push_back(Cell(row, c));
        if (possibility_set.size() == 0) {
          return false;
        }
      }
    }
  }
  // Add constraints to columns
  for (int r = 0; r < 9; ++r) {
    if (r != row) {
      PossibilitySet &possibility_set = possibility_map[r][col];
      PossibilitySet::iterator it = possibility_set.find(val);
      if (it != possibility_set.end()) {
        possibility_set.erase(it);
        constrained_cells.push_back(Cell(r, col));
        if (possibility_set.size() == 0) {
          return false;
        }
      }
    }
  }
  // Add constraints to submatrices
  int sub_r = (int)(row / 3);
  int sub_c = (int)(col / 3);
  for (int r = sub_r * 3; r < sub_r * 3 + 3; ++r) {
    for (int c = sub_c * 3; c < sub_c * 3 + 3; ++c) {
      if (r != row && c != col) {
        PossibilitySet &possibility_set = possibility_map[r][c];
        PossibilitySet::iterator it = possibility_set.find(val);
        if (it != possibility_set.end()) {
          possibility_set.erase(it);
          constrained_cells.push_back(Cell(r, c));
          if (possibility_set.size() == 0) {
            return false;
          }
        }
      }
    }
  }
  return true;
}

/*
 * Removes a constraint from the possibility_map essentially
 * undoing an action
 */
void RemoveConstraints(PossibilityMap &possibility_map, 
                       const int row,
                       const int col,
                       const int val,
                       const CellList &constrained_cells) {
  for (CellList::const_iterator it = constrained_cells.begin();
      it != constrained_cells.end(); ++it) {
    assert(it->r != row || it->c != col);
    possibility_map[it->r][it->c].insert(val);
  }
}

/*
 * Returns how many contraints are placed on other variables by
 * a given assignment. Used to find least constraining value.
 */
int NumConstraintsAdded(const PossibilityMap &possibility_map, 
                        const int row,
                        const int col,
                        const int val) {
  // Check constraints on rows
  int constraints = 0;
  for (int c = 0; c < 9; ++c) {
    const PossibilitySet &possibility_set = possibility_map[row][c];
    if (c != col && possibility_set.find(val) != possibility_set.end()) {
      ++constraints;
    }
  }
  // Check constraints on columns
  for (int r = 0; r < 9; ++r) {
    const PossibilitySet &possibility_set = possibility_map[r][col];
    if (r != row && possibility_set.find(val) != possibility_set.end()) {
      ++constraints;
    }
  }
  // Check constraints on submatrices
  int sub_r = (int)(row / 3);
  int sub_c = (int)(col / 3);
  for (int r = sub_r * 3; r < sub_r * 3 + 3; ++r) {
    for (int c = sub_c * 3; c < sub_c * 3 + 3; ++c) {
      const PossibilitySet &possibility_set = possibility_map[r][c];
      // Don't double count elements of the submatrix
      if (r != row && c != col && 
          possibility_set.find(val) != possibility_set.end()) {
        ++constraints;
      }
    }
  }
  return constraints;
}

/*
 * Returns a list of the most constrained variables for a puzzle
 */
void MostConstrained(const Puzzle &puzzle, 
                     const PossibilityMap &possibility_map, 
                     const int (&row_variables)[9],
                     const int (&column_variables)[9],
                     CellList &cells) {
  int max_constraints = 0;
  // For each variable check its number of possibilities
  for (int r = 0; r < 9; ++r) {
    for (int c = 0; c < 9; ++c) {
      // Only count cells which are open
      if (puzzle[r][c] == 0) {
        assert(possibility_map[r][c].size() > 0);
        int constraints = 9 - possibility_map[r][c].size();
        if (constraints > max_constraints) {
          // We have a new most constraining
          max_constraints = constraints;
          cells.clear();
          cells.push_back(Cell(r, c));
        } else if (constraints == max_constraints) {
          // We have a tie for most constraining
          cells.push_back(Cell(r, c));  
        }
      }
    }
  }
}

/*
 * Returns the most contraining cell in a list of cells
 */
Cell MostConstraining(const Puzzle &puzzle,
                      const int (&row_variables)[9],
                      const int (&column_variables)[9],
                      const CellList &cells) {
  Cell most_constraining(-1, -1);
  int max_constraints = -1;
  for (CellList::const_iterator it = cells.begin(); it != cells.end(); ++it) {
    // Count submatrix variables seperately since they could double count
    // variables otherwise
    int submatrix_variables = 0;
    int sub_r = (int)(it->r / 3);
    int sub_c = (int)(it->c / 3);
    for (int r = sub_r * 3; r < sub_r * 3 + 3; ++r) {
      for (int c = sub_c * 3; c < sub_c * 3 + 3; ++c) {
        // Only count open variables and don't variables along row or column
        if (r != it->r && c != it->c && puzzle[r][c] == 0) {
          ++submatrix_variables;
        }
      }
    }
    // Subtract 1 since it shouldn't count itself
    // even though this wouldn't affect correctness since ordering
    // would be preserved
    int constraints = (row_variables[it->r] - 1) + 
                      (column_variables[it->c] - 1) +
                      submatrix_variables;
    if (constraints > max_constraints) {
      max_constraints = constraints;
      most_constraining = *it; 
    }
  }
  return most_constraining;
}

/*
 * Return the least constraining value for a cell given it's possible values
 */
int LeastConstrainingValue(const PossibilityMap &possibility_map,
                           const Cell &selected_cell) {
  int least_constraining = -1;
  int min_constraints_added = numeric_limits<int>::max();
  assert(possibility_map[selected_cell.r][selected_cell.c].size() > 0);
  // For each possibility check which creates the least constraints
  const PossibilitySet &possibility_set = 
      possibility_map[selected_cell.r][selected_cell.c];
  for (PossibilitySet::const_iterator it = possibility_set.begin(); 
      it != possibility_set.end(); ++it) {
    int constraints_added = 
        NumConstraintsAdded(possibility_map, selected_cell.r, selected_cell.c, *it);
    // If this value adds fewer constraints select it
    if (constraints_added < min_constraints_added) {
      min_constraints_added = constraints_added;
      least_constraining = *it;
    }
  }
  return least_constraining;
}

/*
 * Solves a given Sudoku puzzle recursively
 * Call wrapper for easy setup
 */
bool SolvePuzzle(Puzzle &puzzle, 
                 PossibilityMap &possibility_map,
                 int (&row_variables)[9],  
                 int (&column_variables)[9],
                 int open_count,
                 int &variable_assignments,
                 int depth,
                 bool print_depth,
                 bool print_intermediate) {
  // If there are no variables left to assign
  if (open_count == 0) {
    return true;
  }

  if (print_depth) {
    cout << depth << endl;
  }
  if (print_intermediate) {
    cout << puzzle << endl << endl;
  }

  // Select most constrained variables
  CellList most_constrained_list;
  MostConstrained(puzzle, possibility_map, row_variables, 
      column_variables, most_constrained_list);
  assert(most_constrained_list.size() > 0);
  
  // Select most constraining variable as a tie-breaker
  Cell selected_cell = MostConstraining(puzzle, row_variables, 
      column_variables, most_constrained_list); 
  assert(selected_cell.r != -1 && selected_cell.c != -1);

  // Decrease number of open row and column variables
  --row_variables[selected_cell.r];
  --column_variables[selected_cell.c];

  // Select least constraining value for the selected variable
  PossibilitySet &possibility_set = possibility_map[selected_cell.r][selected_cell.c];
  PossibilitySet values_tried;
  
  unsigned int value_count = possibility_set.size();
  for (unsigned int i = 0; i < value_count; ++i) {
    // Early exit if we expand beyond the threshold
    if (variable_assignments == FAILURE_THRESHOLD) {
      return false;
    }
    ++variable_assignments;
    
    // Select least constraining value for the chosen cell
    int least_constraining = LeastConstrainingValue(possibility_map, selected_cell);
    assert(least_constraining > 0);
    
    // Update puzzle
    puzzle[selected_cell.r][selected_cell.c] = least_constraining;

    // Add constraints and perform forward checking
    // NOTE: I assume we don't care if I apply forwarding checking here or
    // in the next step
    CellList constrained_cells;
    bool passes_forward_checking = AddConstraints(possibility_map, 
        selected_cell.r, selected_cell.c, least_constraining, constrained_cells);
    
    if (passes_forward_checking && SolvePuzzle(puzzle, 
                                               possibility_map, 
                                               row_variables, 
                                               column_variables, 
                                               open_count - 1, 
                                               variable_assignments,
                                               depth + 1,
                                               print_depth,
                                               print_intermediate)) {
      return true;
    } else { // Backtrack
      // Undo puzzle constraints
      puzzle[selected_cell.r][selected_cell.c] = 0;
      RemoveConstraints(possibility_map, selected_cell.r, selected_cell.c,
          least_constraining, constrained_cells);
    
      // Update available values
      possibility_set.erase(least_constraining);
      values_tried.insert(least_constraining);
    }
  } 

  // Undo decrease to number of open variables
  ++row_variables[selected_cell.r];
  ++column_variables[selected_cell.c];

  // Set the available variables back
  swap(possibility_set, values_tried);

  return false;
}

/*
 * Wrapper function for solving a given Sudoku puzzle
 * handles initialization of puzzle specific variables
 */
void SolvePuzzle(Puzzle &puzzle, 
                 bool print_assignments,
                 bool print_solution,
                 bool print_depth,
                 bool print_intermediate) {
  // Used to keep track of available values for each cell
  PossibilityMap possibility_map;

  // Count of variables to use when deciding most constraining 
  int row_variables[9] = {0};
  int column_variables[9] = {0};
  int open_count = 0;

  // Populate possibility map
  for (int r = 0; r < 9; ++r) {
    for (int c = 0; c < 9; ++c) {
      if (puzzle[r][c] != 0) {
        possibility_map[r][c].insert(puzzle[r][c]);
      } else {
        for (int i = 1; i <= 9; ++i) {
          possibility_map[r][c].insert(i);
        }
      }
    }
  }

  // Apply inital constraints to possibility map
  // and populate variable counters
  for (int r = 0; r < 9; ++r) {
    for (int c = 0; c < 9; ++c) {
      if (puzzle[r][c] != 0) {
        CellList constrained_cells;
        if (!AddConstraints(possibility_map, r, c, puzzle[r][c], constrained_cells)) {
          if (print_assignments) {
            cout << 0 << endl;
          }
          cout << "Puzzle could not be solved with initial constraints" << endl;
          return;
        }
      } else {
        ++open_count;
        ++row_variables[r];
        ++column_variables[c];
      }
    }
  }

  int variable_assignments = 0;

  // Try to solve the puzzle
  if (SolvePuzzle(puzzle, 
                  possibility_map, 
                  row_variables, 
                  column_variables, 
                  open_count, 
                  variable_assignments, 
                  1,
                  print_depth,
                  print_intermediate)) {
    if (print_assignments) {
      cout << variable_assignments << endl;
    }
    // Ensure puzzle is legal
    assert(IsLegal(puzzle));
    if (print_solution) {
      cout << puzzle << endl;
    }
  } else {
    if (print_assignments) {
      cout << variable_assignments << endl;
    }
    cout << "Puzzle could not be solved" << endl;
  }
}

int main(int argc, char **argv) {
  if (argc < 2) {
    Usage();
    return 1;
  }

  // Default arguments 
  string datafile(argv[1]);
  bool print_assignments = false;
  bool print_solution = true;
  bool print_depth = false;
  bool print_intermediate = false;

  // Extract flags
  for (int i = 2; i < argc; ++i) {
    if (strcmp(argv[i], "-pa") == 0) {
      print_assignments = true;
    } else if (strcmp(argv[i], "-ss") == 0) {
      print_solution = false;
    } else if (strcmp(argv[i], "-pd") == 0) {
      print_depth = true;
    } else if (strcmp(argv[i], "-pi") == 0) {
      print_intermediate = true;
    }
  }

  // Load
  Puzzle puzzle;
  if (!LoadPuzzle(datafile, puzzle)) {
    cout << "Failed to load " << datafile << endl;
    return 1;
  }

  if (!IsLegal(puzzle)) {
    cout << "Initial puzzle was not legal" << endl;
    return 1;
  }

  // Solve
  SolvePuzzle(puzzle, print_assignments, print_solution, 
      print_depth, print_intermediate);

  return 0;
}
