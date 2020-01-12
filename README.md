Sudoku-Solver
============

Sudoku Solver using entropy heuristic and back-tracking search.

Language: C++

---

The algorithm starts by applying the row, column and submatrix contraints
to generate a list of all possible values for each unknown cell.

It then performs the following:

1. Pick the **most constrained cell** (cell with fewest possibilities) so we
   enouncter an unsolvable situation as quickly as possible.
   
2. Pick the **least contraining value** for cell (the value which applies the
   fewest contraints to all other cells) so we try leave open as many
   solutions in this solution branch as possible.
   
3. Apply contraints created by the new cell on all other cells.
   
   1. Check if the puzzle is solved. If it is return the solution. 
   
   2. Perform **forward checking** to see if the puzzle is now unsolvable
      (some cell has no values) then **propagate the unsolvable value
      backwards** as long as it's still unsolvable. If we hit the top
      of the return stack then exit, otherwise try new values when
      available.
