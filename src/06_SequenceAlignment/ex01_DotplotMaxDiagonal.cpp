#include "06_SequenceAlignment.h"
#include <array>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <vector>

/* ex3 - Write a test functions that returns the length and the starting
 * position of the longest identity diagonal. It has to be noted that the
 * exercise only seem to ask only for diagonals going into the "forward"
 * direction because the function should return only the starting position and
 * length of the diagonal. This means it doesn't look at reversed alignments (a
 * sequence that is the reverse of the other for example).
 */

using namespace std;

int main() {
  cout << endl;
  const string seq1{"TTTACGTC"};
  const string seq2{"TTTAGCGTC"};
  DotPlot dp{seq1, seq2};
  cout << "Seq 1 = " << seq1 << endl;
  cout << "Seq 2 = " << seq2 << endl;
  dp.print();
  cout << "\nBasic comparison ----------------------------" << endl;
  dp.compare();
  dp.print();

  cout << "\nDiagonals -----------------------------------" << endl;
  vector<diagonal> v = dp.max_diagonal();
  for (diagonal &x : v) {
    cout << "row = " << x.row << "\tcol = " << x.col
         << "\tlength = " << x.length << endl;
  }
}
