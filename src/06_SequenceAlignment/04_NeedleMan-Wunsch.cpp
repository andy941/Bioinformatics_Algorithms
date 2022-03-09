#include "06_SequenceAlignment.h"
#include "Eigen/Core"
#include "Eigen/src/Core/Matrix.h"
#include <iostream>
#include <unordered_map>

using namespace std;

int main() {
  Eigen::Matrix<int, 4, 4> m;
  m(0, 0) = 4444;
  const string seq1{"GGCC"};
  const string seq2{"GGTTTCC"};
  needleman_Wunsch nW{3, -1, "ATCG"};
  nW.print();
  nW.align_sequences(seq1, seq2, 0);
  nW.print();
  cout << m << endl;
}
