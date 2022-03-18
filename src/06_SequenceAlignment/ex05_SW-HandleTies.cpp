#include "06_SequenceAlignment.h"
#include <iostream>
#include <unordered_map>

/*
 * ex04 - Handle the situation where there are multiple optimal alignments usign
 * the Smith-Waterman algorithm. The improvement that has to be done here
 * concerns the multiple "starting points" where there are mutliple best scores
 * in the S matrix. I chose the two sequences that result both in multiple
 * starting points and tied alignment scores.
 */

using namespace std;

int main() {
  string seq1{"PHSRRWWWPHHHGWWAG"};
  string seq2{"PHHHGWWWWG"};
  // string seq1{"PHHSWWG"};
  // string seq2{"PHGWAG"};
  smith_Waterman sW{"data/BLOSUM62.csv"};
  sW.align_sequences_withties(seq1, seq2, -8);
  sW.trace_back_withties();
  sW.print_withties();
  cout << endl;

  // seq1 = "PHSWGGAAPHKKKKKRRKSKSPHRWHR";
  // seq2 = "HGWAGGAAGPPHKKSKSPHSPHR";
  // cout << "seq1: " << seq1 << endl;
  // cout << "seq2: " << seq2 << endl;
  // cout << endl;
  // cout << "\nGap cost = -8 "
  //        "------------------------------------------------------"
  //     << endl;
  // sW.align_sequences(seq1, seq2, -8);
  // sW.trace_back();
  // sW.print();
}
