#include "06_SequenceAlignment.h"
#include <iostream>
#include <unordered_map>

/*
 * ex04 - Handle the situation where there are multiple optimal alignments usign
 * the Needleman-Wunsch algorithm.
 */

using namespace std;

int main() {
  string seq1{"PHHSWWG"}; // There are 2 branches
  string seq2{"PHGWAG"};
  needleman_Wunsch nW{"data/BLOSUM62.csv"};
  nW.align_sequences_withties(seq1, seq2, -8);
  nW.trace_back_withties();
  nW.print_withties();
  cout << endl;
  cout << "\nGap cost = -8 "
          "------------------------------------------------------"
       << endl;

  /* Easy to check: bottom right cell in T should correspond eactly to the
   * score of the optimal alignment.
   */

  seq1 = "PHSWGGAAPHKKKKKRRKSKSPHRWHR";
  seq2 = "HGWAGGAAGPPHKKSKSPHSPHR";
  cout << "seq1: " << seq1 << endl;
  cout << "seq2: " << seq2 << endl;
  cout << endl;
  cout << "\nGap cost = -8 "
          "------------------------------------------------------"
       << endl;
  nW.align_sequences_withties(seq1, seq2, -8);
  nW.trace_back_withties();
  nW.print_withties();
}
