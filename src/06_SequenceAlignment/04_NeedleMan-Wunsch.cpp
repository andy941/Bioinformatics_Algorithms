#include "06_SequenceAlignment.h"
#include <iostream>
#include <unordered_map>

using namespace std;

int main() {
  string seq1{"PHSWG"};
  string seq2{"HGWAG"};
  needleman_Wunsch nW{"data/BLOSUM62.csv"};
  nW.align_sequences(seq1, seq2, -8);
  nW.trace_back();
  cout << "\nSequences From the Book "
          "-------------------------------------------"
       << endl;
  nW.print();
  auto sm{read_submat("data/BLOSUM62.csv")};
  cout << endl;
  // Easy to check: bottom right cell in T should correspond eactly to the score
  // of the optimal alignment, in this case 9.
  cout << "Best alignment score = " << score_align(nW.aln.a, nW.aln.b, sm, -8)
       << endl;

  cout << "\nGap Test -----"
          "------------------------------------------------------"
       << endl;
  seq1 = "PHSWGGAAPHKKRRKSKSPHRWAAPHKKRRTLLDWSPHRLL";
  seq2 = "HGWAGGAAGPPHKKKKKSKSPHRWAAPHKRTLLDWSKSPHR";
  cout << "seq1: " << seq1 << endl;
  cout << "seq2: " << seq2 << endl;
  cout << endl;
  cout << "\nGap cost = -8 "
          "------------------------------------------------------"
       << endl;
  nW.align_sequences(seq1, seq2, -8);
  nW.trace_back();
  nW.aln.print();
  cout << endl;
  cout << "Best alignment score = " << score_align(nW.aln.a, nW.aln.b, sm, -8)
       << endl;

  cout << "\nGap cost = -120 "
          "----------------------------------------------------"
       << endl;
  nW.align_sequences(seq1, seq2, -120);
  nW.trace_back();
  nW.aln.print();
  cout << endl;
  cout << "Best alignment score = " << score_align(nW.aln.a, nW.aln.b, sm, -120)
       << endl;
}
