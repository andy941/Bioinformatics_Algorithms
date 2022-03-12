#include "06_SequenceAlignment.h"
#include <iostream>
#include <unordered_map>

using namespace std;

int main() {
  string seq1{"PHSWG"};
  string seq2{"HGWAG"};
  smith_Waterman sW{"data/BLOSUM62.csv"};
  sW.align_sequences(seq1, seq2, -8);
  sW.trace_back();
  cout << "\nSequences From the Book "
          "-------------------------------------------"
       << endl;
  sW.print();
  auto sm{read_submat("data/BLOSUM62.csv")};
  cout << endl;
  // Easy to check: bottom right cell in T should correspond eactly to the score
  // of the optimal alignment, in this case 9.
  cout << "Best alignment score = " << score_align(sW.aln.a, sW.aln.b, sm, -8)
       << endl;

  cout << "\nGap Test -----"
          "------------------------------------------------------"
       << endl;
  seq1 = "PFRFFFAASKSPHRWAAPHKKMVCCVMTAAHR";
  seq2 = "HGWAGGAASKSPHRWAAPHKKKRTLLDWSKSPHR";
  cout << "seq1: " << seq1 << endl;
  cout << "seq2: " << seq2 << endl;
  cout << endl;
  cout << "\nGap cost = -8 "
          "------------------------------------------------------"
       << endl;
  sW.align_sequences(seq1, seq2, -4);
  sW.trace_back();
  sW.aln.print();
  cout << endl;
  cout << "Best alignment score = " << score_align(sW.aln.a, sW.aln.b, sm, -4)
       << endl;

  cout << "\nGap cost = -20 "
          "----------------------------------------------------"
       << endl;
  sW.align_sequences(seq1, seq2, -20);
  sW.trace_back();
  sW.aln.print();
  cout << endl;
  cout << "Best alignment score = " << score_align(sW.aln.a, sW.aln.b, sm, -20)
       << endl;
}
