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
  nW.print();
  auto sm{read_submat("data/BLOSUM62.csv")};
  cout << endl;
  // Easy to check: bottom right cell in T should correspond eactly to the score
  // of the optimal alignment, in this case 9.
  cout << "Best alignment score = " << score_align(nW.aln.a, nW.aln.b, sm, -8)
       << endl;
  nW.print();
  seq1 = "PHSWGGAAPHKKRRKSKSPHRW";
  seq2 = "HGWAGGAAGPPHKKKKKSKSPHRW";
  nW.align_sequences(seq1, seq2, -8);
  nW.trace_back();
  nW.aln.print();
  cout << endl;
  cout << "Best alignment score = " << score_align(nW.aln.a, nW.aln.b, sm, -8)
       << endl;

  seq1 = "PHSWGGAAPHKKRRKSKSPHRW";
  seq2 = "HGWAGGAAGPPHKKKKKSKSPHRW";
  nW.align_sequences(seq1, seq2, -1);
  nW.trace_back();
  nW.aln.print();
  cout << endl;
  cout << "Best alignment score = " << score_align(nW.aln.a, nW.aln.b, sm, -1)
       << endl;
}
