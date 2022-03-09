#include "06_SequenceAlignment.h"
#include <iostream>
#include <unordered_map>

using namespace std;

int main() {
  const string seq1{"PHSWG"};
  const string seq2{"HGWAG"};
  needleman_Wunsch nW{"data/BLOSUM62.csv"};
  nW.align_sequences(seq1, seq2, -8);
  nW.trace_back();
  nW.print();
  auto sm{read_submat("data/BLOSUM62.csv")};
  cout << endl;
  cout << "Best alignment score = " << score_align(nW.aln.a, nW.aln.b, sm, -8)
       << endl;
}
