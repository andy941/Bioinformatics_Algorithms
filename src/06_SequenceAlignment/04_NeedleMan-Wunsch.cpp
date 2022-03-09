#include "06_SequenceAlignment.h"
#include <iostream>
#include <unordered_map>

using namespace std;

int main() {
  const string seq1{"GGCC"};
  const string seq2{"GGTTTCC"};
  needleman_Wunsch nW{3, -1, "ATCG"};
  nW.print();
  nW.align_sequences(seq1, seq2, 0);
  nW.print();
}
