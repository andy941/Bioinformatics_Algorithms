#include "06_SequenceAlignment.h"
#include <iostream>
#include <unordered_map>

using namespace std;

int main() {
  const string seq1{"GATCCATGATGA"};
  const string seq2{"AGTCCATGATGGTATGAT"};
  needleman_Wunsch nW{2, -2, "ATCG"};
  nW.print();
  nW.align_sequences(seq1, seq2, -2);
  nW.print();
}
