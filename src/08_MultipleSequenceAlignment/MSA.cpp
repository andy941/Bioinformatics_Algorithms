#include "08_MultipleSequenceAlignment.h"
#include <string>
#include <utility>
#include <vector>

int main() {
  /*
   * Should give back:
   * AT-AGC
   * A--ACC
   * ATGAC-
   */

  std::vector<std::pair<std::string, std::string>> sequences{
      std::make_pair("seq1", "ATAGC"), // std::make_pair("seq2", "AACC"),
      std::make_pair("seq3", "ATGAC")};
  msa::msa msa{1, -1, "ATCG"};
  msa.align_sequences(sequences, -1);
  msa.print(3);
}
