#include "08_MultipleSequenceAlignment.h"
#include "06_SequenceAlignment.h"
#include "07_DatabaseSearch.h"

int MSA::multiple_alignment::get_score() {
  int score = 0;
  if (alignments.empty())
    return 0;
  for (unsigned int i = 0; i < alignments[0].second.size(); i++) {
    for (auto &seq : alignments) {

      seq.second[i];
    }
  }
};
std::string MSA::multiple_alignment::get_consensus();
void MSA::multiple_alignment::print();

MSA::msa::msa(const int match, const int mismatch, const std::string &alphabet)
    : nw{match, mismatch, alphabet} {
  sm = create_submat(match, mismatch, alphabet);
};
MSA::msa::msa(const std::string &file) : nw{file} { sm = read_submat(file); };
