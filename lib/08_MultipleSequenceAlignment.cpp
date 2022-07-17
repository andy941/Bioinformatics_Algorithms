#include "08_MultipleSequenceAlignment.h"
#include "06_SequenceAlignment.h"
#include "07_DatabaseSearch.h"

int MSA::sum_of_pairs(const std::string &str,
                      std::unordered_map<std::string, int> &sm,
                      const int gap_cost) {
  int score = 0;
  for (unsigned int i = 0; i < str.size() - 1; i++) {
    for (unsigned int j = i + 1; j < str.size(); j++) {
      score += score_pos(str[i], str[j], sm, gap_cost);
    }
  }
  return score;
};

int MSA::multiple_alignment::get_score(std::unordered_map<std::string, int> &sm,
                                       const int gap_cost) {
  int score = 0;
  if (aln_pos.empty())
    return 0;
  for (auto &x : aln_pos) {
    score += sum_of_pairs(x, sm, gap_cost);
  }
  return score;
};

std::string MSA::multiple_alignment::get_consensus() {
  std::string cons_aln{};
  std::vector<std::pair<std::string, unsigned int>> letter_count{};
  if (aln_pos.empty())
    return cons_aln;
  cons_aln.reserve(aln_pos.size());
  for (auto &x : aln_pos) {
    find()
  }
  return cons_aln;
};

void MSA::multiple_alignment::print(){};

MSA::msa::msa(const int match, const int mismatch, const std::string &alphabet)
    : nw{match, mismatch, alphabet} {
  sm = create_submat(match, mismatch, alphabet);
};
MSA::msa::msa(const std::string &file) : nw{file} { sm = read_submat(file); };
