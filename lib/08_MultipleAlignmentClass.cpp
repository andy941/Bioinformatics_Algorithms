#include "08_MultipleAlignmentClass.h"
#include "06_SequenceAlignment.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <ostream>

int msa::sum_of_pairs(const std::string &str,
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

void msa::multiple_alignment::add_sequence(const std::string &name,
                                           const std::string &seq) {
  if (aln_pos.size() == 0) {
    seq_names.push_back(name);
    for (size_t i = 0; i < seq.size(); i++) {
      aln_pos.push_back(std::string{seq[i]});
    }
    return;
  }
  if (seq.size() != aln_pos.size()) {
    std::cerr << "trying to insert a sequence of size " << seq.size()
              << " in alignment of size " << aln_pos.size() << std::endl;
    return;
  }
  for (size_t i = 0; i < seq.size(); i++) {
    aln_pos[i].push_back(seq[i]);
  }
  seq_names.push_back(name);
};

int msa::multiple_alignment::get_score(std::unordered_map<std::string, int> &sm,
                                       const int gap_cost) {
  int score = 0;
  if (aln_pos.empty())
    return 0;
  for (auto &x : aln_pos) {
    score += sum_of_pairs(x, sm, gap_cost);
  }
  return score;
};

std::string msa::multiple_alignment::get_consensus() {
  std::string cons_aln{};
  if (aln_pos.empty())
    return cons_aln;
  for (auto &x : aln_pos) {
    std::vector<char> s_x;
    s_x.reserve(x.size());
    std::partial_sort_copy(x.begin(), x.end(), s_x.begin(), s_x.end());
    unsigned int count{0}, curr_count{0};
    char c = {*x.begin()};
    for (char &curr_c : x) {
      if (curr_c == '-')
        continue;
      if (curr_c != c) {
        curr_count = 0;
      }
      curr_count++;
      if (curr_count > count) {
        count = curr_count;
        c = curr_c;
      }
    }
    cons_aln.push_back(c);
  }
  return cons_aln;
};

void msa::multiple_alignment::print(size_t line_size) {
  for (size_t i = 0; i < aln_pos.size(); i += line_size) {
    for (size_t k = 0; k < seq_names.size(); k++) {
      if (seq_names[k].size() > 10)
        std::cout << seq_names[k].substr(0, 10) << '\t';
      else
        std::cout << std::setw(10) << seq_names[k] << '\t';

      for (size_t j = i; j < std::min(i + line_size, aln_pos.size()); j++) {
        std::cout << aln_pos[j][k];
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
};
