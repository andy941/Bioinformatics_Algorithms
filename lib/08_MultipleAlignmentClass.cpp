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
  if (aln_col.size() == 0) {
    seq_names.push_back(name);
    for (size_t i = 0; i < seq.size(); i++) {
      aln_col.push_back(std::string{seq[i]});
    }
    return;
  }
  if (seq.size() != aln_col.size()) {
    std::cerr << "trying to insert a sequence of size " << seq.size()
              << " in alignment of size " << aln_col.size() << std::endl;
    return;
  }
  for (size_t i = 0; i < seq.size(); i++) {
    aln_col[i].push_back(seq[i]);
  }
  seq_names.push_back(name);
};

int msa::multiple_alignment::get_score(std::unordered_map<std::string, int> &sm,
                                       const int gap_cost) {
  int score = 0;
  if (aln_col.empty())
    return 0;
  for (auto &x : aln_col) {
    score += sum_of_pairs(x, sm, gap_cost);
  }
  return score;
};

std::string msa::multiple_alignment::get_consensus() {
  std::string cons_aln{};
  if (aln_col.empty())
    return cons_aln;
  for (auto &x : aln_col) {
    std::string s_x;
    s_x.reserve(x.size());
    std::string letters;
    std::vector<unsigned int> occurrs;
    occurrs.reserve(x.size());
    std::partial_sort_copy(x.begin(), x.end(), s_x.begin(), s_x.end());
    letters.push_back(*x.begin());
    occurrs.push_back(1);
    for (auto c = x.begin() + 1; c < x.end(); c++) {
      if (*c == '-')
        continue;
      if (*c != letters.back()) {
        letters.push_back(*c);
        occurrs.push_back(0);
      }
      occurrs.back()++;
    }
    char l = letters.front();
    unsigned int count = occurrs.front();
    for (size_t i = 1; i < letters.size(); i++) {
      if (occurrs[i] > count) {
        l = letters[i];
      }
    }
    cons_aln.push_back(l);
  }
  return cons_aln;
};

void msa::multiple_alignment::print(size_t line_size) {
  for (size_t i = 0; i < aln_col.size(); i += line_size) {
    for (size_t k = 0; k < seq_names.size(); k++) {
      if (seq_names[k].size() > 10)
        std::cout << seq_names[k].substr(0, 10) << '\t';
      else
        std::cout << std::setw(10) << seq_names[k] << '\t';

      for (size_t j = i; j < std::min(i + line_size, aln_col.size()); j++) {
        std::cout << aln_col[j][k];
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
};
