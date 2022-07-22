#include "08_MultipleSequenceAlignment.h"
#include "06_SequenceAlignment.h"
#include "07_DatabaseSearch.h"
#include <algorithm>
#include <array>

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

void msa::multiple_alignment::add_sequence(std::string &name,
                                           std::string &seq) {
  if (seq.size() != aln_pos.size()) {
    std::cerr << "trying to insert a sequence of size " << seq.size()
              << " in alignment of size " << aln_pos.size() << std::endl;
    return;
  }
  seq_names.push_back(name);
  for (size_t i = 0; i < aln_pos.size(); i++) {
    aln_pos[i].push_back(seq[i]);
  }
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
      for (size_t j = 0; j < std::min(line_size, aln_pos.size() - i); j++) {
      }
    }
  }
};

msa::msa::msa(const int match, const int mismatch, const std::string &alphabet)
    : nw{match, mismatch, alphabet} {
  sm = create_submat(match, mismatch, alphabet);
};
msa::msa::msa(const std::string &file) : nw{file} { sm = read_submat(file); };
