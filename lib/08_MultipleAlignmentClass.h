#pragma once

#include <string>
#include <unordered_map>
#include <vector>

namespace msa {

int sum_of_pairs(const std::string &str,
                 std::unordered_map<std::string, int> &sm, const int gap_cost);

struct multiple_alignment {
  std::vector<std::string> aln_col{};
  std::vector<std::string> seq_names{};

  multiple_alignment() = default;
  void add_sequence(const std::string &name, const std::string &seq);
  int get_score(std::unordered_map<std::string, int> &sm, const int gap_cost);
  std::string get_consensus();
  void print(size_t line_size = 80);
};

} // namespace msa
