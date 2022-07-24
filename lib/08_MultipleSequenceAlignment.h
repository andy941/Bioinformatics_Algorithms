#pragma once

#include "06_SequenceAlignment.h"
#include <Eigen/Core>
#include <string>
#include <vector>

namespace msa {

typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> similarity_matrix;

int sum_of_pairs(const std::string &str,
                 std::unordered_map<std::string, int> &sm, const int gap_cost);

struct multiple_alignment {
  std::vector<std::string> aln_pos{};
  std::vector<std::string> seq_names{};
  int score{0};

  multiple_alignment() = default;
  void add_sequence(const std::string &name, const std::string &seq);
  int get_score(std::unordered_map<std::string, int> &sm, const int gap_cost);
  std::string get_consensus();
  void print(size_t line_size = 80);
};

class msa {
  multiple_alignment m_aln{};
  std::unordered_map<std::string, int> sm;
  needleman_Wunsch nw;
  int gap_cost{0};

public:
  msa(const int match, const int mismatch, const std::string &alphabet);
  msa(const std::string &file);

public:
  void align_sequences(const std::vector<std::pair<std::string, std::string>> &,
                       int gap_cost);
  void print(size_t line_size = 80);
};

similarity_matrix
pairwise_aln_scores(const std::vector<std::pair<std::string, std::string>> &);

std::vector<int> upgma_order(similarity_matrix &);

} // namespace msa
