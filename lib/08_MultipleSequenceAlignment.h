#pragma once

#include "06_SequenceAlignment.h"
#include <Eigen/Core>
#include <string>
#include <vector>

namespace MSA {

typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> similarity_matrix;

struct multiple_alignment {
  std::vector<std::pair<std::string, std::string>> alignments{};
  int score{0};

  multiple_alignment() = default;
  int get_score();
  std::string get_consensus();
  void print();
};

class MSA {
  std::vector<std::pair<std::string, std::string>> m_seq;
  multiple_alignment m_aln{};
  std::unordered_map<std::string, int> sm;
  needleman_Wunsch nw;
  int gap_cost;

public:
  MSA(const int &match, const int &mismatch, const std::string &alphabet);
  MSA(const std::string &file);

public:
  void align_sequences(const std::vector<std::pair<std::string, std::string>> &,
                       const int gap_cost);
};

similarity_matrix
upgma_matrix(const std::vector<std::pair<std::string, std::string>> &);

std::vector<int> upgma_order(similarity_matrix &);

} // namespace MSA
