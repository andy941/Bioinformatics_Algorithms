#pragma once

#include "06_SequenceAlignment.h"
#include "08_MultipleAlignmentClass.h"
#include <Eigen/Core>
#include <list>
#include <string>
#include <vector>

namespace msa {

typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> similarity_matrix;

class msa {
  multiple_alignment m_aln{};
  std::unordered_map<std::string, int> sm;
  needleman_Wunsch nw;
  int gap_cost{0};

public:
  msa(const int match, const int mismatch, const std::string &alphabet);
  msa(const std::string &file);

public:
  void
  align_sequences(const std::vector<std::pair<std::string, std::string>> &seqs,
                  int gap_cost);
  void print(size_t line_size = 80);
};

similarity_matrix pairwise_aln_scores(
    const std::vector<std::pair<std::string, std::string>> &seqs);

std::array<size_t, 2> upgma_min(similarity_matrix &sim_mat);
void upgma_reduce(similarity_matrix *sim_mat, similarity_matrix *condensed_mat,
                  std::array<size_t, 2> &rowcol, std::list<size_t> &labels);
std::vector<size_t> upgma_order(similarity_matrix);

} // namespace msa
