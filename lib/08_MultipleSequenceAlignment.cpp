#include "08_MultipleSequenceAlignment.h"
#include "06_SequenceAlignment.h"
#include "07_DatabaseSearch.h"
#include <algorithm>
#include <array>
#include <iomanip>
#include <vector>

void msa::msa::print(size_t line_size) { m_aln.print(line_size); };

msa::msa::msa(const int match, const int mismatch, const std::string &alphabet)
    : nw{match, mismatch, alphabet} {
  sm = create_submat(match, mismatch, alphabet);
};
msa::msa::msa(const std::string &file) : nw{file} { sm = read_submat(file); };

msa::similarity_matrix msa::pairwise_aln_scores(
    const std::vector<std::pair<std::string, std::string>> &seqs) {
  needleman_Wunsch nw{0, -1, "ATCG"};
  similarity_matrix sim_mat;
  sim_mat = similarity_matrix::Zero(seqs.size(), seqs.size());
  for (size_t i = 0; i < seqs.size(); i++) {
    for (size_t j = i + 1; j < seqs.size(); j++) {
      nw.align_sequences(seqs[i].second, seqs[j].second, -1);
      nw.trace_back();
      nw.print();
      sim_mat(j, i) = -nw.get_score();
    }
  }
  return sim_mat;
}

std::vector<int> upgma_order(msa::similarity_matrix &);

void msa::msa::align_sequences(
    const std::vector<std::pair<std::string, std::string>> &sequences,
    int gap_cost) {
  gap_cost = gap_cost;
  for (auto &seq : sequences) {
    m_aln.add_sequence(seq.first, seq.second);
  }
  std::cout << m_aln.get_consensus() << std::endl;
  std::cout << m_aln.get_score(sm, gap_cost) << std::endl;
  auto sim_mat = pairwise_aln_scores(sequences);
  std::cout << sim_mat << std::endl;
}
