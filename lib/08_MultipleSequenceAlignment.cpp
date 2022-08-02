#include "08_MultipleSequenceAlignment.h"
#include "06_SequenceAlignment.h"
#include "07_DatabaseSearch.h"
#include <algorithm>
#include <array>
#include <iomanip>
#include <list>
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
      sim_mat(j, i) = -nw.get_score();
    }
  }
  return sim_mat;
}

std::array<size_t, 2> msa::upgma_min(similarity_matrix &sim_mat) {
  std::array<size_t, 2> rowcol{1, 0};
  unsigned int min{static_cast<unsigned int>(sim_mat(rowcol[0], rowcol[1]))};
  for (size_t i = 1; i < sim_mat.rows(); i++) { // rows
    for (size_t j = 0; j < i; j++) {            // columns
      if (min > sim_mat(i, j)) {
        min = sim_mat(i, j);
        rowcol[0] = i;
        rowcol[1] = j;
      }
    }
  }
  return rowcol;
}

void msa::upgma_reduce(similarity_matrix *sim_mat,
                       similarity_matrix *condensed_mat,
                       std::array<size_t, 2> &rowcol,
                       std::list<size_t> &labels){

};

std::vector<size_t> msa::upgma_order(similarity_matrix sim_mat) {
  std::vector<size_t> res;
  similarity_matrix sim_condensed;
  sim_condensed = similarity_matrix::Zero(sim_mat.rows(), sim_mat.cols());
  similarity_matrix *psm = &sim_mat;
  similarity_matrix *psc = &sim_condensed;
  similarity_matrix *tmp = nullptr;
  std::list<size_t> labels; // "condensed" cols/rows are < 0
  for (size_t i = 0; i < sim_mat.rows(); i++) {
    labels.push_back(i);
  }
  while (res.size() < sim_mat.rows()) {
    std::array<size_t, 2> rowcol{std::move(upgma_min(sim_mat))};
    if (rowcol[0] != -1)
      res.push_back(rowcol[0]);
    if (rowcol[1] != -1)
      res.push_back(rowcol[1]);
    upgma_reduce(psm, psc, rowcol, labels);
    // swap pointers and zero matrix for next cycle
    tmp = psm;
    psm = psc;
    psc = tmp;
    *psc = similarity_matrix::Zero(sim_mat.rows(), sim_mat.cols());
  }
  return res;
};

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
  std::vector<size_t> order{std::move(upgma_order(sim_mat))};
  std::cout << sim_mat << std::endl;
}
