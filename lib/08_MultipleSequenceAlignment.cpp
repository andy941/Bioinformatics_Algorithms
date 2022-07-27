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

void msa::msa::align_sequences(
    const std::vector<std::pair<std::string, std::string>> &sequences,
    int gap_cost) {
  gap_cost = gap_cost;
  for (auto &seq : sequences) {
    m_aln.add_sequence(seq.first, seq.second);
  }
  std::cout << m_aln.get_consensus() << std::endl;
  std::cout << m_aln.get_score(sm, gap_cost) << std::endl;
}
