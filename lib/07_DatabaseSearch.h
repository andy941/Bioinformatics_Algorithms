#include "06_SequenceAlignment.h"
#include <string>
#include <unordered_map>
#include <vector>

inline void print_fasta(const std::string &fasta, unsigned int linesize = 80);

std::vector<std::pair<std::string, std::string>>
read_fasta(const std::string &filename);

inline std::vector<std::string> return_kmers(unsigned int kmer_size,
                                             const std::string &seq);

struct BLAST_hit {
  std::string name;
  std::string seq;
  alignment aln;
  unsigned int begin{0};
  unsigned int end{0};
  BLAST_hit() = default;
};

class BLAST_db {
  std::vector<std::pair<std::string, std::string>> db;
  std::unordered_map<std::string, int> sm;
  std::vector<std::string> qkmers;
  std::vector<BLAST_hit> bhits;
  unsigned int ksize{17};

protected:
  void build_index(unsigned short int kmer_size, int threshold, int gap_cost);

public:
  BLAST_db() = delete;
  BLAST_db(const std::string &filename_db, const std::string &filename_blosum);
  BLAST_db(const std::string &filename_db, const int &match,
           const int &mismatch, const std::string &alphabet);

public:
  void find_sequence(const std::string &query);
};
