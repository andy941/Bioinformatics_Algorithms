#pragma once

#include "06_SequenceAlignment.h"
#include <Eigen/Core>
#include <string>
#include <unordered_map>
#include <vector>

typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> Matrix_hits;

inline void print_fasta(const std::string &fasta, unsigned int linesize = 80);

std::vector<std::pair<std::string, std::string>>
read_fasta(const std::string &filename);

std::unordered_map<std::string, int> read_blosum(const std::string &file,
                                                 std::string &alphabet);

struct BLAST_hit {
  double E_score{0};
  alignment aln;
  unsigned int query_start{0};
  unsigned int seq_start{0};
  unsigned int length{0};
  unsigned int db_position{0};
  int score{0};

  BLAST_hit() = default;
  BLAST_hit(unsigned int qs, unsigned int ss, unsigned int l, unsigned int dp)
      : query_start{qs}, seq_start{ss}, length{l}, db_position{dp} {};
};

double compute_bitscore(const int score, const double lambda, const double K);
void compute_Escore(BLAST_hit &hit, double bitscore, unsigned long int db_size);

class BLAST_db {

  double lambda{1.94643};    // From search on BLAST of Arabidopsis sequence
  double K{0.805949};        // From search on BLAST of Arabidopsis sequence
  double Escore_limit{0.01}; // From search on BLAST of Arabidopsis sequence
  std::vector<std::pair<std::string, std::string>> db;
  std::unordered_map<std::string, int> sm;
  std::vector<BLAST_hit> bhits_result;
  unsigned int kmer_size{0};
  unsigned int min_score{0};
  std::string letters;
  std::string query_sequence;
  unsigned long int db_size{0};

private:
  int match_score(const std::string &s1, const std::string &s2);
  std::unordered_map<std::string, std::vector<unsigned int>>
  extract_kmers(const std::string seq);
  void
  print_kmers(std::unordered_map<std::string, std::vector<unsigned int>> &mat);
  Matrix_hits
  find_hits(std::unordered_map<std::string, std::vector<unsigned int>> &mat,
            const unsigned int db_position);
  std::vector<BLAST_hit> collapse_hits(Matrix_hits &hmat,
                                       const unsigned int db_position,
                                       const unsigned int collapse_limit);
  void extend_hits(std::vector<BLAST_hit> &hits);
  std::vector<BLAST_hit> get_HSP(Matrix_hits &hmat, unsigned int db_position,
                                 const unsigned int collapse_limit);
  void compute_alignment(BLAST_hit &hit);

public:
  BLAST_db() = delete;
  BLAST_db(const std::string &filename_db, const std::string filename_blosum,
           unsigned int k = 3, unsigned int ms = 12);
  BLAST_db(const std::string &filename_db, const int match, const int mismatch,
           const std::string &alphabet, unsigned int k = 11,
           unsigned int ms = 4);

public:
  void blast_sequence(std::string &seq, double Escore_limit = 0.01);
  void print_report(unsigned int results = 10);
};
