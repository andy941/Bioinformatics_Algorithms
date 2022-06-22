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
  unsigned int query_start{0};
  unsigned int seq_start{0};
  unsigned int length{0};
  std::string seq_name;
  alignment aln;
  unsigned int score{0};
  unsigned int E_score{0};
  BLAST_hit() = default;
  BLAST_hit(unsigned int qs, unsigned int ss, unsigned int l, std::string sn)
      : query_start{qs}, seq_start{ss}, length{l}, seq_name{sn} {};
};

class BLAST_db {

  std::vector<std::pair<std::string, std::string>> db;
  std::unordered_map<std::string, int> sm;
  std::unordered_map<std::string,
                     std::pair<unsigned int, std::vector<unsigned int>>>
      kmers;
  std::vector<alignment> bhits_aln;
  unsigned int kmer_size{0};
  unsigned int min_score{0};
  std::string letters;
  std::string query_sequence;

private:
  int match_score(const std::string &s1, const std::string &s2);
  std::unordered_map<std::string, std::vector<unsigned int>>
  extract_kmers(const std::string seq);
  void
  print_kmers(std::unordered_map<std::string, std::vector<unsigned int>> &mat);
  Matrix_hits
  find_hits(std::unordered_map<std::string, std::vector<unsigned int>> &mat,
            std::string &seq, unsigned int len_db_seq);
  void collapse_hits(Matrix_hits &hmat);
  std::vector<BLAST_hit> get_HSP(const Matrix_hits &hmat,
                                 const std::pair<std::string, std::string> &seq,
                                 const unsigned int dist);

public:
  BLAST_db() = delete;
  BLAST_db(const std::string &filename_db, const std::string filename_blosum,
           unsigned int k = 3, unsigned int ms = 12);
  BLAST_db(const std::string &filename_db, const int match, const int mismatch,
           const std::string &alphabet, unsigned int k = 11,
           unsigned int ms = 4);

public:
  void blast_sequence(std::string &seq);
};
