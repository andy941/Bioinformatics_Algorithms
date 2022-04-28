#pragma once

#include "06_SequenceAlignment.h"
#include <string>
#include <unordered_map>
#include <vector>

inline void print_fasta(const std::string &fasta, unsigned int linesize = 80);

std::unordered_map<std::string, int> read_blosum(const std::string &file,
                                                 std::string &alphabet);

std::vector<std::pair<std::string, std::string>>
read_fasta(const std::string &filename);

class BLAST_db {

  std::vector<std::pair<std::string, std::string>> db;
  std::unordered_map<std::string, int> sm;
  std::vector<alignment> bhits_aln;
  unsigned int kmer_size{0};
  std::string letters;

private:
  std::unordered_map<std::string, std::vector<unsigned int>>
  extract_kmers(const std::string &seq);

public:
  BLAST_db() = delete;
  BLAST_db(const std::string &filename_db, const std::string filename_blosum,
           unsigned int k = 3);
  BLAST_db(const std::string &filename_db, const int match, const int mismatch,
           const std::string &alphabet, unsigned int k = 11);
};
