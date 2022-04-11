#pragma once

#include "06_SequenceAlignment.h"
#include <string>
#include <unordered_map>
#include <vector>

inline void print_fasta(const std::string &fasta, unsigned int linesize = 80);

std::vector<std::pair<std::string, std::string>>
read_fasta(const std::string &filename);

inline std::unordered_map<std::string, std::vector<std::string::const_iterator>>
return_kmers(unsigned int kmer_size, const std::string &seq);

struct BLAST_hit {

  std::string::const_iterator q; // query hit
  std::string::const_iterator b; // seq hit
  std::string::const_iterator e; // seq hit
  std::string seq_name;
  std::string seq;
  alignment aln;
  int score{0};
  BLAST_hit(std::string::const_iterator qq, std::string::const_iterator bb,
            std::string::const_iterator ee, std::string sn, int sc)
      : q{qq}, b{bb}, e{ee}, seq_name{sn}, score{sc} {};
};

std::vector<BLAST_hit> find_hits_seq(
    std::unordered_map<std::string, std::vector<std::string::const_iterator>>
        &kmers,
    const std::string &seq, unsigned int kmer_size, unsigned int dist);

void extend_hits(const std::string &query, const std::string &seq,
                 std::vector<BLAST_hit> &hits, unsigned int kmer_size);

class BLAST_db {
  std::vector<std::pair<std::string, std::string>> db;
  std::unordered_map<std::string, int> sm;
  std::vector<BLAST_hit> bhits;
  smith_Waterman sW;
  std::vector<alignment> bhits_aln;

private:
  std::vector<BLAST_hit> find_hits_seq(
      std::unordered_map<std::string, std::vector<std::string::const_iterator>>
          &kmers,
      const std::pair<std::string, std::string> &seq, unsigned int kmer_size);
  void extend_hits(const std::string &query, const std::string &seq,
                   std::vector<BLAST_hit> &hits, unsigned int kmer_size);

public:
  BLAST_db() = delete;
  BLAST_db(const std::string &filename_db, const std::string &filename_blosum);
  BLAST_db(const std::string &filename_db, const int &match,
           const int &mismatch, const std::string &alphabet);

public:
  void find_sequence(const std::string &query, unsigned int ksize = 17);
  void print_report();
};
