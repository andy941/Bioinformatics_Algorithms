#include "07_DatabaseSearch.h"
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

std::vector<std::pair<const std::string, const std::string>>
read_fasta(const std::string &filename) {
  std::vector<std::pair<const std::string, const std::string>> fa_seqs;
  std::string name;
  std::string seq;
  std::ifstream ifs{filename};
  if (!ifs)
    std::cerr << "Can't find file!" << std::endl;
  char c;
  while (!ifs.eof()) {
    ifs >> c;
    if (c != '>')
      std::cerr << "File not in FASTA format" << std::endl;
    getline(ifs, name);
    ifs >> c;
    while (c != '>') {
      seq.push_back(c);
    }
    fa_seqs.push_back(std::make_pair(name, seq));
  }
  return fa_seqs;
};

inline void print_fasta(const std::string &fasta, unsigned int linesize = 80) {
  int limit = 0;
  for (auto &x : fasta) {
    if (limit < 80)
      std::cout << x;
    else
      std::cout << x << std::endl;
  }
};

BLAST_db::BLAST_db(const std::string &filename_db,
                   const std::string &filename_blosum, unsigned int gap_cost);
BLAST_db::BLAST_db(const std::string &filename_db, const int &match,
                   const int &mismatch, const std::string &alphabet,
                   unsigned int gap_cost);
