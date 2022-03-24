#include "07_DatabaseSearch.h"
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

inline void print_fasta(const std::string &fasta, unsigned int linesize) {
  int limit = 0;
  for (auto &x : fasta) {
    if (limit < linesize) {
      std::cout << x;
      limit++;

    } else {
      std::cout << x << std::endl;
      limit = 0;
    }
  }
};

std::vector<std::pair<std::string, std::string>>
read_fasta(const std::string &filename) {
  std::vector<std::pair<std::string, std::string>> fa_seqs;
  std::string name;
  std::string seq;
  std::ifstream ifs{filename};
  if (!ifs)
    std::cerr << "Can't find file!" << std::endl;

  for (std::string line; !ifs.eof(); getline(ifs, line)) {
    std::istringstream iss{line};
    char c;
    iss >> c;
    if (c == '>' && name != "") {
      fa_seqs.push_back(std::make_pair(name, seq));
      name = "";
      seq = "";
    } else if (c == '>') {
      getline(iss, name);

    } else {
      seq += line;
    }
  }
  return fa_seqs;
};

inline std::vector<std::string> return_kmers(unsigned int kmer_size,
                                             const std::string &seq){

};

/*
 * Implement a very simple BLAST database (from python examples in the book
 * chapter)
 */

BLAST_db::BLAST_db(const std::string &filename_db,
                   const std::string &filename_blosum) {
  db = read_fasta(filename_db);
  sm = read_submat(filename_blosum);
};

BLAST_db::BLAST_db(const std::string &filename_db, const int &match,
                   const int &mismatch, const std::string &alphabet) {
  db = read_fasta(filename_db);
  sm = create_submat(match, mismatch, alphabet);
};

void BLAST_db::build_index(unsigned short int kmer_size, int threshold,
                           int gap_cost){

};

void BLAST_db::find_sequence(const std::string &query){};
