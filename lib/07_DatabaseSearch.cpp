#include "07_DatabaseSearch.h"
#include "Tools.h"
#include <algorithm>
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

inline std::unordered_map<std::string, std::vector<std::string::const_iterator>>
return_kmers(unsigned int kmer_size, const std::string &seq) {
  std::unordered_map<std::string, std::vector<std::string::const_iterator>>
      kmers;
  kmers.reserve(seq.size() - kmer_size);

  for (int i = 0; i <= seq.size() - kmer_size; i++) {
    std::string ks = seq.substr(i, kmer_size);
    auto pks = kmers.find(ks);
    if (pks == kmers.end())
      kmers[ks] = std::vector<std::string::const_iterator>{seq.begin() + i};
    else
      kmers[ks].emplace_back(seq.begin() + i);
  }

  return kmers;
};

std::vector<BLAST_hit> find_hits_seq(
    std::unordered_map<std::string, std::vector<std::string::const_iterator>>
        &kmers,
    const std::string &seq, unsigned int kmer_size) {

  std::vector<BLAST_hit> hits;
  for (int i = 0; i <= seq.size() - kmer_size; i++) {
    std::string ks = seq.substr(i, kmer_size);
    auto pks = kmers.find(ks);
    if (pks != kmers.end()) {
      for (auto &x : pks->second) {
        hits.push_back(
            BLAST_hit(x, seq.begin() + i, seq.begin() + i + kmer_size));
      }
    }
  }
  return hits;
}

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

void BLAST_db::find_sequence(const std::string &query, unsigned int ksize) {
  std::unordered_map<std::string, std::vector<std::string::const_iterator>>
      kmers = return_kmers(ksize, query);

  for (auto &x : db) {
    std::vector<BLAST_hit> hits = find_hits_seq(kmers, x.second, ksize);
    if (hits.size() != 0) {
      std::cout << x.first << " = " << hits.size() << std::endl;
      // std::cout << x.second << std::endl;
      // print_pattern_hits(x.second, std::string(hits[0].q, hits[0].q + ksize),
      //                   hits[0].b);
      // print_pattern_hits(x.second, std::string(hits[1].q, hits[1].q + ksize),
      //                   hits[1].b);
      // print_pattern_hits(x.second, std::string(hits[2].q, hits[2].q + ksize),
      //                   hits[2].b);
      // for (auto &y : hits) {
      // std::cout << std::string(y.q, y.q + ksize) << "\n"
      //          << std::string(y.b, y.e) << "\n"
      //          << std::endl;
      //}
    }
  }
}
