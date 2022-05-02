#include "07_DatabaseSearch_improved.h"
#include "Tools.h"
#include <algorithm>
#include <fstream>
#include <numeric>
#include <ostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

std::unordered_map<std::string, int> read_blosum(const std::string &file,
                                                 std::string &alphabet) {

  std::unordered_map<std::string, int> submat;
  std::ifstream ifs{file};
  if (!ifs)
    std::cerr << "Can't find file!" << std::endl;
  std::string header;
  std::string header_line;
  getline(ifs, header_line);
  for (auto &x : header_line) {
    if (isalpha(x))
      header.push_back(x);
  }

  alphabet = header;
  submat.reserve(header.size() * header.size());

  auto x{header.begin()};
  for (std::string line; std::getline(ifs, line);) {
    std::istringstream iss{line};
    for (auto y : header) {
      std::string key{*x};
      key += y;
      int value;
      iss >> value;
      submat[key] = value;
      char comma;
      iss >> comma;
      if (comma != ',')
        std::cerr << "unexpected input, not a csv file." << std::endl;
    }
    advance(x, 1);
  }
  return submat;
};

void print_fasta(const std::string &fasta, unsigned int linesize) {
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

int BLAST_db::match_score(std::string &s) {
  int score{0};

  for (auto c_it = s.begin(); c_it < s.end(); c_it++) {
    std::string s{*c_it, *c_it};
    score += sm[s];
  }

  return score;
}

// BLAST_db -----------------------------------------------------------------

BLAST_db::BLAST_db(const std::string &filename_db,
                   const std::string filename_blosum, unsigned int k)
    : kmer_size{k} {
  db = read_fasta(filename_db);
  sm = read_blosum(filename_blosum, letters);
  extract_kmers(db[0].second);
};

BLAST_db::BLAST_db(const std::string &filename_db, const int match,
                   const int mismatch, const std::string &alphabet,
                   unsigned int k)
    : kmer_size{k} {
  db = read_fasta(filename_db);
  sm = create_submat(match, mismatch, alphabet);
  letters = alphabet;
};

std::unordered_map<std::string, std::vector<unsigned int>>
BLAST_db::extract_kmers(const std::string seq) {

  std::unordered_map<std::string, std::vector<unsigned int>> kmers;
  kmers.reserve(kmer_size * letters.size());

  for (unsigned int i = 0; i <= seq.size() - kmer_size; i++) {
    std::string ks = seq.substr(i, kmer_size);
    for (auto c_it = ks.begin(); c_it < ks.end(); c_it++) {
      char c_orig = *c_it;
      for (auto &l : letters) {
        *c_it = l;
        if (match_score(ks) > 12) {
          auto pks = kmers.find(ks);
          if (pks == kmers.end())
            kmers[ks] = std::vector<unsigned int>{i};
          else
            kmers[ks].push_back(i);
        }
      }
      *c_it = c_orig;
    }
  }

  return kmers;
};
