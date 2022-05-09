#include "07_DatabaseSearch_improved.h"
#include "Eigen/src/Core/Matrix.h"
#include "Eigen/src/Core/util/Constants.h"
#include "Tools.h"
#include <algorithm>
#include <fstream>
#include <iostream>
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
    if (x != ',')
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
    if (line == "")
      continue;
    if (name == "") {
      std::cout << line << std::endl;
      name = std::string{line.begin() + 1, line.end()};
      continue;
    }
    char c = line[0];
    if (c == '>') {
      fa_seqs.push_back(std::make_pair(name, seq));
      name = std::string{line.begin() + 1, line.end()};
      seq = "";
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
                   const std::string filename_blosum, unsigned int k,
                   unsigned int ms)
    : kmer_size{k}, min_score(ms) {
  db = read_fasta(filename_db);
  sm = read_blosum(filename_blosum, letters);
}

BLAST_db::BLAST_db(const std::string &filename_db, const int match,
                   const int mismatch, const std::string &alphabet,
                   unsigned int k, unsigned int ms)
    : kmer_size{k}, min_score(ms) {
  db = read_fasta(filename_db);
  sm = create_submat(match, mismatch, alphabet);
  letters = alphabet;
};

std::unordered_map<std::string, std::vector<unsigned int>>
BLAST_db::extract_kmers(const std::string seq) {

  std::unordered_map<std::string, std::vector<unsigned int>> kmers;
  kmers.reserve(pow(letters.size(), kmer_size));

  for (unsigned int i = 0; i <= seq.size() - kmer_size; i++) {
    std::string ks = seq.substr(i, kmer_size);
    for (auto c_it = ks.begin(); c_it < ks.end(); c_it++) {
      char c_orig = *c_it;
      for (auto &l : letters) {
        *c_it = l;
        if (match_score(ks) > min_score) {
          auto pks = kmers.find(ks);
          if (pks == kmers.end())
            kmers[ks] = std::vector<unsigned int>{i};
          else
            pks->second.push_back(i);
        }
      }
      *c_it = c_orig;
    }
  }
  return kmers;
};

void BLAST_db::print_kmers(
    std::unordered_map<std::string, std::vector<unsigned int>> &mat) {
  for (auto &x : mat) {
    std::cout << x.first << " = ";
    for (auto &y : x.second) {
      std::cout << y << ", ";
    }
    std::cout << std::endl;
  }
  std::cout << "size = " << mat.size()
            << " max_size = " << pow(kmer_size, letters.size()) << std::endl;
  std::cout << "letter space = " << letters << " (" << letters.size() << ")"
            << std::endl;
}

Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic> BLAST_db::find_hits(
    std::unordered_map<std::string, std::vector<unsigned int>> &kmers,
    std::string &seq, unsigned int len_query) {

  Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic> hmat(
      len_query, seq.length());

  for (auto i = 0; i <= seq.size() - kmer_size; i++) {
    std::string ks{seq.begin() + i, seq.begin() + i + kmer_size};
    auto pks = kmers.find(ks);
    if (pks != kmers.end()) {
      for (auto &qhit : pks->second) {
        hmat(qhit, i) = true;
      }
    }
  }
  return hmat;
}

void BLAST_db::blast_sequence(std::string &query) {
  auto kmers_map = extract_kmers(query);
  std::string db_seq_name = db[11].first;
  std::string db_seq = db[11].second;
  auto hits_mat = find_hits(kmers_map, db_seq, query.size());
  std::cout << hits_mat << std::endl;
  std::cout << hits_mat.size() << " == " << query.size() * db_seq.size()
            << std::endl;
  std::cout << "ID = " << db_seq_name << std::endl;
  std::cout << db_seq << std::endl;
  std::cout << "Query" << std::endl;
  std::cout << query << std::endl;

  // print_kmers(mat);
}
