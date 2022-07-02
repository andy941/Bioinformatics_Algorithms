#include "07_DatabaseSearch_improved.h"
#include "06_SequenceAlignment.h"
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

int BLAST_db::match_score(const std::string &s1, const std::string &s2) {
  int score{0};
  if (s1.size() != s2.size())
    std::cerr << "Strings not the same length in match_score()";

  for (int i = 0; i < s1.size(); i++) {
    std::string s{s1[i], s2[i]};
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
  kmers.reserve(letters.size() * kmer_size * query_sequence.size());

  for (unsigned int i = 0; i <= seq.size() - kmer_size; i++) {
    std::string ks = seq.substr(i, kmer_size);
    std::string ks_orig = ks;
    for (auto c_it = ks.begin(); c_it < ks.end(); c_it++) {
      char c_orig = *c_it;
      for (auto &l : letters) {
        *c_it = l;
        if (match_score(ks, ks_orig) > min_score) {
          auto pks = kmers.find(ks);
          if (pks == kmers.end())
            kmers[ks] = std::vector<unsigned int>{i};
          else if (pks->second.back() != i)
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
  std::cout << "size = " << mat.size() << " max_size = "
            << kmer_size * letters.size() * query_sequence.size() << std::endl;
  std::cout << "letter space = " << letters << " (" << letters.size() << ")"
            << std::endl;
}

Matrix_hits BLAST_db::find_hits(
    std::unordered_map<std::string, std::vector<unsigned int>> &kmers,
    const unsigned int db_position) {

  auto &seq = db[db_position];
  unsigned int rows = query_sequence.size() - kmer_size + 1;
  unsigned int cols = seq.second.size() - kmer_size + 1;

  Matrix_hits hmat(rows, cols);
  hmat = Matrix_hits::Zero(rows, cols);

  for (auto i = 0; i < cols; i++) {
    std::string ks{seq.second.begin() + i, seq.second.begin() + i + kmer_size};
    auto pks = kmers.find(ks);
    if (pks != kmers.end()) {
      for (auto &qhit : pks->second) {
        hmat(qhit, i) = true;
      }
    }
  }
  return hmat;
}

std::vector<BLAST_hit>
BLAST_db::collapse_hits(Matrix_hits &hmat, const unsigned int db_position,
                        const unsigned int collapse_limit) {

  std::vector<BLAST_hit> hits;

  auto &seq = db[db_position];
  int qstart = 0;
  int sstart = 0;
  unsigned int length = 0;
  unsigned int nohit_length = 0;

  for (unsigned int i = 0; i < hmat.cols(); i++) {
    unsigned int j = 0;
    qstart = j;
    sstart = i;
    length = nohit_length = 0;

    while (j < hmat.rows() && j + i < hmat.cols()) {
      const bool b = hmat(j, j + i);

      if (b == true && length == 0) {
        qstart = j;
        sstart = j + i;
        nohit_length = 0;
      }
      length += b;
      nohit_length += !b;
      nohit_length *= !b;

      if (nohit_length > collapse_limit && length != 0) {
        hits.push_back(BLAST_hit(qstart, sstart,
                                 j - qstart - nohit_length + kmer_size,
                                 db_position));
        nohit_length = 0;
        length = 0;
      }
      j++;
    }
    if (length != 0) {
      hits.push_back(BLAST_hit(qstart, sstart,
                               j - qstart - nohit_length + kmer_size - 1,
                               db_position));
    }
  }

  for (unsigned int i = 1; i < hmat.rows(); i++) {
    unsigned int j = 0;
    qstart = i;
    sstart = j;
    length = nohit_length = 0;

    while (j < hmat.cols() && j + i < hmat.rows()) {
      const bool b = hmat(j + i, j);

      if (b == true && length == 0) {
        qstart = j + i;
        sstart = j;
        nohit_length = 0;
      }
      length += b;
      nohit_length += !b;
      nohit_length *= !b;

      if (nohit_length > collapse_limit && length != 0) {
        hits.push_back(BLAST_hit(qstart, sstart,
                                 j - sstart - nohit_length + kmer_size,
                                 db_position));
        nohit_length = 0;
        length = 0;
      }
      j++;
    }
    if (length != 0) {
      hits.push_back(BLAST_hit(qstart, sstart,
                               j - sstart - nohit_length + kmer_size - 1,
                               db_position));
    }
  }

  return hits;
}

// Keep extending until the running value of the score dips below 0.
void BLAST_db::extend_hits(std::vector<BLAST_hit> hits) {
  for (auto &hit : hits) {
    int score = match_score(query_sequence[], const std::string &s2);
    auto &seq = db[hit.db_position];
    while (hit.query_start > 0 && hit.seq_start > 0) {
      char c1 = query_sequence[hit.query_start - 1];
      char c2 = seq.second[hit.seq_start - 1];
      if (score += score_pos(c1, c2, sm, 0) < 0) {
        break;
      }
      hit.query_start--;
      hit.seq_start--;
      hit.length++;
    }
    score = 0;
    while (hit.length + hit.query_start < query_sequence.length() &&
           hit.length + hit.seq_start < seq.second.length()) {
      char c1 = query_sequence[hit.query_start + hit.length];
      char c2 = seq.second[hit.seq_start + hit.length];
      if (score += score_pos(c1, c2, sm, 0) < 0) {
        break;
      }
      hit.length++;
    }
  }
}
void BLAST_db::compute_Escore(BLAST_hit &hit) {}

std::vector<BLAST_hit> BLAST_db::get_HSP(Matrix_hits &hmat,
                                         const unsigned int db_position,
                                         const unsigned int collapse_limit) {

  std::vector<BLAST_hit> hits =
      collapse_hits(hmat, db_position, collapse_limit);
  extend_hits(hits); // Could be improved (see boundary alignment)

  return hits;
}

void BLAST_db::blast_sequence(std::string &query) {
  std::vector<BLAST_hit> hits;
  query_sequence = query;
  auto kmers_map = extract_kmers(query);
  std::string db_seq_name = db[3].first;
  std::string db_seq = db[3].second;
  auto hits_mat = find_hits(kmers_map, 3);
  hits = get_HSP(hits_mat, 3, 3);

  std::cout << hits_mat << std::endl;
}
