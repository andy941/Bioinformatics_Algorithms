#include "06_SequenceAlignment.h"
#include "Tools.h"
#include <algorithm>
#include <array>
#include <ctype.h>
#include <fstream>
#include <iostream>
#include <istream>
#include <numeric>
#include <sstream>
#include <vector>

// 03 Objective function
std::unordered_map<std::string, int>
create_submat(const int &match, const int &mismatch,
              const std::string &alphabet) {

  std::unordered_map<std::string, int> submat;
  for (auto &x1 : alphabet) {
    for (auto &x2 : alphabet) {
      std::string s;
      s += x1;
      s += x2;
      if (x1 == x2)
        submat[s] = match;
      else
        submat[s] = mismatch;
    }
  }
  return submat;
}

std::unordered_map<std::string, int> read_submat(const std::string &file) {

  std::unordered_map<std::string, int> submat;
  std::ifstream ifs{file};
  if (!ifs)
    std::cerr << "Can't find file!" << std::endl;
  std::vector<char> header;
  std::string header_line;
  getline(ifs, header_line);
  for (auto &x : header_line) {
    if (isalpha(x) || x == '*')
      header.push_back(x);
  }
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

int score_pos(const char &c1, const char &c2,
              std::unordered_map<std::string, int> &sm, const int &gap_cost) {
  if (c1 == '-' || c2 == '-')
    return gap_cost;
  std::string s{c1};
  s += c2;
  return sm[s];
}

int score_align(const std::string &s1, const std::string &s2,
                std::unordered_map<std::string, int> &sm, const int &gap_cost) {
  int res = 0;
  for (int i = 0; i < s1.size(); i++) {
    res += score_pos(s1[i], s2[i], sm, gap_cost);
  }
  return res;
}

int score_align_gapaff(const std::string &s1, const std::string &s2,
                       std::unordered_map<std::string, int> &sm,
                       const int &gap_cost, const int &gap_aff) {
  int res = 0;
  bool gap1{false};
  bool gap2{false};

  for (int i = 0; i < s1.size(); i++) {
    if (s1[i] == '-')
      if (gap1 == true)
        res += gap_aff;
      else {
        res += gap_cost;
        gap1 = true;
      }
    else if (s2[i] == '-')
      if (gap2 == true)
        res += gap_aff;
      else {
        res += gap_cost;
        gap2 = true;
      }
    else {
      res += score_pos(s1[i], s2[i], sm, gap_cost);
      gap1 = false;
      gap2 = false;
    }
  }
  return res;
}

// 04 NeedleMan-Wunsch
// NOTE: This is a very simple implementation that DOESN'T take into account
// multiple optimal alignments. The book doesn't delve into that for now. It
// would require more symbols in the T matrix and some simple path finder
// algorithm.

void alignment::print() {
  std::cout << "seq1: " << a << std::endl;
  std::cout << "      ";
  for (int i = 0; i < a.size(); i++) {
    if (a[i] == b[i])
      std::cout << '|';
    else if (a[i] != '-' && b[i] != '-')
      std::cout << '*';
    else
      std::cout << ' ';
  }
  std::cout << std::endl;
  std::cout << "seq2: " << b << std::endl;
};

void alignment::add(const char &ac, const char &bc) {
  a.push_back(ac);
  b.push_back(bc);
}

void alignment::flip() {
  reverse(a.begin(), a.end());
  reverse(b.begin(), b.end());
}

double alignment::identity() {
  double res = 0;
  int max_length = std::max(a.size(), b.size());
  for (int i = 0; i < a.size(); i++) {
    if (a[i] == b[i])
      res++;
  }
  return res / max_length;
}

std::string alignment::identical_subseq() {
  std::string res;
  res.reserve(std::max(a.size(), b.size()));
  for (int i = 0; i < a.size(); i++) {
    if (a[i] == b[i])
      res.push_back(a[i]);
  }
  return res;
};

needleman_Wunsch::needleman_Wunsch(const int &match, const int &mismatch,
                                   const std::string &alphabet) {
  sm = create_submat(match, mismatch, alphabet);
};

needleman_Wunsch::needleman_Wunsch(const std::string &file) {

  sm = read_submat(file);
};

needleman_Wunsch::~needleman_Wunsch() {
  delete[] S;
  delete[] T;
  S = nullptr;
  T = nullptr;
};

void needleman_Wunsch::align_sequences(const std::string &seq1,
                                       const std::string &seq2,
                                       const int &gap_cost) {
  reset();
  s1 = seq1;
  s2 = seq2;
  dim1 = s1.size() + 1; // because of the gap row and column
  dim2 = s2.size() + 1;
  S = new int[dim1 * dim2](); // initialize all zeroes
  T = new int[dim1 * dim2]();

  // NOTE: In T matrix 0 = go diagonal; 1 = go left; 2 = go up; start from
  // bottom right.

  // Gap row
  int cost = 0;
  for (int i = 0; i < dim2; i++) {
    S[i] = cost;
    T[i] = 1;
    cost += gap_cost;
  }

  // Gap column
  cost = 0;
  for (int i = 0; i < dim1; i++) {
    S[i * dim2] = cost;
    T[i * dim2] = 2;
    cost += gap_cost;
  }

  T[0] = 0; // top right has to be 0;

  // calculate the recursive costs for each cell of S and record the best path
  // choices in T.
  for (int i = 1; i < dim1; i++) { // start from 1, gap col/row calculated
    for (int j = 1; j < dim2; j++) {
      int pos = j + i * dim2;
      std::array<int, 3> arr;
      arr[0] =
          S[pos - 1 - dim2] + score_pos(s1[i - 1], s2[j - 1], sm, gap_cost);
      arr[1] = S[pos - 1] + gap_cost;
      arr[2] = S[pos - dim2] + gap_cost;
      int max = max_arr(arr);

      S[pos] = arr[max];
      T[pos] = max;
    }
  }

  best_score = S[dim2 * dim1 - 1];
};

// walk T from bottom right corner and reconstruct the alignment
void needleman_Wunsch::trace_back() {
  unsigned int px = dim2 - 1;
  unsigned int py = dim1 - 1;
  while (px != 0 || py != 0) {
    switch (T[px + py * dim2]) {
    case 0:
      aln.add(s1[py - 1], s2[px - 1]);
      px--;
      py--;
      break;
    case 1:
      aln.add('-', s2[px - 1]);
      px--;
      break;
    case 2:
      aln.add(s1[py - 1], '-');
      py--;
      break;
    default:
      std::cerr << "error tracing back the optimal alignment, check T matrix."
                << std::endl;
      throw 1;
    }
  }
  aln.flip();
};

void needleman_Wunsch::reset() {
  s1 = "";
  s2 = "";
  delete[] S;
  delete[] T;
  S = nullptr;
  T = nullptr;
  dim1 = 0;
  dim2 = 0;
  aln = alignment();
};

void needleman_Wunsch::print() {
  std::cout << "Sequences" << std::endl;
  std::cout << "seq1: " << s1 << std::endl;
  std::cout << "seq2: " << s2 << std::endl;
  std::string gap_s1 = '-' + s1;
  std::string gap_s2 = '-' + s2;
  std::cout << std::endl;
  std::cout << "Alignment" << std::endl;
  aln.print();
  std::cout << std::endl;
  std::cout << "S matrix" << std::endl;
  if (S != nullptr)
    print_matrix(gap_s1, gap_s2, S);
  std::cout << std::endl;
  std::cout << "T matrix" << std::endl;
  if (T != nullptr)
    print_matrix(gap_s1, gap_s2, T);
};

// 05 Smith-Waterman
// This is a variant of the algorithm tailored for local alignments.
void smith_Waterman::align_sequences(const std::string &seq1,
                                     const std::string &seq2,
                                     const int &gap_cost) {
  reset();
  s1 = seq1;
  s2 = seq2;
  dim1 = s1.size() + 1; // because of the gap row and column
  dim2 = s2.size() + 1;
  S = new int[dim1 * dim2](); // initialize all zeroes
  T = new int[dim1 * dim2]();

  // NOTE: In T matrix 0 = go diagonal; 1 = go left; 2 = go up; 3 = stop.
  // Gap row terminators
  int cost = 0;
  for (int i = 0; i < dim2; i++) {
    T[i] = 3;
  }
  // Gap column terminators
  for (int i = 0; i < dim1; i++) {
    T[i * dim2] = 3;
  }

  // calculate the recursive costs for each cell of S and record the best path
  // choices in T. The local alignment requires to evaluate 0 as alignment
  // termination which will be encoded as '3' in the T matrix.
  for (int i = 1; i < dim1; i++) { // start from 1, gap col/row calculated
    for (int j = 1; j < dim2; j++) {
      int pos = j + i * dim2;
      std::array<int, 4> arr;
      arr[0] =
          S[pos - 1 - dim2] + score_pos(s1[i - 1], s2[j - 1], sm, gap_cost);
      arr[1] = S[pos - 1] + gap_cost;
      arr[2] = S[pos - dim2] + gap_cost;
      arr[3] = 0;
      int max = max_arr(arr);

      S[pos] = arr[max];
      T[pos] = max;
    }
  }
};

// walk T from bottom right corner and reconstruct the alignment
void smith_Waterman::trace_back() {
  int pos = 0;
  int score_max = 0;
  for (int i = 0; i < dim1 * dim2; i++) { // This whole thing is ugly
    if (S[i] >= score_max) {
      pos = i;
      score_max = S[i];
    }
  }
  unsigned int px = pos % dim2;
  unsigned int py = pos / dim2;

  while (T[px + py * dim2] != 3) {
    switch (T[px + py * dim2]) {
    case 0:
      aln.add(s1[py - 1], s2[px - 1]);
      px--;
      py--;
      break;
    case 1:
      aln.add('-', s2[px - 1]);
      px--;
      break;
    case 2:
      aln.add(s1[py - 1], '-');
      py--;
      break;
    default:
      std::cerr << "error tracing back the optimal alignment, check T matrix."
                << std::endl;
      throw 1;
    }
  }
  aln.flip();
};
