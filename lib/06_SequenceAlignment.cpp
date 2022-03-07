#include "06_SequenceAlignment.h"
#include <ctype.h>
#include <fstream>
#include <iostream>
#include <istream>
#include <sstream>
#include <vector>

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
