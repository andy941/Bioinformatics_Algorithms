#include "Tools.h"
#include <cstdint>
#include <iostream>
#include <iterator>
#include <unordered_map>
#include <vector>

using namespace std;

class BoyerMoore {

  string alphabet;
  string pattern;
  unordered_map<char, int> occ;
  vector<int> f;
  vector<int> s;

  void preprocess_bcr();
  void preprocess_gsr();

public:
  BoyerMoore(string alphabet, string pattern);
};

BoyerMoore::BoyerMoore(string alphabet, string pattern)
    : alphabet{alphabet}, pattern{pattern} {
  preprocess_bcr();
  // for (auto x : occ)
  //  cout << x.first << " " << x.second << endl;
  preprocess_gsr();
};

void BoyerMoore::preprocess_bcr() {
  for (char s : alphabet)
    occ[s] = -1;
  for (int i = 0; i < pattern.length(); i++)
    occ[pattern[i]] = i;
}

void BoyerMoore::preprocess_gsr() {
  f = vector<int>(pattern.length() + 1, 0);
  s = vector<int>(pattern.length() + 1, 0);
  int i = pattern.length();
  int j = pattern.length();
  f[i] = j;
  while (i > 0) {
  }
}

auto naive_find_all(const string &str, const string &pattern) {
  vector<string::const_iterator> res;
  if (pattern.begin() == pattern.end() || str.begin() == str.end()) {
    res.push_back(str.end());
    return res;
  }
  auto p = pattern.begin();
  for (auto x = str.begin(); x != str.end(); advance(x, 1)) {
    auto tmp = x;
    while (*tmp == *p) {
      if (p + 1 == pattern.end()) {
        res.push_back(x);
        break;
      }
      advance(p, 1);
      advance(tmp, 1);
    }
    p = pattern.begin();
  }
  return res;
}

int main() {
  const string seq = random_seq(80, 'D');
  const string pattern{"GA"};
  cout << "pattern :\t" << pattern << endl;
  cout << "sequence:\t" << seq << endl;

  auto pa = naive_find_all(seq, pattern);
  cout << "all     :\t";
  print_pattern_hits(seq, pattern, pa);

  BoyerMoore("ACGT", pattern);
}
