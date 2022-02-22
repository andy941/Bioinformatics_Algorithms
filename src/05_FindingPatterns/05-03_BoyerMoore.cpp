#include "Tools.h"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <unordered_map>
#include <vector>

using namespace std;

class DFA {

  string alphabet;
  string pattern;
  unordered_map<char, int> occ;
  vector<int> f;
  vector<int> s;

public:
  void preprocess_bcr();
  void preprocess_gsr();
  auto find_all(const string &text);

public:
  DFA(string alphabet, string pattern);
};

DFA::DFA(string alphabet, string pattern)
    : alphabet{alphabet}, pattern{pattern} {
  preprocess_bcr();
  preprocess_gsr();
};

void DFA::preprocess_bcr() {
  for (char s : alphabet)
    occ[s] = -1;
  for (int i = 0; i < pattern.length(); i++)
    occ[pattern[i]] = i;
}

void DFA::preprocess_gsr() {
  f = vector<int>(pattern.length() + 1, 0);
  s = vector<int>(pattern.length() + 1, 0);
  int i = pattern.length();
  int j = pattern.length() + 1;
  f[i] = j;
  while (i > 0) {
    while (j <= pattern.length() && pattern[i - 1] != pattern[j - 1]) {
      if (s[j] == 0)
        s[j] = j - 1;
      j = f[j];
    }
    --i;
    --j;
    f[i] = j;
  }
  j = f[0];
  for (int i = 0; i < pattern.length(); i++) {
    if (s[i] == 0)
      s[i] = j;
    if (i == j)
      j = f[j];
  }
}

auto DFA::find_all(const string &text) {

  Timer t;
  vector<string::const_iterator> res;
  int i = 0;
  int j = pattern.length() - 1;
  char c;

  while (i <= text.length() - pattern.length()) {

    j = pattern.length() - 1;
    while (j >= 0 && pattern[j] == text[j + i]) {
      --j;
    }
    if (j < 0) {
      res.push_back(text.begin() + i);
      i += s[0];
    } else {
      c = text[j + i];
      i += max(s[j + 1], j - occ[c]);
    }
  }
  return res;
}

int main() {
  const string seq = random_seq(1e9, 'D');
  const string pattern{"GATCC"};
  // cout << "pattern :\t" << pattern << endl;
  // cout << "sequence:\t" << seq << endl;

  DFA bm{"ATCG", pattern};
  auto pa = bm.find_all(seq);

  cout << "hits     :\t" << pa.size() << endl;
  // cout << "hits    :\t";
  //  print_pattern_hits(seq, pattern, pa);
}
