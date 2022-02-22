#include "Tools.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

int overlap(const string &s1, const string &s2) {
  int maxov = min(s1.length(), s2.length());
  for (int i = maxov; i > 0; i--) {
    string ss1{s1.substr(maxov - i)};
    string ss2{s2.substr(0, i)};
    if (strcmp(ss1.c_str(), ss2.c_str()) == 0) {
      return i;
    }
  }
  return 0;
}

class DFA {
  string alphabet;
  string pattern;
  unsigned long int numstates;
  unordered_map<string, int> transition_table;

  char next_state(const int &, const char &);
  void build_transition_table(const string &);

public:
  DFA(string alphabet, string pattern);
  void print_automata();
  vector<int> apply_seq(const string &);
  auto occurrences_pattern(const string &);
};

void DFA::build_transition_table(const string &pattern) {
  for (int q = 0; q < numstates; q++) {
    for (auto a : alphabet) {
      string prefix{pattern.substr(0, q) + a};
      string key{to_string(q) + a};
      transition_table[key] = overlap(prefix, pattern);
    }
  }
}

DFA::DFA(string alphabet, string pattern)
    : alphabet{alphabet}, pattern{pattern} {
  numstates = pattern.length() + 1;
  build_transition_table(pattern);
};

void DFA::print_automata() {
  cout << "States   : " << numstates << endl;
  cout << "Alphabet : " << alphabet << endl;
  cout << "Tr_table : " << endl;
  for (auto x : transition_table) {
    cout << x.first << " -> " << x.second << endl;
  }
};

char DFA::next_state(const int &q, const char &c) {
  return transition_table[to_string(q) + c];
};

vector<int> DFA::apply_seq(const string &seq) {
  int q = 0;
  vector<int> res{q};
  for (char c : seq) {
    q = next_state(q, c);
    res.push_back(q);
  };
  return res;
}

auto DFA::occurrences_pattern(const string &text) {
  Timer t;
  int q = 0;
  vector<string::const_iterator> res;
  for (int i = 0; i < text.length(); i++) {
    q = next_state(q, text[i]);
    if (q == numstates - 1)
      res.push_back(text.begin() + i - numstates + 2);
  }
  return res;
};

int main() {
  const string seq = random_seq(1e9, 'D');
  const string pattern{"GATCC"};

  DFA dfa{"ATCG", pattern};
  auto pa = dfa.occurrences_pattern(seq);

  cout << "hits     :\t" << pa.size() << endl;

  dfa.print_automata();
  // cout << "pattern :\t " << pattern << endl;
  // cout << "sequence:\t " << seq << endl;
  // auto ps = dfa.apply_seq(seq);
  // cout << "states  :\t";
  // for (auto x : ps)
  //   cout << x;
  // cout << endl;
  // cout << "hits    :\t ";
  // print_pattern_hits(seq, pattern, pa);
}
