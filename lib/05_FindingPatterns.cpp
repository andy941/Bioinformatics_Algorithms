#include "05_FindingPatterns.h"
#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

// 02 Naive Find
std::string::const_iterator naive_find(const std::string &str,
                                       const std::string &pattern) {
  if (pattern.begin() == pattern.end() || str.begin() == str.end())
    return str.end();

  auto p = pattern.begin();
  for (auto x = str.begin(); x != str.end(); advance(x, 1)) {
    auto tmp = x;
    while (*tmp == *p) {
      if (p + 1 == pattern.end())
        return x;
      advance(p, 1);
      advance(tmp, 1);
    }
    p = pattern.begin();
  }
  return str.end();
}

std::vector<std::string::const_iterator>
naive_find_all(const std::string &str, const std::string &pattern) {
  std::vector<std::string::const_iterator> res;
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

// 03 Boyer-Moore
BoyerMoore::BoyerMoore(std::string alphabet, std::string pattern)
    : alphabet{alphabet}, pattern{pattern} {
  preprocess_bcr();
  preprocess_gsr();
};

void BoyerMoore::preprocess_bcr() {
  for (char s : alphabet)
    occ[s] = -1;
  for (int i = 0; i < pattern.length(); i++)
    occ[pattern[i]] = i;
}

void BoyerMoore::preprocess_gsr() {
  f = std::vector<int>(pattern.length() + 1, 0);
  s = std::vector<int>(pattern.length() + 1, 0);
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

std::vector<std::string::const_iterator>
BoyerMoore::find_all(const std::string &text) {

  std::vector<std::string::const_iterator> res;
  res.reserve(text.size() / pattern.size());
  int i = 0;
  int j = pattern.length() - 1;
  char c;

  while (i <= text.length() - pattern.length()) {

    j = pattern.length() - 1;
    while (j >= 0 && pattern[j] == text[j + i]) {
      --j;
    }
    if (j < 0) {
      res.emplace_back(text.begin() + i);
      i += s[0];
    } else {
      c = text[j + i];
      i += fmax(s[j + 1], j - occ[c]);
    }
  }
  return res;
}

// 04 Deterministic Finite Automata
int overlap(const std::string &s1, const std::string &s2) {
  int maxov = std::min(s1.length(), s2.length());
  for (int i = maxov; i > 0; i--) {
    std::string ss1{s1.substr(maxov - i)};
    std::string ss2{s2.substr(0, i)};
    if (strcmp(ss1.c_str(), ss2.c_str()) == 0) {
      return i;
    }
  }
  return 0;
}

void DFA::build_transition_table(const std::string &pattern) {
  for (int q = 0; q < numstates; q++) {
    for (auto a : alphabet) {
      std::string prefix{pattern.substr(0, q) + a};
      std::string key{std::to_string(q) + a};
      transition_table[key] = overlap(prefix, pattern);
    }
  }
}

DFA::DFA(std::string alphabet, std::string pattern)
    : alphabet{alphabet}, pattern{pattern} {
  numstates = pattern.length() + 1;
  build_transition_table(pattern);
};

void DFA::print_automata() {
  std::cout << "States   : " << numstates << std::endl;
  std::cout << "Alphabet : " << alphabet << std::endl;
  std::cout << "Tr_table : " << std::endl;
  for (auto x : transition_table) {
    std::cout << x.first << " -> " << x.second << std::endl;
  }
};

char DFA::next_state(const int &q, const char &c) {
  return transition_table[std::to_string(q) + c];
};

std::vector<int> DFA::apply_seq(const std::string &seq) {
  int q = 0;
  std::vector<int> res{q};
  for (char c : seq) {
    q = next_state(q, c);
    res.push_back(q);
  };
  return res;
}

std::vector<std::string::const_iterator>
DFA::occurrences_pattern(const std::string &text) {
  int q = 0;
  std::vector<std::string::const_iterator> res;
  for (int i = 0; i < text.length(); i++) {
    q = next_state(q, text[i]);
    if (q == numstates - 1)
      res.push_back(text.begin() + i - numstates + 2);
  }
  return res;
};
