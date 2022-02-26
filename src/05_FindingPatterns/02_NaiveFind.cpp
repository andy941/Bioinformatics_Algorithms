#include "Tools.h"
#include <cstdint>
#include <iostream>
#include <iterator>
#include <vector>

using namespace std;

auto naive_find(const string &str, const string &pattern) {
  Timer t;
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

auto naive_find_all(const string &str, const string &pattern) {
  Timer t;
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
  const string seq = random_seq(1e9, 'D');
  const string pattern{"GATCC"};
  // cout << "pattern :\t" << pattern << endl;
  // cout << "sequence:\t" << seq << endl;
  const string seq2{"AGTCCGATGGAGAAATG"};

  auto p = naive_find(seq, pattern);

  if (p == seq.end()) {
    cout << "Not Found!\n";
    return 0;
  }

  cout << "first   :\t";
  if (p == seq.end())
    cout << "Not Found!" << endl;
  else
    cout << "Found!" << endl;
  // print_pattern_hits(seq, pattern, p);
  auto pa = naive_find_all(seq, pattern);
  cout << "hits     :\t" << pa.size() << endl;
  // print_pattern_hits(seq, pattern, pa);
}
