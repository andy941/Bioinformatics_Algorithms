#include "Tools.h"
#include <chrono>
#include <iostream>
#include <random>
#include <stdio.h>
#include <string>

using namespace std;

Timer::~Timer() {
  auto end = std::chrono::high_resolution_clock::now();
  auto elapsed =
      std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
  printf("Time: \t%.10f seconds.\n", elapsed.count() * 1e-9);
}

string random_seq(size_t n, char kind) {
  string result;
  string space;
  unsigned short int sz;
  result.reserve(n);
  switch (kind) {
  case 'R':
    space = "AUGC";
    break;
  case 'P':
    space = "ARNDCEQGHILKMFPSTWYV";
    break;
  case 'D':
    space = "ATGC";
    break;
  default:
    cerr << kind << " <- Not supported" << endl;
    return "";
  }
  sz = space.length();

  srand(33);
  for (auto i = 0; i < n; i++)
    result.push_back(space[rand() % sz]);
  return result;
}

void print_pattern_hits(const string &seq, const string &pattern,
                        vector<string::const_iterator> &hits) {
  auto beg = seq.begin();
  auto end = seq.end();
  for (auto x : hits) {
    cout << string(beg, x);
    cout << RED << pattern << RESET;
    beg = x;
    advance(beg, pattern.size());
  }
  cout << string(beg, end) << endl;
}

void print_pattern_hits(const string &seq, const string &pattern,
                        string::const_iterator &hit) {
  cout << string(seq.begin(), hit);
  cout << RED << pattern << RESET;
  cout << string(hit + pattern.size(), seq.end()) << endl;
}
