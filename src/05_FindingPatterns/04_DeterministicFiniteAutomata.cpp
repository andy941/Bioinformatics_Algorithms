#include "05_FindingPatterns.h"
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

int main() {
  const string seq = random_seq(1e9, 'D');
  const string pattern{"GATCGATCGATCGATCCCCCGATCC"};

  DFA dfa{"ATCG", pattern};
  {
    Timer t;
    auto pa = dfa.occurrences_pattern(seq);
  }

  // cout << "hits     :\t" << pa.size() << endl;

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
