#include "05_FindingPatterns.h"
#include "Tools.h"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <unordered_map>
#include <vector>

using namespace std;

int main() {
  const string seq = random_seq(1e9, 'D');
  const string pattern{"GATCC"};
  // cout << "pattern :\t" << pattern << endl;
  // cout << "sequence:\t" << seq << endl;

  BoyerMoore bm{"ATCG", pattern};
  auto pa = bm.find_all(seq);

  cout << "hits     :\t" << pa.size() << endl;
  // cout << "hits    :\t";
  //  print_pattern_hits(seq, pattern, pa);
}
