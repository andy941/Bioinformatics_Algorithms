#include "05_FindingPatterns.h"
#include "Tools.h"
#include <cstdint>
#include <iostream>
#include <iterator>
#include <vector>

using namespace std;

int main() {
  const string seq = random_seq(1e9, 'D');
  const string pattern{"GATCGATCGATCGATCCCCCGATCC"};
  // cout << "pattern :\t" << pattern << endl;
  // cout << "sequence:\t" << seq << endl;
  // const string seq2{"AGTCCGATGGAGAAATG"};

  // auto p = naive_find(seq, pattern);

  // if (p == seq.end()) {
  //   cout << "Not Found!\n";
  //   return 0;
  // }

  // cout << "first   :\t";
  // if (p == seq.end())
  //   cout << "Not Found!" << endl;
  // else
  //   cout << "Found!" << endl;
  //  print_pattern_hits(seq, pattern, p);
  {
    Timer t;
    auto pa = naive_find_all(seq, pattern);
  }
  // cout << "hits     :\t" << pa.size() << endl;
  //  print_pattern_hits(seq, pattern, pa);
}
