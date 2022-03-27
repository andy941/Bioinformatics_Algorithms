#include "05_FindingPatterns.h"
#include "Tools.h"
#include <cstdint>
#include <iostream>
#include <iterator>
#include <vector>

using namespace std;

int main() {
  const string seq = random_seq(1e9, 'D');
  const string pattern{"GATCCGATCC"};

  {
    Timer t;
    auto pa = naive_find_all(seq, pattern);
    // print_pattern_hits(seq, pattern, pa);
    cout << "hits     :\t" << pa.size() << endl;
  }
}
