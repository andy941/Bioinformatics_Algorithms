#include "05_FindingPatterns.h"
#include "Tools.h"

using namespace std;

int main() {
  const string seq = random_seq(1e9, 'D');
  const string pattern{"GATCGATCGATCGATCCCCCGATCC"};
  // cout << "pattern :\t" << pattern << endl;
  // cout << "sequence:\t" << seq << endl;

  BoyerMoore bm{"ATCG", pattern};
  {
    Timer t;
    auto pa = bm.find_all(seq);
  }

  // cout << "hits     :\t" << pa.size() << endl;
  //  cout << "hits    :\t";
  //   print_pattern_hits(seq, pattern, pa);
}
