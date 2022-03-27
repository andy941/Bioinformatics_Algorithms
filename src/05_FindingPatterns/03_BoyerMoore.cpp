#include "05_FindingPatterns.h"
#include "Tools.h"

using namespace std;

int main() {
  const string seq = random_seq(1e9, 'D');
  const string pattern{"GATCCGATCC"};

  // const string seq = random_seq(1e4, 'D');
  // const string pattern{"CGAA"};

  // cout << "pattern :\t" << pattern << endl;
  // cout << "sequence:\t" << seq << endl;
  // cout << "hits     :\t" << pa.size() << endl;
  //  cout << "hits    :\t";
  //
  BoyerMoore bm{"ATCG", pattern};
  {
    Timer t;
    auto pa = bm.find_all(seq);
    // print_pattern_hits(seq, pattern, pa);
    cout << "hits     :\t" << pa.size() << endl;
  }
}
