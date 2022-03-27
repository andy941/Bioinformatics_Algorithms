#include "05_FindingPatterns.h"
#include "Tools.h"
#include <iostream>
#include <string>

using namespace std;

/*
 * Small bug somewhere, finds 954 of 957 total.
 * Possibly because it doesn't find overlapping patterns?
 */

int main() {
  const string seq = random_seq(1e9, 'D');
  const string pattern{"GATCCGATCC"};

  DFA dfa{"ATCG", pattern};
  // dfa.print_automata();
  {
    Timer t;
    auto pa = dfa.occurrences_pattern(seq);
    // print_pattern_hits(seq, pattern, pa);
    cout << "hits     :\t" << pa.size() << endl;
  }
}
