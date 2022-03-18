#include "06_SequenceAlignment.h"
#include <iostream>
#include <string>
#include <utility>
#include <vector>

/*
 * Given two lists of sequences return a list of index pairs referring to the
 * the most similar sequences in both lists. Similarity is defined as sequence
 * identity.
 */

using namespace std;

int main() {
  vector<string> l1{"AATGT", "TGGG", "TAGT", "GATC", "GGGT"};
  vector<string> l2{"ATTGT", "TAGT", "GGTC", "GGG"};

  cout << "\nSequence mates  -------------------------------------" << endl;

  vector<int> v = find_mate(l1, l2);
  for (int i = 0; i < l1.size(); i++) {
    cout << l1[i] << "\t -> \t" << l2[v[i]] << "\tat index " << v[i] << endl;
  }
}
