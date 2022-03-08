#include "06_SequenceAlignment.h"
#include <iostream>
#include <unordered_map>

using namespace std;

int main() {
  unordered_map<string, int> submat = create_submat(2, -2, "ATCG");
  cout << endl;
  cout << "Substitution Matrix --------------------------------" << endl;
  cout << endl;
  for (auto &x : submat)
    cout << x.first << " : " << x.second << endl;

  cout << endl;
  cout << "Substitution Matrix BLOSUM62 -----------------------" << endl;
  cout << endl;
  unordered_map<string, int> submat_blosum = read_submat("./data/BLOSUM62.csv");
  cout << "BLOSUM62: " << submat_blosum.size() << " entries." << endl;

  const string d1{"-CAGTGCATG-ACATA"};
  const string d2{"TCAG-GC-TCTACATA"};
  int dna_score = score_align(d1, d2, submat, -3); // = -8
  cout << endl;
  cout << "Score DNA ------------------------------------------" << endl;
  cout << endl;
  cout << "seq1: " << d1 << endl;
  cout << "seq2: " << d2 << endl;
  cout << "Score = " << dna_score << endl;

  const string p1{"LGPSSGCASRIWTKSA"};
  const string p2{"TGPS-G--S-IWSKSG"};
  int prot_score = score_align(p1, p2, submat_blosum, -3);
  int prot_score_aff = score_align_gapaff(p1, p2, submat_blosum, -3, -1);
  cout << endl;
  cout << "Score Protein --------------------------------------" << endl;
  cout << endl;
  cout << "seq1: " << p1 << endl;
  cout << "seq2: " << p2 << endl;
  cout << "Score = " << prot_score << endl;
  cout << "Score Gaps Affinity = " << prot_score_aff << endl;
}
