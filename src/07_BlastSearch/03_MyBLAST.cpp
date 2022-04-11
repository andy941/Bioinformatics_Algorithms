#include "05_FindingPatterns.h"
#include "07_DatabaseSearch.h"
#include "Tools.h"
#include <iostream>
#include <string>

/*
 * Given two lists of sequences return a list of index pairs referring to the
 * the most similar sequences in both lists. Similarity is defined as sequence
 * identity. Do not consider gaps in this implementation (Chapter example).
 */

using namespace std;

int main() {

  // AT1G02750.2
  string query{"TTGTGTCACCATATCGACGAAGAACATCGTCATGAAGCTAACAATGGGATATGTCCTGTA"};

  BLAST_db db{"data/Ath_Chr1_1-1000000_CDS.fasta", 1, 0, "ATCG"};

  {
    Timer t;
    db.find_sequence(query, 31);
    cout << endl;
    cout << "\nSearch time -----------------------------------------" << endl;
    cout << endl;
  }

  cout << endl;
  cout << "\nSequence hits ---------------------------------------" << endl;
  cout << endl;

  db.print_report();
}
