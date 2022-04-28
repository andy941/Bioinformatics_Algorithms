#include "07_DatabaseSearch_improved.h"
//#include "Tools.h"
#include <iostream>
#include <string>

/*
 * - Overview of the implemented algorithm (see:
 * https://en.wikipedia.org/wiki/BLAST_(biotechnology)#Process) .
 * 2. Make a k-letter word list of the query sequence.
 * 3. List the possible matching words.
 * 4. Organize the remaining high-scoring words into an efficient search tree.
 * 5. Repeat 3-4 for each k-letter word in the query sequence
 * 6. Scan the database sequences for exact matches with the remaining
 * high-scoring words.
 * 7. Extend the exact matches to high-scoring segment pair (HSP).
 * 8. List all of the HSPs in the database whose score is high enough to be
 * considered.
 *
 * - I am not implementing the E-score calculation.
 * 1. Remove low-complexity region or sequence repeats in the query sequence.
 * 9. Evaluate the significance of the HSP score.
 * 10. Make two or more HSP regions into a longer alignment.
 * 11. Show the gapped Smith-Waterman local alignments of the query and each of
 * the matched database sequences.
 * 12. Report every match whose expect score is lower than a threshold parameter
 * E.
 *
 */

using namespace std;

int main() {
  // AT1G02750.2
  string query{"GAAGAACATCGTCATGAAGCTAACAATGGGATATGTCCTGTA"};

  BLAST_db db{"./data/Ath_Chr1_1-1000000_Peptide.fasta", "./data/BLOSUM62.csv"};
}
