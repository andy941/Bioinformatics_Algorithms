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
 * 9. Evaluate the significance of the HSP score.
 * 10. Make two or more HSP regions into a longer alignment.
 * 12. Report every match whose expect score is lower than a threshold parameter
 * E.
 *
 * - I am not implementing the following steps:
 * 1. Remove low-complexity region or sequence repeats in the query sequence.
 * 11. Show the gapped Smith-Waterman local alignments of the query and each of
 * the matched database sequences.
 *
 * I am using a threshold score of 12 for protein kmers (3 AA).
 */

using namespace std;

int main() {
  // For testing
  //// >AT1G01060.3
  // string query{"MDTSGEELLAELLAKA"};
  // BLAST_db db{"./data/Ath_small-test_BLAST_pep.fasta",
  // "./data/BLOSUM62.csv"}; db.blast_sequence(query, 0.0001);
  // db.print_report(5);

  // >AT1G01040.1
  // There are deletions and substitutions
  string query{"IVFGQKSHPTYGIRGAVASFDVVRASGLLPVRDAEVEEDLSKGKLMMADGCMVAEDLIGK"
               "IVTAAHSGKRGFYDSICYDMSAETSFPTHHHYLGPLEYNTYADYYKQKYGVDLNCKQQ"};
  BLAST_db db{"./data/Ath_Chr1_1-1000000_Peptide.fasta", "./data/BLOSUM62.csv"};
  db.blast_sequence(query, 0.01);
  db.print_report(3);
}
