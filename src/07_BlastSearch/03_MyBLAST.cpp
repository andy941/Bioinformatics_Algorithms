#include "05_FindingPatterns.h"
#include "07_DatabaseSearch.h"
#include <iostream>
#include <string>

/*
 * Given two lists of sequences return a list of index pairs referring to the
 * the most similar sequences in both lists. Similarity is defined as sequence
 * identity. Do not consider gaps in this implementation (Chapter example).
 */

using namespace std;

int main() {

  // AT1G01220.4
  string query{
      "GTATTGTGGTCTCCATGACAATCCAAAGAACTCAATTCATAAAGATGGAACTTTTTGCGGTAAACCCTTGGA"
      "GAAGGTATTGTTTGATCTTGGCATTGAGGAAAGCGACCTCTGGAGCTCGTATGTTGCACAAGATAGATGTTT"
      "GTGGAATGCAAAACTGTTCCCGATTCTTACGTATAGTGAAATGCTGAAGTTAGCGTCGTGGTTGATGGGTTT"
      "AGATGATAGTAGAAACAAGGAGAAGATTAAGTTGTGGAGAAGCTCACAACGTGTAAGCTTAGAAGAGTTGCA"
      "TGGATCAATCAACTTTCCTGAGATGTGCAATGGTTCCAGCAATCATCAAGCTGATCTTGCGGGTGGAATCGC"
      "TAAAGCATGTATGAACTATGGTATGCTTGGGCGTAATTTGTCTCAGCTGTGCCATGAGATTTTACAGAAAGA"
      "GTCATTAGGATTGGAAATATGCAAGAATTTTCTGGATCAATGTCCCAAATTTCAGGAGCAGAACTCCAAAAT"
      "TCTTCCAAAGAGTCGAGCATA"};

  cout << "\nSequence hits ---------------------------------------" << endl;

  BLAST_db db{"data/Ath_Chr1_1-1000000_CDS.fasta", 1, 0, "ATCG"};
  db.find_sequence(query, 5);
}
