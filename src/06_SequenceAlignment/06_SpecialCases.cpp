#include "06_SequenceAlignment.h"
#include <iostream>
#include <unordered_map>

using namespace std;
int main() {

  cout << "\nSequences for analysis -------------------------------------"
       << endl;
  string seq1 = "AATGTCGCAAGGGCTAGTATCGACC";
  string seq2 = "AATGTGTTAAGGGCTTGTAGCATGA";

  cout << endl;
  cout << "seq1: " << seq1 << endl;
  cout << "seq2: " << seq2 << endl;

  smith_Waterman sW{1, 0, "ATCG"};

  cout << "\nSequence Identity ------------------------------------------"
       << endl;
  needleman_Wunsch nWi{1, 0, "ATCG"};
  nWi.align_sequences(seq1, seq2, 0);
  nWi.trace_back();
  nWi.aln.print();
  cout << ">>> Identity = " << nWi.aln.identity() * 100 << '%' << endl;

  cout << "\nSequence Distance ------------------------------------------"
       << endl;
  needleman_Wunsch nWd{0, -1, "ATCG"};
  nWd.align_sequences(seq1, seq2, -1);
  nWd.trace_back();
  nWd.aln.print();
  cout << ">>> Distance = " << -nWd.get_score() << endl;

  cout << "\nLongest common subsequence ---------------------------------"
       << endl;
  needleman_Wunsch nWs{1, 0, "ATCG"};
  nWs.align_sequences(seq1, seq2, 0);
  nWs.trace_back();
  nWs.aln.print();
  cout << ">>> Identical subseq = " << nWs.aln.identical_subseq() << endl;

  cout << "\nLongest identical contiguous -------------------------------"
       << endl;
  int pen = -(max(seq1.size(), seq2.size()));
  smith_Waterman sWs{1, pen, "ATCG"};
  sWs.align_sequences(seq1, seq2, pen);
  sWs.trace_back();
  sWs.aln.print();
  cout << ">>> Identical contiguous = " << sWs.aln.identical_subseq() << endl;
}
