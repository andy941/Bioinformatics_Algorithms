#include "06_SequenceAlignment.h"
#include <array>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <vector>

using namespace std;

int main() {
  const string seq1{"TTGAGAAATGGTGA"};
  const string seq2{"TGGAGAAATGGTATGATA"};
  DotPlot dp{seq1, seq2};
  cout << "Seq 1 = " << seq1 << endl;
  cout << "Seq 2 = " << seq2 << endl;
  dp.print();
  cout << "\nBasic comparison ----------------------------" << endl;
  dp.compare();
  dp.print();
  cout << "\nWindow denoise ------------------------------" << endl;
  dp.reset();
  dp.denoise(4, 3);
  dp.print();
}
