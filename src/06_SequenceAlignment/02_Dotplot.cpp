#include "06_SequenceAlignment.h"
#include <array>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <vector>

using namespace std;

int main() {
  const string seq1{"GATCCCAATGAGAAATGATGGGGATGATGATGA"};
  const string seq2{"AGTCCGATGGAGAAATGATGAAGATGATGGTATGATGATGTGAGATGATGAAA"};
  DotPlot dp{seq1, seq2};
  cout << "Seq 1 = " << seq1 << endl;
  cout << "Seq 2 = " << seq2 << endl;
  dp.print();
  cout << "Basic comparison ----------------------------" << endl;
  dp.compare();
  dp.print();
  cout << "Window denoise ------------------------------" << endl;
  dp.reset();
  dp.denoise(6, 5);
  dp.print();
}
