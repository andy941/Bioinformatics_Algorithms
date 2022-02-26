#include "Tools.h"
#include <array>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <vector>

using namespace std;

class DotPlot {
  const string s1;
  const string s2;
  bool *mat{nullptr};
  int dim{0};

public:
  DotPlot() = delete;
  ~DotPlot() { delete[] mat; };
  DotPlot(const string, const string);

public:
  void print();
  void compare();
  void denoise();
};

DotPlot::DotPlot(const string seq1, const string seq2) : s1{seq1}, s2{seq2} {
  if (seq1.length() != seq2.length())
    cerr << "sequences must have same size" << endl;
  dim = seq1.length();
  mat = new bool[dim * dim];
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      mat[i + j * dim] = 0;
    }
  }
};

void DotPlot::print() {
  cout << ' ' << s2 << endl;
  for (int i = 0; i < dim; i++) {
    cout << s1[i];
    for (int j = 0; j < dim; j++) {
      if (mat[i + j * dim])
        cout << "*";
      else
        cout << "|";
    }
    cout << endl;
  }
}

void DotPlot::compare() {
  for (int i = 0; i < dim; i++) {
    if (s1[i] == s2[i])
      mat[i + i * dim] = 1;
    else
      mat[i + i * dim] = 0;
  }
}
void DotPlot::denoise() {}

int main() {
  const string seq1{"GATCCCAATGAGAAATG"};
  const string seq2{"AGTCCGATGGAGAAATG"};
  DotPlot dp{seq1, seq2};
  dp.print();
  cout << endl;
  dp.compare();
  dp.print();
}
