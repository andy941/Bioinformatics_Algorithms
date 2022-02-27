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
  int dim1{0};
  int dim2{0};

public:
  DotPlot() = delete;
  ~DotPlot() { delete[] mat; };
  DotPlot(const string, const string);

public:
  void print();
  void reset();
  void compare();
  void denoise(unsigned int window, unsigned int stringency);
};

DotPlot::DotPlot(const string seq1, const string seq2) : s1{seq1}, s2{seq2} {
  dim1 = seq1.length();
  dim2 = seq2.length();
  mat = new bool[dim1 * dim2];
  for (int i = 0; i < dim1; i++) {
    for (int j = 0; j < dim2; j++) {
      mat[i + j * dim1] = 0;
    }
  }
};

void DotPlot::reset() {
  for (int i = 0; i < dim1; i++) {
    for (int j = 0; j < dim2; j++) {
      mat[i + j * dim1] = 0;
    }
  }
}

void DotPlot::print() {
  cout << ' ';
  for (int i = 0; i < dim2; i++)
    cout << ' ' << s2[i];
  cout << endl;
  for (int i = 0; i < dim1; i++) {
    cout << s1[i];
    for (int j = 0; j < dim2; j++) {
      if (mat[i + j * dim1])
        cout << " *";
      else
        cout << " |";
    }
    cout << endl;
  }
}

void DotPlot::compare() {
  for (int i = 0; i < dim1; i++) {
    for (int j = 0; j < dim2; j++) {
      if (s1[i] == s2[j])
        mat[i + j * dim1] = 1;
      else
        mat[i + j * dim1] = 0;
    }
  }
}
inline unsigned int match(string::const_iterator b1, string::const_iterator e1,
                          string::const_iterator b2,
                          string::const_iterator e2) {
  unsigned int res{0};
  while (b1 != e1 && b2 != e2) {
    if (*b1++ == *b2++)
      res++;
  }
  return res;
}
void DotPlot::denoise(unsigned int window, unsigned int stringency) {
  unsigned int matches{0};
  unsigned int start{window / 2};
  for (unsigned int i = start; i < dim1 - start; i++) {
    for (unsigned int j = start; j < dim2 - start; j++) {
      matches = match(s1.begin() + i - start, s1.begin() + i + start + 1,
                      s2.begin() + j - start, s2.begin() + j + start + 1);
      if (matches >= stringency)
        mat[i + j * dim1] = 1;
    }
  }
}
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
