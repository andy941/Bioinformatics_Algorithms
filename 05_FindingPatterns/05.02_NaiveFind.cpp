#include "../lib/Tools.h"
#include <iostream>
using namespace std;

int main() {
  string seq = random_seq(100, 'D');
  Timer t;
  cout << seq << endl;
}
