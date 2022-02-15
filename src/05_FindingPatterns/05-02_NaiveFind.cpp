#include "Tools.h"
#include <cstdint>
#include <iostream>
#include <iterator>
#include <vector>

using namespace std;

string::iterator naive_find(string &str, string &pattern) {
  if (pattern == "")
    return str.end();
  long long int count = 0;
  char c = pattern[0];
  for (string::iterator x = str.begin(); x != str.end(); advance(x, 1)) {
    if (*x == c) {
      count += 1;
    } else
      count = 0;
    if (count == pattern.size())
      return x;
  }
  return str.end();
}

int main() {
  string seq = random_seq(100, 'D');
  string pattern{"TGG"};
  cout << "pattern: " << pattern << endl;
  cout << "sequence: " << seq << endl;

  auto p = naive_find(seq, pattern);
  if (p == seq.end()) {
    cout << "Not Found!\n";
    return 0;
  }
  cout << "\n--\nsequence: ";
  cout << string(seq.begin(), p - 1);
  cout << RED << pattern << RESET;
  cout << string(p + pattern.size(), seq.end()) << endl;
  cout << "Found!\n";
  Timer t;
}
