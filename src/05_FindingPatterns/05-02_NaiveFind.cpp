#include "Tools.h"
#include <cstdint>
#include <iostream>
#include <iterator>
#include <vector>

using namespace std;

auto naive_find(const string &str, const string &pattern) {
  if (pattern == "" || str == "")
    return str.end();

  auto p = pattern.begin();
  for (auto x = str.begin(); x != str.end(); advance(x, 1)) {
    auto tmp = x;
    while (*tmp == *p) {
      if (p + 1 == pattern.end())
        return x;
      advance(p, 1);
      advance(tmp, 1);
    }
    p = pattern.begin();
  }
  return str.end();
}

int main() {
  const string seq = random_seq(100, 'D');
  const string pattern{"TGG"};
  cout << "pattern: " << pattern << endl;
  cout << "sequence: " << seq << endl;

  auto p = naive_find(seq, pattern);
  if (p == seq.end()) {
    cout << "Not Found!\n";
    return 0;
  }
  cout << "\n--\nsequence: ";
  cout << string(seq.begin(), p);
  cout << RED << pattern << RESET;
  cout << string(p + pattern.size(), seq.end()) << endl;
  cout << "Found!\n";
  Timer t;
}
