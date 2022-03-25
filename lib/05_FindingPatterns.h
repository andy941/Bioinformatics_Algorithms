#include <string>
#include <unordered_map>
#include <vector>

// 02 Naive Find
std::string::const_iterator naive_find(const std::string &str,
                                       const std::string &pattern);
std::vector<std::string::const_iterator>
naive_find_all(const std::string &str, const std::string &pattern);

// 03 Boyer-Moore
class BoyerMoore {

  std::string alphabet;
  std::string pattern;
  std::unordered_map<char, int> occ;
  std::vector<int> f;
  std::vector<int> s;

public:
  void preprocess_bcr();
  void preprocess_gsr();
  std::vector<std::string::const_iterator> find_all(const std::string &text);

public:
  BoyerMoore(std::string alphabet, std::string pattern);
};

// 04 Deterministic Finite Automata
int overlap(const std::string &s1, const std::string &s2);

class DFA {
  std::string alphabet;
  std::string pattern;
  unsigned long int numstates;
  std::unordered_map<std::string, int> transition_table;

  char next_state(const int &, const char &);
  void build_transition_table(const std::string &);

public:
  DFA(std::string alphabet, std::string pattern);
  void print_automata();
  std::vector<int> apply_seq(const std::string &);
  std::vector<std::string::const_iterator>
  occurrences_pattern(const std::string &);
};
