#include <iostream>
#include <string>
#include <type_traits>
#include <unordered_map>

// 03 Objective function
std::unordered_map<std::string, int> create_submat(const int &match,
                                                   const int &mismatch,
                                                   const std::string &alphabet);

std::unordered_map<std::string, int> read_submat(const std::string &file);

int score_pos(const char &c1, const char &c2,
              std::unordered_map<std::string, int> &sm, const int &gap_cost);

int score_align(const std::string &s1, const std::string &s2,
                std::unordered_map<std::string, int> &sm, const int &gap_cost);

int score_align_gapaff(const std::string &s1, const std::string &s2,
                       std::unordered_map<std::string, int> &sm,
                       const int &gap_cost, const int &gap_aff);

// 04 NeedleMan-Wunsch
int max3t(const int &a, const int &b, const int &c);

class alignment {
  std::string a;
  std::string b;

public:
  alignment() = default;
  alignment(const std::string s1, const std::string s2) : a{s1}, b{s2} {};
  void print();
  void add(const char &ac, const char &bc);
};

class needleman_Wunsch {
  std::string s1;
  std::string s2;
  int *S{nullptr};
  int *T{nullptr};
  int dim1{0};
  int dim2{0};
  std::unordered_map<std::string, int> sm;
  alignment aln;

public:
  needleman_Wunsch() = delete;
  ~needleman_Wunsch();
  needleman_Wunsch(const int &match, const int &mismatch,
                   const std::string &alphabet);
  needleman_Wunsch(const std::string &file);

  void align_sequences(const std::string &s1, const std::string &s2);
  void recover_alignment();
  void reset();
  void print();
};
