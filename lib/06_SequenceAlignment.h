#include <iostream>
#include <string>
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
template <class T, size_t S> int max_arr(const std::array<T, S> arr) {
  int max = arr[0];
  int r = 0;
  for (int i = 1; i < arr.size(); i++)
    if (arr[i] > max)
      r = i;
  return r;
};

struct alignment {
  std::string a;
  std::string b;

  alignment() = default;
  alignment(const std::string s1, const std::string s2) : a{s1}, b{s2} {};
  void print();
  void add(const char &ac, const char &bc);
  void flip();
  double identity();
  std::string identical_subseq();
};

class needleman_Wunsch {
protected:
  std::string s1;
  std::string s2;
  int *S{nullptr};
  int *T{nullptr};
  int dim1{0};
  int dim2{0};
  std::unordered_map<std::string, int> sm;
  int best_score{0};

public:
  alignment aln;
  needleman_Wunsch() = delete;
  ~needleman_Wunsch();
  needleman_Wunsch(const int &match, const int &mismatch,
                   const std::string &alphabet);
  needleman_Wunsch(const std::string &file);

  virtual void align_sequences(const std::string &seq1, const std::string &seq2,
                               const int &gap_cost);
  virtual void trace_back();
  void reset();
  void print();
  int get_score() { return best_score; };
};

// 05 Smith-Waterman
class smith_Waterman : public needleman_Wunsch {

public:
  using needleman_Wunsch::needleman_Wunsch;
  void align_sequences(const std::string &seq1, const std::string &seq2,
                       const int &gap_cost) override;
  void trace_back() override;
};
