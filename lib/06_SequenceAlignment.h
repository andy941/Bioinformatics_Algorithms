#include <string>
#include <unordered_map>

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
