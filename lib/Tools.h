#include <chrono>
#include <string>
#include <vector>

// Define some useful ANSI colors
#define RED "\033[31m"
#define RESET "\033[0m"

class Timer {

  std::chrono::high_resolution_clock::time_point begin =
      std::chrono::high_resolution_clock::now();

public:
  Timer() = default;
  ~Timer();
};

std::string random_seq(size_t, char);

void print_pattern_hits(const std::string &, const std::string &,
                        std::string::const_iterator &);
void print_pattern_hits(const std::string &, const std::string &,
                        std::vector<std::string::const_iterator> &);
