#include <chrono>
#include <string>
#include <vector>

// Define some useful ANSI colors
#define RED	"\033[31m"
#define RESET "\033[0m"

std::string random_seq(size_t, char);

class Timer {

  std::chrono::high_resolution_clock::time_point begin =
      std::chrono::high_resolution_clock::now();

public:
  Timer() = default;
  ~Timer();
};
