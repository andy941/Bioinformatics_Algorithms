#include<string>
std::string random_seq(size_t, char);

class Timer {

  std::chrono::high_resolution_clock::time_point begin =
      std::chrono::high_resolution_clock::now();

public:
  Timer() = default;
  ~Timer();
};

