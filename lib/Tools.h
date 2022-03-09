#include <chrono>
#include <iostream>
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

template <typename T>
void print_matrix(std::string &rownames, std::string &colnames, T *mat);

template <class T>
void print_matrix(std::string &rownames, std::string &colnames, T *mat) {
  int dim1 = rownames.size();
  int dim2 = colnames.size();
  std::cout << ' ';
  for (int i = 0; i < dim2; i++)
    std::cout << ' ' << colnames[i];
  std::cout << std::endl;
  for (int i = 0; i < dim1; i++) {
    std::cout << rownames[i];
    for (int j = 0; j < dim2; j++)
      std::cout << ' ' << mat[j + i * dim2];
    std::cout << std::endl;
  }
};
