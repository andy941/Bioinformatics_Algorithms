#pragma once

#include "06_SequenceAlignment.h"
#include <string>
#include <vector>

class multiple_alignment {
  std::vector<std::string> alignments{};
  int score{0};

public:
  multiple_alignment() = default;
  int get_score();
  std::string get_consensus();
};

class MSA {
  multiple_alignment m_aln{};
  needleman_Wunsch nw;

public:
  MSA();
};
