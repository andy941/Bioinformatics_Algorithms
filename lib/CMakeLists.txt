# Chapter 05 - 06 -07
add_library(
  lib
  "Tools.cpp"
  "05_FindingPatterns.cpp"
  "06_SequenceAlignment.cpp"
  "07_DatabaseSearch.cpp"
)
target_include_directories(
  lib
  PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_SOURCE_DIR}/ext/Eigen/"
)

# This is the correct, non lazy way of doing things
# TODO: port the rest of the code

# Chapter 07 ------------------------------------------------------------------
add_library(
  lib_06_SequenceAlignment
  "06_SequenceAlignment.cpp"
)
target_include_directories(
  lib_06_SequenceAlignment
  PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

add_library(
  lib_07_DatabaseSearch_improved
  "07_DatabaseSearch_improved.cpp"
)
target_link_libraries(
  lib_07_DatabaseSearch_improved
  lib_06_SequenceAlignment
)
target_include_directories(
  lib_07_DatabaseSearch_improved
  PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_SOURCE_DIR}/ext/Eigen/"
)

# Chapter 08 ------------------------------------------------------------------
