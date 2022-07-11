# Tools -----------------------------------------------------------------------
add_library(lib_Tools "Tools.cpp")
target_include_directories(lib_Tools PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

# Chapter 05 ------------------------------------------------------------------
add_library(lib_05_FindingPatterns "05_FindingPatterns.cpp")
target_link_libraries(lib_05_FindingPatterns lib_Tools)
target_include_directories(lib_05_FindingPatterns PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

# Chapter 06 ------------------------------------------------------------------
add_library(lib_06_SequenceAlignment "06_SequenceAlignment.cpp")
target_include_directories(lib_06_SequenceAlignment PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

# Chapter 07 ------------------------------------------------------------------
add_library(lib_07_DatabaseSearch "07_DatabaseSearch.cpp")
add_library(lib_07_DatabaseSearch_improved "07_DatabaseSearch_improved.cpp")

target_link_libraries(lib_07_DatabaseSearch lib_06_SequenceAlignment lib_Tools)
target_link_libraries(lib_07_DatabaseSearch_improved lib_06_SequenceAlignment)

target_include_directories(
  lib_07_DatabaseSearch PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}"
)
target_include_directories(
  lib_07_DatabaseSearch_improved PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_SOURCE_DIR}/ext/Eigen/"
)

# Chapter 08 ------------------------------------------------------------------
add_library(lib_08_MultipleSequenceAlignment "08_MultipleSequenceAlignment.cpp")
target_link_libraries(lib_08_MultipleSequenceAlignment lib_06_SequenceAlignment)
target_include_directories(
  lib_08_MultipleSequenceAlignment PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_SOURCE_DIR}/ext/Eigen/"
)
