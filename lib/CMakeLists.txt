add_library(
  lib
  "Tools.cpp"
  "06_SequenceAlignment.cpp"
  "07_DatabaseSearch.cpp"
)
target_include_directories(lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
