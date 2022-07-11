# See: https://cmake.org/examples/

cmake_minimum_required(VERSION 3.17)
project(Algos)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)


# Build libraries -------------------------------------------------------------
add_subdirectory(lib)

## Chapter 05 ------------------------------------------------------------------
#add_executable(02_NaiveFind.o "${CMAKE_SOURCE_DIR}/src/05_FindingPatterns/02_NaiveFind.cpp")
#add_executable(04_DeterministicFiniteAutomata.o "${CMAKE_SOURCE_DIR}/src/05_FindingPatterns/04_DeterministicFiniteAutomata.cpp")
#add_executable(03_BoyerMoore.o "${CMAKE_SOURCE_DIR}/src/05_FindingPatterns/03_BoyerMoore.cpp")
#
#target_link_libraries(02_NaiveFind.o lib_05_FindingPatterns)
#target_link_libraries(04_DeterministicFiniteAutomata.o lib_05_FindingPatterns)
#target_link_libraries(03_BoyerMoore.o lib_05_FindingPatterns)
#
#install(
#  TARGETS
#  02_NaiveFind.o
#  04_DeterministicFiniteAutomata.o
#  03_BoyerMoore.o
#  DESTINATION "${CMAKE_SOURCE_DIR}/src/05_FindingPatterns/"
#)
#
## Chapter 06 ------------------------------------------------------------------
#add_executable(ex05_SW-HandleTies.o "${CMAKE_SOURCE_DIR}/src/06_SequenceAlignment/ex05_SW-HandleTies.cpp")
#add_executable(ex04_NW-HandleTies.o "${CMAKE_SOURCE_DIR}/src/06_SequenceAlignment/ex04_NW-HandleTies.cpp")
#add_executable(06_SpecialCases.o "${CMAKE_SOURCE_DIR}/src/06_SequenceAlignment/06_SpecialCases.cpp")
#add_executable(ex06_FindMate.o "${CMAKE_SOURCE_DIR}/src/06_SequenceAlignment/ex06_FindMate.cpp")
#add_executable(05_Smith-Waterman.o "${CMAKE_SOURCE_DIR}/src/06_SequenceAlignment/05_Smith-Waterman.cpp")
#add_executable(04_NeedleMan-Wunsch.o "${CMAKE_SOURCE_DIR}/src/06_SequenceAlignment/04_NeedleMan-Wunsch.cpp")
#add_executable(03_ObjectiveFunction.o "${CMAKE_SOURCE_DIR}/src/06_SequenceAlignment/03_ObjectiveFunction.cpp")
#add_executable(02_Dotplot.o "${CMAKE_SOURCE_DIR}/src/06_SequenceAlignment/02_Dotplot.cpp")
#add_executable(ex03_DotplotMaxDiagonal.o "${CMAKE_SOURCE_DIR}/src/06_SequenceAlignment/ex03_DotplotMaxDiagonal.cpp")
#
#target_link_libraries(ex05_SW-HandleTies.o lib_06_SequenceAlignment)
#target_link_libraries(ex04_NW-HandleTies.o lib_06_SequenceAlignment)
#target_link_libraries(06_SpecialCases.o lib_06_SequenceAlignment)
#target_link_libraries(ex06_FindMate.o lib_06_SequenceAlignment)
#target_link_libraries(05_Smith-Waterman.o lib_06_SequenceAlignment)
#target_link_libraries(04_NeedleMan-Wunsch.o lib_06_SequenceAlignment)
#target_link_libraries(03_ObjectiveFunction.o lib_06_SequenceAlignment)
#target_link_libraries(02_Dotplot.o lib_06_SequenceAlignment)
#target_link_libraries(ex03_DotplotMaxDiagonal.o lib_06_SequenceAlignment)
#
#install(
#  TARGETS
#  ex05_SW-HandleTies.o
#  ex04_NW-HandleTies.o
#  06_SpecialCases.o
#  ex06_FindMate.o
#  05_Smith-Waterman.o
#  04_NeedleMan-Wunsch.o
#  03_ObjectiveFunction.o
#  02_Dotplot.o
#  ex03_DotplotMaxDiagonal.o
#  DESTINATION "${CMAKE_SOURCE_DIR}/src/06_SequenceAlignment/"
#)
#
## Chapter 07 ------------------------------------------------------------------
#add_executable(03_MyBLAST.o "${CMAKE_SOURCE_DIR}/src/07_BlastSearch/03_MyBLAST.cpp")
#add_executable(My_improved_BLAST.o "${CMAKE_SOURCE_DIR}/src/07_BlastSearch/My_improved_BLAST.cpp")
#
#target_link_libraries(03_MyBLAST.o lib_07_DatabaseSearch)
#target_link_libraries(My_improved_BLAST.o lib_07_DatabaseSearch_improved)
#
#install(
#  TARGETS
#  03_MyBLAST.o
#  My_improved_BLAST.o
#  DESTINATION
#  "${CMAKE_SOURCE_DIR}/src/07_BlastSearch/"
#)

# Chapter 08 ------------------------------------------------------------------
add_executable(MSA.o "${CMAKE_SOURCE_DIR}/src/08_MultipleSequenceAlignment/MSA.cpp")
target_link_libraries(MSA.o lib_08_MultipleSequenceAlignment)
install(TARGETS MSA.o DESTINATION "${CMAKE_SOURCE_DIR}/src/08_MultipleSequenceAlignment/")
