# See: https://cmake.org/examples/

cmake_minimum_required(VERSION 3.17)
project(Algos)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)


# Chapter 05-06-07 ------------------------------------------------------------
add_subdirectory(lib)

list(
  APPEND dirs
  05_FindingPatterns
  06_SequenceAlignment
  #07_BlastSearch
)

foreach(dir ${dirs})
  file(GLOB files "${CMAKE_SOURCE_DIR}/src/${dir}/*.cpp")
  foreach(file ${files})
    get_filename_component(name ${file} NAME)
    get_filename_component(basename ${name} NAME_WE)
    add_executable(${basename}.o ${file})
    target_link_libraries(${basename}.o lib)
    install(TARGETS ${basename}.o DESTINATION ${CMAKE_SOURCE_DIR}/src/${dir}/)
  endforeach()
endforeach()

# Chapter 07 ------------------------------------------------------------------
add_executable(My_improved_BLAST.o "${CMAKE_SOURCE_DIR}/src/07_BlastSearch/My_improved_BLAST.cpp")

target_link_libraries(
  My_improved_BLAST.o
  lib_07_DatabaseSearch_improved
)
install(TARGETS My_improved_BLAST.o DESTINATION "${CMAKE_SOURCE_DIR}/src/07_BlastSearch/")

# Chapter 08 ------------------------------------------------------------------
