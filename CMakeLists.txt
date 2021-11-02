cmake_minimum_required(VERSION 3.17)
project(Algorithms)
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

# easy way to add object files, better qould be to compile into a static add_library
# but I can't easily achieve that. I'll look into it when I have some time.
# It works for now.
file(GLOB libs "${CMAKE_SOURCE_DIR}/lib/*.cpp")
add_library( mysharedlib OBJECT ${libs} )

list(APPEND dirs 05_FindingPatterns)

foreach(dir ${dirs})
	file(GLOB files "${CMAKE_SOURCE_DIR}/src/${dir}/*.cpp")
	foreach(file ${files})
		get_filename_component(name ${file} NAME)
		get_filename_component(basename ${name} NAME_WE)
		add_executable( ${basename}.o src/${dir}/${name} $<TARGET_OBJECTS:mysharedlib>)
		target_include_directories(${basename}.o PUBLIC lib)
		install( TARGETS ${basename}.o DESTINATION ${CMAKE_SOURCE_DIR}/src/${dir}/)
	endforeach()
endforeach()
