cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# project name and enable C++ support
project(wilson_line_reader CXX C)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE Release)


find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})


add_subdirectory(src)


set_target_properties(wline PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
	set_target_properties(test PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
