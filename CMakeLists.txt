cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
cmake_policy(SET CMP0042 NEW)
cmake_policy(SET CMP0048 NEW)

# project title
project(sct VERSION 0.1.0)

set(CMAKE_BUILD_TYPE RELEASE)

# build options
option(BUILD_SHARED_LIBS "Build libsct.so" ON)
option(BUILD_BINARIES "Build executables" OFF)
option(BUILD_TEST "Build test suite binaries" OFF)
option(COMPARE_STGLAUBER "Build comparison routines with StFastGlauber. Requires GTest" OFF)

if (COMPARE_STGLAUBER AND NOT BUILD_TEST AND NOT BUILD_BINARIES)
message(STATUS "Can not run StGlauber comparisons without GTest and Boost: must enable binaries and tests")
set(COMPARE_STGLAUBER OFF)
endif (COMPARE_STGLAUBER AND NOT BUILD_TEST)

# compiler & system checks
include(cmake/SysCheck.cmake)

# add module path to cmake
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

# helper functions
include(cmake/Utils.cmake)

# dependencies
include(cmake/Dependencies.cmake)

# build options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2 -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")

# forcing proper include lookup paths - must take configured files in build dir
# before empty files in the source dir
include_directories(BEFORE ${PROJECT_SOURCE_DIR})
include_directories(BEFORE ${PROJECT_BINARY_DIR})

# finally, add the source tree
add_subdirectory(sct)

# add an uninstall option, since cmake does not provide one by default
if(NOT TARGET uninstall)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in
                 ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
                 IMMEDIATE
                 @ONLY)

  add_custom_target(uninstall
                    COMMAND ${CMAKE_COMMAND} -P
                            ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
