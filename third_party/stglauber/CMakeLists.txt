cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
cmake_policy(SET CMP0042 NEW)
cmake_policy(SET CMP0048 NEW)

# project title
project(stglauber VERSION 0.1.0)

set(CMAKE_BUILD_TYPE RELEASE)

# build options
option(BUILD_SHARED_LIBS "Build libstglauber.so" ON)

# compiler & system checks
include(cmake/SysCheck.cmake)

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

# exported variables
set(STGLAUBER_LIBS)
set(STGLAUBER_INCLUDE)

# finally, add the source tree
add_subdirectory(stglauber)

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

# export the project variables if not root directory
if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	set(STGLAUBER_INCLUDE_DIRS ${STGLAUBER_INCLUDE} ${ROOT_INCLUDE_DIRS} PARENT_SCOPE)
	set(STGLAUBER_LIBS ${STGLAUBER_LIBRARIES} PARENT_SCOPE)
  set(STGLAUBER_DEPENDENCY_LIBS ${STGLAUBER_DEPENDENCY_LIBS} PARENT_SCOPE)
endif(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
