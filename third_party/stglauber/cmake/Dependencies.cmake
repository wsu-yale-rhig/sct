# packages are either found or built, depending on if they are statically or
# dynamically linked

set(STGLAUBER_DEPENDENCY_LIBS "")

# ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED
             COMPONENTS MathCore
                        RIO
                        Hist
                        Tree
                        Net)
                    list(APPEND GLAUBER_DEPENDENCY_LIBS ${ROOT_LIBRARIES})
include(${ROOT_USE_FILE})
message(STATUS "Found ROOT")

