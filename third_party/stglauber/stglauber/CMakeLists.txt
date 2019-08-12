# declare our source files GLAUBER_SRCS is used to build libglauber

set(STGLAUBER_SRCS)
set(STGLAUBER_HEADERS)

set(STGLAUBER_INCLUDE)

# add all subdirectories
add_subdirectory(StFastGlauberMcMaker)
add_subdirectory(StGlauberUtilities)
add_subdirectory(StCentralityMaker)
add_subdirectory(StGlauberTree)
add_subdirectory(StGlauberAnalysisMaker)

# compile libraries
add_library(stglauber ${STGLAUBER_SRCS})
target_link_libraries(stglauber PRIVATE ${GLAUBER_DEPENDENCY_LIBS})
target_include_directories(stglauber INTERFACE $<INSTALL_INTERFACE:include>)
target_compile_options(stglauber INTERFACE "--std=c++11")
install(TARGETS stglauber EXPORT stglauberTargets DESTINATION lib)
link_directories(${CMAKE_INSTALL_PREFIX}/lib)

set(STGLAUBER_LIBS stglauber)

set(STGLAUBER_LIBS ${STGLAUBER_LIBS} PARENT_SCOPE)
set(STGLAUBER_INCLUDE ${STGLAUBER_INCLUDE} PARENT_SCOPE)
