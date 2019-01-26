## declare our source files
## SCT_SRCS is used to build libsct
## SCT_TEST_SRCS contains all test sources

set(SCT_SRCS)
set(SCT_HEADERS)
set(SCT_TEST_SRCS)
set(SCT_BENCH_SRCS)
set(SCT_BINARY_SRCS)

## add all subdirectories
add_subdirectory(lib)
add_subdirectory(glauber)
add_subdirectory(utils)
add_subdirectory(analysis)
add_subdirectory(centrality)
add_subdirectory(systematics)

## handle any files that have to be configured
configure_file(${PROJECT_SOURCE_DIR}/sct/lib/macros.h.in
               ${PROJECT_BINARY_DIR}/sct/lib/macros.h)
install(FILES ${PROJECT_BINARY_DIR}/sct/lib/macros.h
        DESTINATION include/sct/lib)

## compile libraries
add_library(sct ${SCT_SRCS})
target_link_libraries(sct PRIVATE ${SCT_DEPENDENCY_LIBS})
target_include_directories(sct INTERFACE $<INSTALL_INTERFACE:include>)
target_compile_options(sct INTERFACE "--std=c++11")
install(TARGETS sct EXPORT SCTTargets DESTINATION lib)
link_directories(${CMAKE_INSTALL_PREFIX}/lib)

set(SCT_LIBS sct)

## build sct binaries
if (BUILD_BINARIES)
  foreach(bin_src ${SCT_BINARY_SRCS})
    get_filename_component(bin_name ${bin_src} NAME_WE)
    add_executable(${bin_name} "${bin_src}")
    add_dependencies(${bin_name} ${SCT_LIBS})
    target_link_libraries(${bin_name} ${SCT_LIBS}
                        ${SCT_DEPENDENCY_LIBS})
    if (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.0)
      target_compile_features(${bin_name} PRIVATE cxx_range_for)
    endif()
    add_test(NAME ${bin_name} COMMAND $<TARGET_FILE:${bin_name}>)
    install(TARGETS ${bin_name} DESTINATION bin)
  endforeach()
endif()

## test binaries
if (BUILD_TEST)
  foreach(test_src ${SCT_TEST_SRCS})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} "${test_src}")
    add_dependencies(${test_name} ${SCT_LIBS})
    target_link_libraries(${test_name} ${SCT_LIBS}
                    ${SCT_DEPENDENCY_LIBS} gtest_main)
    if (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.0)
      target_compile_features(${test_name} PRIVATE cxx_range_for)
    endif()
    add_test(NAME ${test_name} COMMAND $<TARGET_FILE:${test_name}>)
    install(TARGETS ${test_name} DESTINATION test)
  endforeach()
endif()

## benchmark binaries
if (BUILD_TEST)
  foreach(bm_src ${SCT_BENCH_SRCS})
    get_filename_component(bm_name ${bm_src} NAME_WE)
    add_executable(${bm_name} "${bm_src}")
    add_dependencies(${bm_name} ${SCT_LIBS})
    target_link_libraries(${bm_name} ${SCT_LIBS}
                  ${SCT_DEPENDENCY_LIBS})
    if (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.0)
      target_compile_features(${bm_name} PRIVATE cxx_range_for)
    endif()
    add_test(NAME ${bm_name} COMMAND $<TARGET_FILE:${bm_name}>)
    install(TARGETS ${bm_name} DESTINATION test)
  endforeach()
endif()
