# get source files for libstglauber
file(GLOB tmp *.cxx)
set(STGLAUBER_SRCS ${STGLAUBER_SRCS} ${tmp})

# get header files for libstglauber
file(GLOB tmp *.h)
set(STGLAUBER_HEADERS ${tmp})

# install headers
install(FILES ${STGLAUBER_HEADERS} DESTINATION include/stglauber/StGlauberUtilities)

#export include path
set(STGLAUBER_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR} ${STGLAUBER_INCLUDE})

# export file lists to parent scope
set(STGLAUBER_SRCS ${STGLAUBER_SRCS} PARENT_SCOPE)
set(STGLAUBER_INCLUDE ${STGLAUBER_INCLUDE} PARENT_SCOPE)