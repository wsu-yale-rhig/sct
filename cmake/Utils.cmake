# anaconda has quite a few packages that can be in conflict with the OS - this
# script attempts to avoid the problem, we always put anaconda include paths
# after system include paths

function(sct_include_directories)
  foreach(dir IN LISTS ARGN)
    if(${dir} MATCHES "/anaconda")
      include_directories(AFTER SYSTEM ${dir})
    else()
      include_directories(BEFORE SYSTEM ${dir})
    endif()
  endforeach()
endfunction()
