if(NOT ${CHEMKIT_BUILD_DEMOS})
  return()
endif()

add_subdirectory(cube-viewer)
add_subdirectory(pdb-viewer)
add_subdirectory(surface-viewer)
add_subdirectory(trajectory-viewer)
