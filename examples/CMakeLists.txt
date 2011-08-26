if(NOT ${CHEMKIT_BUILD_EXAMPLES})
  return()
endif()

add_subdirectory(cube-viewer)
add_subdirectory(pdb-viewer)
add_subdirectory(uff-energy)
