if(NOT ${CHEMKIT_WITH_MD_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io md md-io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  grofileformat.cpp
  gromacsplugin.cpp
  topfileformat.cpp
)

add_chemkit_plugin(gromacs ${SOURCES})
target_link_libraries(gromacs ${CHEMKIT_LIBRARIES})
