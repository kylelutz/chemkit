if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  pdbfileformat.cpp
  pdbmlfileformat.cpp
  pdbplugin.cpp
)

add_chemkit_plugin(pdb ${SOURCES})
target_link_libraries(pdb ${CHEMKIT_LIBRARIES})
