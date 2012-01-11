if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  pdbmlfileformat.cpp
  pdbmlplugin.cpp
)

add_chemkit_plugin(pdbml ${SOURCES})
target_link_libraries(pdbml ${CHEMKIT_LIBRARIES})
