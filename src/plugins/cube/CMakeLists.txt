if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  cubefileformat.cpp
  cubeplugin.cpp
)

add_chemkit_plugin(cube ${SOURCES})
target_link_libraries(cube ${CHEMKIT_LIBRARIES})
