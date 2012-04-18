if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  cubefileformat.cpp
  gaussianplugin.cpp
)

add_chemkit_plugin(gaussian ${SOURCES})
target_link_libraries(gaussian ${CHEMKIT_LIBRARIES})
