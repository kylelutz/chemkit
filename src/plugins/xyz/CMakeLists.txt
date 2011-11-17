if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  xyzfileformat.cpp
  xyzplugin.cpp
)

add_chemkit_plugin(xyz ${SOURCES})
target_link_libraries(xyz ${CHEMKIT_LIBRARIES})
