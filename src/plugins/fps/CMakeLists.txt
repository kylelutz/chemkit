if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  fpsfileformat.cpp
  fpsplugin.cpp
)

add_chemkit_plugin(fps ${SOURCES})
target_link_libraries(fps ${CHEMKIT_LIBRARIES})
