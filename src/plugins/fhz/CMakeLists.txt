if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  fhzfileformat.cpp
  fhzplugin.cpp
)

add_chemkit_plugin(fhz ${SOURCES})
target_link_libraries(fhz ${CHEMKIT_LIBRARIES})
