if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  pqrfileformat.cpp
  pqrplugin.cpp
)

add_chemkit_plugin(pqr ${SOURCES})
target_link_libraries(pqr ${CHEMKIT_LIBRARIES})
