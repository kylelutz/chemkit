if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  txyzfileformat.cpp
  txyzplugin.cpp
)

add_chemkit_plugin(txyz ${SOURCES})
target_link_libraries(txyz ${CHEMKIT_LIBRARIES})
