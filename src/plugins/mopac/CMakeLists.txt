if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  mopacplugin.cpp
  mopcrtfileformat.cpp
  mopinfileformat.cpp
)

add_chemkit_plugin(mopac ${SOURCES})
target_link_libraries(mopac ${CHEMKIT_LIBRARIES})
