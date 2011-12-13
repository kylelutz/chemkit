if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  mdlfileformat.cpp
  mdlplugin.cpp
)

add_chemkit_plugin(mdl ${SOURCES})
target_link_libraries(mdl ${CHEMKIT_LIBRARIES})
