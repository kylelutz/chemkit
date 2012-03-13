find_package(Chemkit)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  tpsadescriptor.cpp
  tpsaplugin.cpp
)

add_chemkit_plugin(tpsa ${SOURCES})
target_link_libraries(tpsa ${CHEMKIT_LIBRARIES})
