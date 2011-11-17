find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  randicindexdescriptor.cpp
  randicindexplugin.cpp
)

add_chemkit_plugin(randicindex ${SOURCES})
target_link_libraries(randicindex ${CHEMKIT_LIBRARIES})
