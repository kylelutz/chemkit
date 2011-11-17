find_package(Chemkit)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  wienerindexdescriptor.cpp
  wienerindexplugin.cpp
)

add_chemkit_plugin(wienerindex ${SOURCES})
target_link_libraries(wienerindex ${CHEMKIT_LIBRARIES})
