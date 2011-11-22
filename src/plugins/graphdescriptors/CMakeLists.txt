find_package(Chemkit)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  graphdescriptors.cpp
  graphdescriptorsplugin.cpp
)

add_chemkit_plugin(graphdescriptors ${SOURCES})
target_link_libraries(graphdescriptors ${CHEMKIT_LIBRARIES})
