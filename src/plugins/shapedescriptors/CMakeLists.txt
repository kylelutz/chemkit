find_package(Chemkit)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  gravitationalindexdescriptor.cpp
  radiusofgyrationdescriptor.cpp
  shapedescriptorsplugin.cpp
)

add_chemkit_plugin(shapedescriptors ${SOURCES})
target_link_libraries(shapedescriptors ${CHEMKIT_LIBRARIES})
