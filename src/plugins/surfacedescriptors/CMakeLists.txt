find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  surfacedescriptors.cpp
  surfacedescriptorsplugin.cpp
)

add_chemkit_plugin(surfacedescriptors ${SOURCES})
target_link_libraries(surfacedescriptors ${CHEMKIT_LIBRARIES})
