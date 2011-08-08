find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES
  cubefileformat.cpp
  cubeplugin.cpp
)

add_chemkit_plugin(cube ${SOURCES})
target_link_libraries(cube chemkit)
