find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES
  wienerindexdescriptor.cpp
  wienerindexplugin.cpp
)

add_chemkit_plugin(wienerindex ${SOURCES})
target_link_libraries(wienerindex chemkit)
