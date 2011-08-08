find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES
  randicindexdescriptor.cpp
  randicindexplugin.cpp
)

add_chemkit_plugin(randicindex ${SOURCES})
target_link_libraries(randicindex chemkit)
