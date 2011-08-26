if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES
  xyzfileformat.cpp
  xyzplugin.cpp
)

add_chemkit_plugin(xyz ${SOURCES})
target_link_libraries(xyz chemkit chemkit-io)
