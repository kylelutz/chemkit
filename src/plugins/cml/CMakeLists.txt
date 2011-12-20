if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  cmlfileformat.cpp
  cmlplugin.cpp
)

add_chemkit_plugin(cml ${SOURCES})
target_link_libraries(cml ${CHEMKIT_LIBRARIES})
