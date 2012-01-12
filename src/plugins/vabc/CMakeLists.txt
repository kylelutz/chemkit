find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  vabcdescriptor.cpp
  vabcplugin.cpp
)

add_chemkit_plugin(vabc ${SOURCES})
target_link_libraries(vabc ${CHEMKIT_LIBRARIES})
