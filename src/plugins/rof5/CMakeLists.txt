find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  rof5plugin.cpp
  ruleoffivedescriptor.cpp
  ruleoffiveviolationsdescriptor.cpp
)

add_chemkit_plugin(rof5 ${SOURCES})
target_link_libraries(rof5 ${CHEMKIT_LIBRARIES})
