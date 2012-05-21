find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  apoldescriptor.cpp
  apolplugin.cpp
)

add_chemkit_plugin(apol ${SOURCES})
target_link_libraries(apol ${CHEMKIT_LIBRARIES})
