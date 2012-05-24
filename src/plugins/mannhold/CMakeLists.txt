find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  mannholdlogpdescriptor.cpp
  mannholdplugin.cpp
)

add_chemkit_plugin(mannhold ${SOURCES})
target_link_libraries(mannhold ${CHEMKIT_LIBRARIES})
