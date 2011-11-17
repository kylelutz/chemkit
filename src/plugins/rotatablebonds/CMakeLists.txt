find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  rotatablebondsdescriptor.cpp
  rotatablebondsplugin.cpp
)

add_chemkit_plugin(rotatablebonds ${SOURCES})
target_link_libraries(rotatablebonds ${CHEMKIT_LIBRARIES})
