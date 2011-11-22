find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  countdescriptors.cpp
  countdescriptorsplugin.cpp
)

add_chemkit_plugin(countdescriptors ${SOURCES})
target_link_libraries(countdescriptors ${CHEMKIT_LIBRARIES})
