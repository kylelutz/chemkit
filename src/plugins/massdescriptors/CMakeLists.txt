find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  massdescriptors.cpp
  massdescriptorsplugin.cpp
)

add_chemkit_plugin(massdescriptors ${SOURCES})
target_link_libraries(massdescriptors ${CHEMKIT_LIBRARIES})
