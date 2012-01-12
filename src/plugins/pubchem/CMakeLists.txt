find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  pubchemfingerprint.cpp
  pubchemplugin.cpp
)

add_chemkit_plugin(pubchem ${SOURCES})
target_link_libraries(pubchem ${CHEMKIT_LIBRARIES})
