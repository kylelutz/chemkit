find_package(Chemkit COMPONENTS io OPTIONAL)
include_directories(${CHEMKIT_INCLUDE_DIRS})

include_directories(../../3rdparty/lemon)

set(SOURCES
  daylightaromaticitymodel.cpp
  kekulizer.cpp
  smilesgraph.cpp
  smileslineformat.cpp
  smilesplugin.cpp
)

add_chemkit_plugin(smiles ${SOURCES})
target_link_libraries(smiles ${CHEMKIT_LIBRARIES})
