find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

include_directories(../../3rdparty/lemon)

set(SOURCES
  kekulizer.cpp
  smilesgraph.cpp
  smileslineformat.cpp
  smilesplugin.cpp
)

add_chemkit_plugin(smiles ${SOURCES})
target_link_libraries(smiles chemkit)

if(${CHEMKIT_WITH_IO})
  target_link_libraries(smiles chemkit-io)
endif()
