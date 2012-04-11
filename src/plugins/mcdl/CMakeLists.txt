find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  mcdllineformat.cpp
  mcdlplugin.cpp
  mcdlreader.cpp
)

add_chemkit_plugin(mcdl ${SOURCES})
target_link_libraries(mcdl ${CHEMKIT_LIBRARIES})
