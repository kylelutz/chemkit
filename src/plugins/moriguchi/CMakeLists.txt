find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  moriguchilogpdescriptor.cpp
  moriguchiplugin.cpp
)

add_chemkit_plugin(moriguchi ${SOURCES})
target_link_libraries(moriguchi ${CHEMKIT_LIBRARIES})
