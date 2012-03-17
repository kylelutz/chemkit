find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  mlogpdescriptor.cpp
  mlogpplugin.cpp
)

add_chemkit_plugin(mlogp ${SOURCES})
target_link_libraries(mlogp ${CHEMKIT_LIBRARIES})
