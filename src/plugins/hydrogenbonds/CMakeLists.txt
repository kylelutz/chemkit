find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  hydrogenbondacceptorsdescriptor.cpp
  hydrogenbonddonorsdescriptor.cpp
  hydrogenbondsplugin.cpp
)

add_chemkit_plugin(hydrogenbonds ${SOURCES})
target_link_libraries(hydrogenbonds ${CHEMKIT_LIBRARIES})
