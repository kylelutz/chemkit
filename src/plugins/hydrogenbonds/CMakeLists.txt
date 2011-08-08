find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES
  hydrogenbondacceptorsdescriptor.cpp
  hydrogenbonddonorsdescriptor.cpp
  hydrogenbondsplugin.cpp
)

add_chemkit_plugin(hydrogenbonds ${SOURCES})
target_link_libraries(hydrogenbonds chemkit)
