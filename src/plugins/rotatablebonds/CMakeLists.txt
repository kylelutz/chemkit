find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES
  rotatablebondsdescriptor.cpp
  rotatablebondsplugin.cpp
)

add_chemkit_plugin(rotatablebonds ${SOURCES})
target_link_libraries(rotatablebonds chemkit)
