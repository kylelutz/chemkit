find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES
  uffatomtyper.cpp
  uffcalculation.cpp
  uffforcefield.cpp
  uffparameters.cpp
  uffplugin.cpp
)

add_chemkit_plugin(uff ${SOURCES})
target_link_libraries(uff chemkit)
