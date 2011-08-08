find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES
  ambercalculation.cpp
  amberforcefield.cpp
  amberparameters.cpp
  amberplugin.cpp
)

add_chemkit_plugin(amber ${SOURCES})
target_link_libraries(amber chemkit)
