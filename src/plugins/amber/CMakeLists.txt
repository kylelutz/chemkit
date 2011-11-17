find_package(Chemkit COMPONENTS md REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  ambercalculation.cpp
  amberforcefield.cpp
  amberparameters.cpp
  amberplugin.cpp
)

add_chemkit_plugin(amber ${SOURCES})
target_link_libraries(amber ${CHEMKIT_LIBRARIES})
