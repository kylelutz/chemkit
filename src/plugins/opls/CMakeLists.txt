find_package(Chemkit COMPONENTS md REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  oplsatomtyper.cpp
  oplscalculation.cpp
  oplsforcefield.cpp
  oplsparameters.cpp
  oplsplugin.cpp
)

add_chemkit_plugin(opls ${SOURCES})
target_link_libraries(opls ${CHEMKIT_LIBRARIES})
