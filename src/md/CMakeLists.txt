if(NOT ${CHEMKIT_WITH_MD})
  return()
endif()

find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(HEADERS
  forcefieldcalculation.h
  forcefieldenergydescriptor.h
  forcefieldenergydescriptor-inline.h
  forcefield.h
  integrator.h
  md.h
  moleculegeometryoptimizer.h
  potential.h
  topology.h
  topologybuilder.h
  trajectory.h
  trajectoryframe.h
)

set(SOURCES
  forcefieldcalculation.cpp
  forcefield.cpp
  integrator.cpp
  md.cpp
  moleculegeometryoptimizer.cpp
  potential.cpp
  topology.cpp
  topologybuilder.cpp
  trajectory.cpp
  trajectoryframe.cpp
)

add_definitions(
  -DCHEMKIT_MD_LIBRARY
)

add_chemkit_library(chemkit-md ${SOURCES})
target_link_libraries(chemkit-md ${CHEMKIT_LIBRARIES})

# install header files
install(FILES ${HEADERS} DESTINATION include/chemkit/)
