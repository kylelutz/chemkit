if(NOT ${CHEMKIT_WITH_MD})
  return()
endif()

find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(HEADERS
  forcefieldatom.h
  forcefieldcalculation.h
  forcefieldcalculation-inline.h
  forcefieldenergydescriptor.h
  forcefieldenergydescriptor-inline.h
  forcefield.h
  forcefieldinteractions.h
  md.h
  moleculegeometryoptimizer.h
  topology.h
  trajectory.h
  trajectoryframe.h
)

set(SOURCES
  forcefieldatom.cpp
  forcefieldcalculation.cpp
  forcefield.cpp
  forcefieldinteractions.cpp
  md.cpp
  moleculegeometryoptimizer.cpp
  topology.cpp
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
