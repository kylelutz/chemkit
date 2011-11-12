if(NOT ${CHEMKIT_WITH_MD})
  return()
endif()

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(HEADERS
  forcefieldatom.h
  forcefieldcalculation.h
  forcefieldcalculation-inline.h
  forcefield.h
  forcefieldinteractions.h
  md.h
  trajectory.h
  trajectoryframe.h
)

set(SOURCES
  forcefieldatom.cpp
  forcefieldcalculation.cpp
  forcefield.cpp
  forcefieldinteractions.cpp
  md.cpp
  trajectory.cpp
  trajectoryframe.cpp
)

# add md io classes if building with io
if(${CHEMKIT_WITH_IO})
  set(HEADERS ${HEADERS}
    trajectoryfile.h
    trajectoryfileformat.h
  )

  set(SOURCES ${SOURCES}
    trajectoryfile.cpp
    trajectoryfileformat.cpp
   )
endif()

add_definitions(
  -DCHEMKIT_MD_LIBRARY
)

include_directories(../../include)

add_chemkit_library(chemkit-md ${SOURCES})
target_link_libraries(chemkit-md chemkit)

if(${CHEMKIT_WITH_IO})
  target_link_libraries(chemkit-md chemkit-io)
endif()

# install header files
install(FILES ${HEADERS} DESTINATION include/chemkit/)
