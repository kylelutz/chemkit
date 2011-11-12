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
)

set(SOURCES
  forcefieldatom.cpp
  forcefieldcalculation.cpp
  forcefield.cpp
  forcefieldinteractions.cpp
  md.cpp
)

add_definitions(
  -DCHEMKIT_MD_LIBRARY
)

include_directories(../../include)

add_chemkit_library(chemkit-md ${SOURCES})
target_link_libraries(chemkit-md chemkit)

# install header files
install(FILES ${HEADERS} DESTINATION include/chemkit/)
