if(NOT ${CHEMKIT_WITH_MD_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io md REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(HEADERS
  topologyfile.h
  topologyfileformat.h
  trajectoryfile.h
  trajectoryfileformat.h
)

set(SOURCES
  topologyfile.cpp
  topologyfileformat.cpp
  trajectoryfile.cpp
  trajectoryfileformat.cpp
)

add_definitions(
  -DCHEMKIT_MD_IO_LIBRARY
)

add_chemkit_library(chemkit-md-io ${SOURCES})
target_link_libraries(chemkit-md-io ${CHEMKIT_LIBRARIES})

# install header files
install(FILES ${HEADERS} DESTINATION include/chemkit/)
