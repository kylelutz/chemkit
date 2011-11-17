find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

set(HEADERS
  genericfile.h
  genericfile-inline.h
  io.h
  moleculefile.h
  moleculefileformat.h
  moleculefileformatadaptor.h
  moleculefileformatadaptor-inline.h
  polymerfile.h
  polymerfileformat.h
)

set(SOURCES
  io.cpp
  moleculefile.cpp
  moleculefileformat.cpp
  polymerfile.cpp
  polymerfileformat.cpp
)

add_definitions(
  -DCHEMKIT_IO_LIBRARY
)

add_chemkit_library(chemkit-io ${SOURCES})
target_link_libraries(chemkit-io ${CHEMKIT_LIBRARIES})

# install header files
install(FILES ${HEADERS} DESTINATION include/chemkit/)
