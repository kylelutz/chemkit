if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

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

include_directories(../../include)

add_chemkit_library(chemkit-io ${SOURCES})
target_link_libraries(chemkit-io chemkit)

# install header files
install(FILES ${HEADERS} DESTINATION include/chemkit/)
