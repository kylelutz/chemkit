find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  atomicnumberatomtyper.cpp
  elementnameatomtyper.cpp
  elementsymbolatomtyper.cpp
  elementtypersplugin.cpp
)

add_chemkit_plugin(elementtypers ${SOURCES})
target_link_libraries(elementtypers ${CHEMKIT_LIBRARIES})
