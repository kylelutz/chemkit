find_package(Chemkit COMPONENTS io OPTIONAL)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  sybylatomtyper.cpp
  sybylplugin.cpp
)

if(${CHEMKIT_WITH_IO})
  set(SOURCES ${SOURCES} mol2fileformat.cpp)
endif()

add_chemkit_plugin(sybyl ${SOURCES})
target_link_libraries(sybyl ${CHEMKIT_LIBRARIES})
