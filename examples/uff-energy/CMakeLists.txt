if(NOT ${CHEMKIT_WITH_IO} OR NOT ${CHEMKIT_WITH_MD})
  return()
endif()

find_package(Chemkit COMPONENTS io md)
include_directories(${CHEMKIT_INCLUDE_DIRS})

add_executable(uff-energy uffenergyexample.cpp)
target_link_libraries(uff-energy ${CHEMKIT_LIBRARIES})
