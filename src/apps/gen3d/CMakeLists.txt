if(NOT ${CHEMKIT_WITH_IO} OR NOT ${CHEMKIT_WITH_MD})
  return()
endif()

find_package(Chemkit COMPONENTS io md REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Boost COMPONENTS program_options REQUIRED)

add_chemkit_executable(gen3d gen3d.cpp)
target_link_libraries(gen3d ${CHEMKIT_LIBRARIES} ${Boost_LIBRARIES})
