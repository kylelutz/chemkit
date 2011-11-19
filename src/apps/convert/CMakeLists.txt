if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Boost COMPONENTS program_options REQUIRED)

add_chemkit_executable(convert convert.cpp)
target_link_libraries(convert ${CHEMKIT_LIBRARIES} ${Boost_LIBRARIES})
