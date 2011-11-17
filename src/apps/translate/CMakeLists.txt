find_package(Chemkit)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Boost COMPONENTS program_options REQUIRED)

add_chemkit_executable(translate translate.cpp)
target_link_libraries(translate ${CHEMKIT_LIBRARIES} ${Boost_LIBRARIES})
