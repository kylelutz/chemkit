find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

include_directories(../../../include)

add_chemkit_executable(translate translate.cpp)
target_link_libraries(translate chemkit ${Boost_LIBRARIES})
