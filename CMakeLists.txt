cmake_minimum_required(VERSION 2.8)

project(chemkit)

set(chemkit_VERSION_MAJOR 0)
set(chemkit_VERSION_MINOR 0)

set(CMAKE_BUILD_TYPE Release)

add_subdirectory(src)
add_subdirectory(examples)

enable_testing()
add_subdirectory(tests)
