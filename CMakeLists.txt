cmake_minimum_required(VERSION 2.8)

project(chemkit)

set(chemkit_VERSION_MAJOR 0)
set(chemkit_VERSION_MINOR 0)

set(CMAKE_MODULE_PATH ${chemkit_SOURCE_DIR}/cmake)

set(CMAKE_BUILD_TYPE Release)

# find eigen3
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_subdirectory(src)
add_subdirectory(examples)

enable_testing()
add_subdirectory(tests)
