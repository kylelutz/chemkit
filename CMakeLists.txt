cmake_minimum_required(VERSION 2.8)

project(Chemkit)

set(CHEMKIT_VERSION_MAJOR 0)
set(CHEMKIT_VERSION_MINOR 0)

set(CMAKE_MODULE_PATH ${Chemkit_SOURCE_DIR}/cmake)

# find dependencies
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)

# add option for CMAKE_BUILD_TYPE
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# uninstall target
configure_file(
	"${CMAKE_MODULE_PATH}/cmake_uninstall.cmake.in"
	"${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
	IMMEDIATE @ONLY)

add_custom_target(uninstall 
	COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

# configure and install ChemkitConfig.cmake file
configure_file(
    "${CMAKE_MODULE_PATH}/ChemkitConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/ChemkitConfig.cmake"
    IMMEDIATE @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/ChemkitConfig.cmake DESTINATION lib/chemkit)

# Create a ChemkitBuildTreeSettings.cmake file for the use from the build tree
configure_file(
    "${CMAKE_MODULE_PATH}/ChemkitBuildTreeSettings.cmake.in"
    "${CMAKE_BINARY_DIR}/ChemkitBuildTreeSettings.cmake"
    IMMEDIATE @ONLY)

# export the Chemkit package
export(PACKAGE Chemkit)

# build options
option(CHEMKIT_WITH_GRAPHICS "Build the chemkit-graphics library." ON)
option(CHEMKIT_WITH_IO "Build the chemkit-io library." ON)
option(CHEMKIT_WITH_MD "Build the chemkit-md library." ON)
option(CHEMKIT_WITH_WEB "Build the chemkit-web library." ON)
option(CHEMKIT_WITH_WIDGETS "Build the chemkit-widgets library." ON)
option(CHEMKIT_BUILD_APPS "Build the chemkit applications" ON)
option(CHEMKIT_BUILD_QT_DESIGNER_PLUGINS "Build the chemkit QtDesigner plugins." OFF)
option(CHEMKIT_BUILD_EXAMPLES "Build the chemkit examples." OFF)
option(CHEMKIT_BUILD_TESTS "Build the chemkit tests." OFF)

add_subdirectory(src)
add_subdirectory(examples)

enable_testing()
add_subdirectory(tests)
