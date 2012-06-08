cmake_minimum_required(VERSION 2.8)

project(Chemkit)

set(CHEMKIT_VERSION_MAJOR 0)
set(CHEMKIT_VERSION_MINOR 1)

set(CMAKE_MODULE_PATH ${Chemkit_SOURCE_DIR}/cmake)

# find dependencies
find_package(Eigen3 REQUIRED)

if(WIN32)
  set(Boost_USE_STATIC_LIBS ON)
endif()

find_package(Boost REQUIRED)
link_directories(${Boost_LIBRARY_DIRS})

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

set(CHEMKIT_IN_BUILD_TREE TRUE)

# so "find_package(Chemkit)" calls will find this build tree's version
set(Chemkit_DIR "${CMAKE_BINARY_DIR}")

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
option(CHEMKIT_WITH_MD_IO "Build the chemkit-md-io library." ON)
option(CHEMKIT_WITH_WEB "Build the chemkit-web library." ON)
option(CHEMKIT_WITH_GUI "Build the chemkit-gui library." ON)
option(CHEMKIT_BUILD_APPS "Build the chemkit applications" ON)
option(CHEMKIT_BUILD_QT_DESIGNER_PLUGINS "Build the chemkit QtDesigner plugins." OFF)
option(CHEMKIT_BUILD_DEMOS "Build the chemkit demos." OFF)
option(CHEMKIT_BUILD_EXAMPLES "Build the chemkit examples." OFF)
option(CHEMKIT_BUILD_TESTS "Build the chemkit tests." OFF)

# compiler options
if(MSVC)
  # disable dll-linkage warnings from MSVC
  add_definitions("/wd4251")

  # disable boost auto-linking
  add_definitions("-DBOOST_ALL_NO_LIB")

  # disable unsafe string function warnings
  add_definitions("-D_CRT_SECURE_NO_WARNINGS")
endif()

# set a variable for the operating system
set(CHEMKIT_OS_UNIX FALSE)
set(CHEMKIT_OS_MAC FALSE)
set(CHEMKIT_OS_WIN32 FALSE)

if(UNIX)
  set(CHEMKIT_OS_UNIX TRUE)
  if(APPLE)
    set(CHEMKIT_OS_MAC TRUE)
  endif()
elseif(WIN32)
  set(CHEMKIT_OS_WIN32 TRUE)
endif()

add_subdirectory(src)
add_subdirectory(bindings)
add_subdirectory(demos)
add_subdirectory(examples)

enable_testing()
add_subdirectory(tests)
