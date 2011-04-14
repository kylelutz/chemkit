
# set compiler warning flags
if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()

add_subdirectory(chemkit)
add_subdirectory(graphics)
add_subdirectory(web)
add_subdirectory(widgets)
add_subdirectory(plugins)
add_subdirectory(qt-designer-plugins)
add_subdirectory(apps)
