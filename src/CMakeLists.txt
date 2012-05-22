# set compiler options
if(CMAKE_COMPILER_IS_GNUCXX AND NOT WIN32)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden -fvisibility-inlines-hidden -pedantic -Wall -Wextra -Wno-long-long")
endif()

macro(add_chemkit_library library_name)
  # add library
  add_library(${library_name} SHARED ${ARGN})

  # add install target
  install(TARGETS ${library_name} DESTINATION lib)

  # copy library into build directory
  get_target_property(library_location ${library_name} LOCATION)
  get_filename_component(library_filename ${library_location} NAME)
  add_custom_command(TARGET ${library_name} POST_BUILD COMMAND ${CMAKE_COMMAND} ARGS -E copy ${library_location} ${CMAKE_BINARY_DIR}/lib/${library_filename})

endmacro(add_chemkit_library)

add_subdirectory(chemkit)
add_subdirectory(graphics)
add_subdirectory(io)
add_subdirectory(md)
add_subdirectory(md-io)
add_subdirectory(web)
add_subdirectory(gui)
add_subdirectory(plugins)
add_subdirectory(qt-designer-plugins)
add_subdirectory(apps)
