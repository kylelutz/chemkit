if(NOT ${CHEMKIT_BUILD_APPS})
  return()
endif()

macro(add_chemkit_executable name)
  add_executable(${name} ${ARGN})
  set_target_properties(${name} PROPERTIES OUTPUT_NAME chemkit-${name})
  install(TARGETS ${name} DESTINATION bin)

  # copy binary
  get_target_property(location ${name} LOCATION)
  get_filename_component(filename ${location} NAME)
  add_custom_command(TARGET ${name} POST_BUILD COMMAND ${CMAKE_COMMAND} ARGS -E copy ${location} ${CMAKE_BINARY_DIR}/bin/chemkit-${name})
endmacro(add_chemkit_executable)

add_subdirectory(convert)
add_subdirectory(gen3d)
add_subdirectory(grep)
add_subdirectory(translate)
