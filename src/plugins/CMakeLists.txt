include_directories(../../include)

macro(add_chemkit_plugin plugin_name)
  # add library
  add_library(${plugin_name} SHARED ${ARGN})

  # add install target
  install(TARGETS ${plugin_name} DESTINATION share/chemkit/plugins/)

  # copy plugin into build directory
  get_target_property(plugin_location ${plugin_name} LOCATION)
  get_filename_component(plugin_filename ${plugin_location} NAME)
  add_custom_command(TARGET ${plugin_name} POST_BUILD COMMAND ${CMAKE_COMMAND} ARGS -E copy ${plugin_location} ${CMAKE_BINARY_DIR}/share/chemkit/plugins/${plugin_filename})

  # plugin data
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/data)
    # copy plugin data into build directory
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data/ DESTINATION ${CMAKE_BINARY_DIR}/share/chemkit/plugins/data/${plugin_name}/)
    
    # install plugin data directory
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data/ DESTINATION share/chemkit/plugins/data/${plugin_name})
  endif()
endmacro(add_chemkit_plugin)

add_subdirectory(amber)
add_subdirectory(cml)
add_subdirectory(cube)
add_subdirectory(fhz)
add_subdirectory(formula)
add_subdirectory(gasteiger)
add_subdirectory(hydrogenbonds)
add_subdirectory(inchi)
add_subdirectory(mcdl)
add_subdirectory(mdl)
add_subdirectory(mmff)
add_subdirectory(mopac)
add_subdirectory(opls)
add_subdirectory(pdb)
add_subdirectory(pdbml)
add_subdirectory(randicindex)
add_subdirectory(rotatablebonds)
add_subdirectory(smiles)
add_subdirectory(sybyl)
add_subdirectory(txyz)
add_subdirectory(uff)
add_subdirectory(wienerindex)
add_subdirectory(xtc)
add_subdirectory(xyz)
