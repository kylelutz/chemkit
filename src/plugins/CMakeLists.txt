macro(add_chemkit_plugin plugin_name)
  # add library
  add_library(${plugin_name} SHARED ${ARGN})

  # add install target
  install(TARGETS ${plugin_name} DESTINATION lib/chemkit/plugins/)

  # remove 'lib' prefix
  set_target_properties(${plugin_name} PROPERTIES PREFIX "")

  # copy plugin into build directory
  get_target_property(plugin_location ${plugin_name} LOCATION)
  get_filename_component(plugin_filename ${plugin_location} NAME)
  add_custom_command(TARGET ${plugin_name} POST_BUILD COMMAND ${CMAKE_COMMAND} ARGS -E copy ${plugin_location} ${CMAKE_BINARY_DIR}/lib/chemkit/plugins/${plugin_filename})

  # plugin data
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/data)
    # copy plugin data into build directory
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data/ DESTINATION ${CMAKE_BINARY_DIR}/lib/chemkit/plugins/data/${plugin_name}/)
    
    # install plugin data directory
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data/ DESTINATION lib/chemkit/plugins/data/${plugin_name})
  endif()
endmacro(add_chemkit_plugin)

add_subdirectory(amber)
add_subdirectory(cas)
add_subdirectory(cml)
add_subdirectory(countdescriptors)
add_subdirectory(cube)
add_subdirectory(fhz)
add_subdirectory(formula)
add_subdirectory(fp2)
add_subdirectory(fps)
add_subdirectory(gasteiger)
add_subdirectory(graphdescriptors)
add_subdirectory(gromacs)
add_subdirectory(hydrogenbonds)
add_subdirectory(inchi)
add_subdirectory(massdescriptors)
add_subdirectory(mcdl)
add_subdirectory(mdl)
add_subdirectory(mlogp)
add_subdirectory(mmff)
add_subdirectory(mopac)
add_subdirectory(opls)
add_subdirectory(pdb)
add_subdirectory(pqr)
add_subdirectory(pubchem)
add_subdirectory(randicindex)
add_subdirectory(rof5)
add_subdirectory(rotatablebonds)
add_subdirectory(smiles)
add_subdirectory(surfacedescriptors)
add_subdirectory(sybyl)
add_subdirectory(tpsa)
add_subdirectory(txyz)
add_subdirectory(uff)
add_subdirectory(vabc)
add_subdirectory(wienerindex)

# only compile xtc plugin on *nix where we have XDR
if(${CHEMKIT_OS_UNIX} AND NOT ${CHEMKIT_OS_MAC})
  add_subdirectory(xtc)
endif()

add_subdirectory(xyz)
