if(NOT ${CHEMKIT_BUILD_BINDINGS_PYTHON})
  return()
endif()

macro(add_chemkit_python_test test_name test_script)
  add_test(${test_name} ${PYTHON_EXECUTABLE} ${test_script})
  set_tests_properties(${test_name} PROPERTIES
                       ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/bindings/python/;CHEMKIT_PLUGIN_PATH=${CMAKE_BINARY_DIR}/lib/chemkit/plugins")
endmacro(add_chemkit_python_test)

add_chemkit_python_test(python.Atom ${CMAKE_CURRENT_SOURCE_DIR}/atomtest.py)
add_chemkit_python_test(python.Bond ${CMAKE_CURRENT_SOURCE_DIR}/bondtest.py)
add_chemkit_python_test(python.Element ${CMAKE_CURRENT_SOURCE_DIR}/elementtest.py)
add_chemkit_python_test(python.Fingerprint ${CMAKE_CURRENT_SOURCE_DIR}/fingerprinttest.py)
add_chemkit_python_test(python.ForceField ${CMAKE_CURRENT_SOURCE_DIR}/forcefieldtest.py)
add_chemkit_python_test(python.Fragment ${CMAKE_CURRENT_SOURCE_DIR}/fragmenttest.py)
add_chemkit_python_test(python.Molecule ${CMAKE_CURRENT_SOURCE_DIR}/moleculetest.py)
add_chemkit_python_test(python.MoleculeFile ${CMAKE_CURRENT_SOURCE_DIR}/moleculefiletest.py)
add_chemkit_python_test(python.Ring ${CMAKE_CURRENT_SOURCE_DIR}/ringtest.py)
