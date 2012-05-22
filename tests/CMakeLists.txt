if(NOT ${CHEMKIT_BUILD_TESTS})
  return()
endif()

macro(add_chemkit_test test_name test_executable)
  add_test(NAME ${test_name} COMMAND ${test_executable})
  set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "CHEMKIT_PLUGIN_PATH=${CMAKE_BINARY_DIR}/lib/chemkit/plugins/")
endmacro(add_chemkit_test)

add_subdirectory(auto)
add_subdirectory(benchmarks)
add_subdirectory(manual)

# copy test data directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
