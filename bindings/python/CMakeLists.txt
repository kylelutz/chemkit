option(CHEMKIT_BUILD_BINDINGS_PYTHON "Built the Python bindings" OFF)
if(NOT ${CHEMKIT_BUILD_BINDINGS_PYTHON})
  return()
endif()

# find chemkit
find_package(Chemkit COMPONENTS io md REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

# find python
find_package(PythonLibs REQUIRED)
find_package(PythonInterp REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

# find cython
find_program(CYTHON cython REQUIRED)

set(SOURCES
  atom.pxd
  atom.pxi
  atomtyper.pxd
  atomtyper.pxi
  bond.pxd
  bond.pxi
  bondpredictor.pxd
  bondpredictor.pxi
  chemkit.pyx
  coordinatepredictor.pxd
  coordinatepredictor.pxi
  element.pxd
  element.pxi
  fingerprint.pxd
  fingerprint.pxi
  forcefield.pxd
  forcefield.pxi
  fragment.pxd
  fragment.pxi
  lineformat.pxd
  lineformat.pxi
  moleculardescriptor.pxd
  moleculardescriptor.pxi
  moleculefile.pxd
  moleculefile.pxi
  molecule.pxd
  molecule.pxi
  moleculegeometryoptimizer.pxd
  moleculegeometryoptimizer.pxi
  partialchargemodel.pxd
  partialchargemodel.pxi
  point3.pxd
  point3.pxi
  ring.pxd
  ring.pxi
  shared_ptr.pxd
  string.pxd
  variant.pxd
  vector3.pxd
  vector3.pxi
)

# run cython
add_custom_command(OUTPUT chemkit.cpp
                   COMMAND ${CYTHON} --cplus --fast-fail
                     -o ${CMAKE_CURRENT_BINARY_DIR}/chemkit.cpp
                     ${CMAKE_CURRENT_SOURCE_DIR}/chemkit.pyx
                   DEPENDS ${SOURCES})

# build library
add_library(chemkit-python SHARED chemkit.cpp)
set_target_properties(chemkit-python PROPERTIES OUTPUT_NAME "chemkit" PREFIX "")
target_link_libraries(chemkit-python ${CHEMKIT_LIBRARIES} ${PYTHON_LIBRARIES})

# find python site-packages directory for installation
execute_process(COMMAND ${PYTHON_EXECUTABLE}
                -c "from distutils.sysconfig import get_python_lib; print get_python_lib()"
                OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)

# install target
install(TARGETS chemkit-python DESTINATION ${PYTHON_SITE_PACKAGES})
