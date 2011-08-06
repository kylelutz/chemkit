# setup config.h
# set a variable for the operating system
if(UNIX)
  set(CHEMKIT_OS_UNIX TRUE)
elseif(WIN32)
  set(CHEMKIT_OS_WIN32 TRUE)
elseif(APPLE)
  set(CHEMKIT_OS_MAC TRUE)
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)

find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

find_package(Boost COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(HEADERS
  alphashape.h
  aminoacid.h
  atom.h
  atom-inline.h
  atomtyper.h
  bond.h
  bond-inline.h
  bondpredictor.h
  chemkit.h
  config.h
  conformer.h
  constants.h
  coordinates.h
  delaunaytriangulation.h
  dynamiclibrary.h
  element.h
  element-inline.h
  forcefield.h
  forcefieldatom.h
  forcefieldcalculation.h
  forcefieldcalculation-inline.h
  forcefieldinteractions.h
  foreach.h
  fragment.h
  fragment-inline.h
  genericfile.h
  genericfile-inline.h
  geometry.h
  geometry-inline.h
  internalcoordinates.h
  lineformat.h
  matrix.h
  moiety.h
  moleculardescriptor.h
  moleculargraph.h
  molecularsurface.h
  molecule.h
  molecule-inline.h
  moleculealigner.h
  moleculefile.h
  moleculefileformat.h
  moleculefileformatadaptor.h
  moleculefileformatadaptor-inline.h
  moleculeobserver.h
  nucleotide.h
  partialchargepredictor.h
  plugin.h
  plugin-inline.h
  pluginmanager.h
  pluginmanager-inline.h
  point3.h
  polymer.h
  polymerchain.h
  polymerfile.h
  polymerfileformat.h
  quaternion.h
  residue.h
  ring.h
  ring-inline.h
  scalarfield.h
  trajectory.h
  trajectoryfile.h
  trajectoryfileformat.h
  trajectoryframe.h
  unitcell.h
  variant.h
  variant-inline.h
  vector3.h
)

set(SOURCES
  alphashape.cpp
  aminoacid.cpp
  atom.cpp
  atomtyper.cpp
  bond.cpp
  bondpredictor.cpp
  chemkit.cpp
  conformer.cpp
  coordinates.cpp
  delaunaytriangulation.cpp
  dynamiclibrary.cpp
  element.cpp
  forcefield.cpp
  forcefieldatom.cpp
  forcefieldcalculation.cpp
  forcefieldinteractions.cpp
  fragment.cpp
  geometry.cpp
  internalcoordinates.cpp
  lineformat.cpp
  matrix.cpp
  moiety.cpp
  moleculardescriptor.cpp
  moleculargraph.cpp
  moleculargraph_rppath.cpp
  moleculargraph_vf2.cpp
  molecularsurface.cpp
  molecule.cpp
  moleculealigner.cpp
  moleculefile.cpp
  moleculefileformat.cpp
  moleculeobserver.cpp
  nucleotide.cpp
  partialchargepredictor.cpp
  plugin.cpp
  pluginmanager.cpp
  polymer.cpp
  polymerchain.cpp
  polymerfile.cpp
  polymerfileformat.cpp
  residue.cpp
  ring.cpp
  scalarfield.cpp
  trajectory.cpp
  trajectoryfile.cpp
  trajectoryfileformat.cpp
  trajectoryframe.cpp
  unitcell.cpp
)

set(MOC_HEADERS
  plugin.h
)

qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

add_definitions(
  -DCHEMKIT_LIBRARY
  -DQT_NO_KEYWORDS
)

add_chemkit_library(chemkit ${SOURCES} ${MOC_SOURCES})
target_link_libraries(chemkit ${Boost_LIBRARIES} ${QT_LIBRARIES})

# install header files
install(FILES ${HEADERS} DESTINATION include/chemkit/)
