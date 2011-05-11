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
  atommapping.h
  atomtyper.h
  bond.h
  bond-inline.h
  bondpredictor.h
  chemkit.h
  commainitializer.h
  commainitializer-inline.h
  conformer.h
  constants.h
  coordinates.h
  delaunaytriangulation.h
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
  genericpoint.h
  genericpoint-inline.h
  genericquaternion.h
  genericquaternion-inline.h
  genericvector.h
  genericvector-inline.h
  geometry.h
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
  moleculewatcher.h
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
  staticmatrix.h
  staticmatrix-inline.h
  staticvector.h
  staticvector-inline.h
  trajectory.h
  trajectoryfile.h
  trajectoryfileformat.h
  trajectoryframe.h
  unitcell.h
  variant.h
  vector3.h
)

set(SOURCES
  alphashape.cpp
  aminoacid.cpp
  atom.cpp
  atommapping.cpp
  atomtyper.cpp
  bond.cpp
  bondpredictor.cpp
  chemkit.cpp
  conformer.cpp
  coordinates.cpp
  delaunaytriangulation.cpp
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
  moleculewatcher.cpp
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
  moleculewatcher.h
  plugin.h
  pluginmanager.h
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
