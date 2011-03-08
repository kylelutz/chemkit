TEMPLATE = lib
TARGET = chemkit
DESTDIR = ../../lib
INCLUDEPATH += ../../include
LIBS += -lblas \
    -llapack
QT -= gui
DEFINES += CHEMKIT_LIBRARY
SOURCES += alphashape.cpp \
    aminoacid.cpp \
    atom.cpp \
    atommapping.cpp \
    atomtyper.cpp \
    biochemicalfile.cpp \
    biochemicalfileformat.cpp \
    bond.cpp \
    bondpredictor.cpp \
    chemicalfile.cpp \
    chemicalfileformat.cpp \
    chemkit.cpp \
    conformer.cpp \
    coordinatematrix.cpp \
    delaunaytriangulation.cpp \
    element.cpp \
    forcefield.cpp \
    forcefieldatom.cpp \
    forcefieldcalculation.cpp \
    forcefieldinteractions.cpp \
    fragment.cpp \
    geometry.cpp \
    lineformat.cpp \
    matrix.cpp \
    moiety.cpp \
    moleculardescriptor.cpp \
    moleculargraph.cpp \
    moleculargraph_rppath.cpp \
    moleculargraph_vf2.cpp \
    molecularsurface.cpp \
    molecule.cpp \
    moleculealigner.cpp \
    moleculewatcher.cpp \
    nucleicacid.cpp \
    nucleicacidchain.cpp \
    nucleotide.cpp \
    partialchargepredictor.cpp \
    plugin.cpp \
    pluginmanager.cpp \
    protein.cpp \
    proteinchain.cpp \
    residue.cpp \
    ring.cpp \
    scalarfield.cpp
HEADERS += alphashape.h \
    aminoacid.h \
    atom.h \
    atom-inline.h \
    atommapping.h \
    atomtyper.h \
    biochemicalfile.h \
    biochemicalfileformat.h \
    blas.h \
    bond.h \
    bond-inline.h \
    bondpredictor.h \
    chemicalfile.h \
    chemicalfileformat.h \
    chemkit.h \
    commainitializer.h \
    commainitializer-inline.h \
    conformer.h \
    constants.h \
    coordinatematrix.h \
    delaunaytriangulation.h \
    element.h \
    element-inline.h \
    forcefield.h \
    forcefieldatom.h \
    forcefieldcalculation.h \
    forcefieldcalculation-inline.h \
    forcefieldinteractions.h \
    fragment.h \
    fragment-inline.h \
    genericmatrix.h \
    genericmatrix-inline.h \
    genericpoint.h \
    genericpoint-inline.h \
    genericquaternion.h \
    genericquaternion-inline.h \
    genericvector.h \
    genericvector-inline.h \
    geometry.h \
    lapack.h \
    lineformat.h \
    matrix.h \
    moiety.h \
    moleculardescriptor.h \
    moleculargraph.h \
    molecularsurface.h \
    molecule.h \
    molecule-inline.h \
    moleculealigner.h \
    moleculewatcher.h \
    nucleicacid.h \
    nucleicacidchain.h \
    nucleotide.h \
    partialchargepredictor.h \
    plugin.h \
    pluginmanager.h \
    point.h \
    point-inline.h \
    protein.h \
    proteinchain.h \
    quaternion.h \
    quaternion-inline.h \
    residue.h \
    ring.h \
    ring-inline.h \
    scalarfield.h \
    staticmatrix.h \
    staticmatrix-inline.h \
    staticvector.h \
    staticvector-inline.h \
    vector.h \
    vector-inline.h
target.path = /usr/lib
INSTALLS += target
headers.files = $$HEADERS
headers.path = /usr/include/chemkit/
INSTALLS += headers
