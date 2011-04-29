find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

include_directories(../../3rdparty/lemon)

set(SOURCES
  kekulizer.cpp
  smilesgraph.cpp
  smileslineformat.cpp
  smilesplugin.cpp
)

set(MOC_HEADERS
  smilesplugin.h
)

qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

add_chemkit_plugin(smiles ${SOURCES} ${MOC_SOURCES})
target_link_libraries(smiles chemkit ${QT_LIBRARIES})
