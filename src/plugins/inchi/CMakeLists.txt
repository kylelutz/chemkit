find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

set(SOURCES
  inchikeylineformat.cpp
  inchilineformat.cpp
  inchiplugin.cpp
)

set(MOC_HEADERS
  inchiplugin.h
)

aux_source_directory(../../3rdparty/inchi/ IUPAC_INCHI_SOURCES)

qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

add_chemkit_plugin(inchi ${SOURCES} ${MOC_SOURCES} ${IUPAC_INCHI_SOURCES})
target_link_libraries(inchi chemkit ${QT_LIBRARIES})
