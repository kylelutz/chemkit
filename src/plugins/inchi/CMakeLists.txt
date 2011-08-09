find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES
  inchikeylineformat.cpp
  inchilineformat.cpp
  inchiplugin.cpp
)

aux_source_directory(../../3rdparty/inchi/ IUPAC_INCHI_SOURCES)

add_chemkit_plugin(inchi ${SOURCES} ${IUPAC_INCHI_SOURCES})
target_link_libraries(inchi chemkit chemkit-io ${QT_LIBRARIES})
