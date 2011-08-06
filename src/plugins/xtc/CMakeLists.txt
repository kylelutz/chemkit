find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

set(SOURCES
  xtcfileformat.cpp
  xtcplugin.cpp
)

aux_source_directory(../../3rdparty/xdrf/ XDRF_SOURCES)

add_chemkit_plugin(xtc ${SOURCES} ${XDRF_SOURCES})
target_link_libraries(xtc chemkit ${QT_LIBRARIES})
