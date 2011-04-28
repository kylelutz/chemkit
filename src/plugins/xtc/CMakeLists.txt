find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

set(SOURCES
  xtcfileformat.cpp
  xtcplugin.cpp
)

set(MOC_HEADERS
  xtcplugin.h
)

aux_source_directory(../../3rdparty/xdrf/ XDRF_SOURCES)

qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

add_chemkit_plugin(xtc ${SOURCES} ${MOC_SOURCES} ${XDRF_SOURCES})
target_link_libraries(xtc ${QT_LIBRARIES})
