if(NOT ${CHEMKIT_WITH_MD_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io md md-io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

set(SOURCES
  xtcfileformat.cpp
  xtcplugin.cpp
)

aux_source_directory(../../3rdparty/xdrf/ XDRF_SOURCES)

add_chemkit_plugin(xtc ${SOURCES} ${XDRF_SOURCES})
target_link_libraries(xtc ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})
