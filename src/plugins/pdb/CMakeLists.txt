find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

set(SOURCES
  pdbfileformat.cpp
  pdbplugin.cpp
)

set(MOC_HEADERS
  pdbplugin.h
)

qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

add_chemkit_plugin(pdb ${SOURCES} ${MOC_SOURCES})
target_link_libraries(pdb chemkit ${QT_LIBRARIES})
