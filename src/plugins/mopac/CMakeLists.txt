find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

set(SOURCES
  mopacplugin.cpp
  mopcrtfileformat.cpp
  mopinfileformat.cpp
)

set(MOC_HEADERS
  mopacplugin.h
)

qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

add_chemkit_plugin(mopac ${SOURCES} ${MOC_SOURCES})
target_link_libraries(mopac chemkit ${QT_LIBRARIES})
