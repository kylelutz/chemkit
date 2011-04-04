find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

set(SOURCES
  fhzfileformat.cpp
  fhzplugin.cpp
)

set(MOC_HEADERS
  fhzplugin.h
)

qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

add_chemkit_plugin(fhz ${SOURCES} ${MOC_SOURCES})
target_link_libraries(fhz ${QT_LIBRARIES})
