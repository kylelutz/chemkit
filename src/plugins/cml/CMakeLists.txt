find_package(Qt4 4.6 COMPONENTS QtCore QtXml REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
set(QT_USE_QTXML TRUE)
include(${QT_USE_FILE})

set(SOURCES
  cmlfileformat.cpp
  cmlplugin.cpp
)

set(MOC_HEADERS
  cmlplugin.h
)

qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

add_chemkit_plugin(cml ${SOURCES} ${MOC_SOURCES})
target_link_libraries(cml ${QT_LIBRARIES})
