find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

set(SOURCES
  mcdllineformat.cpp
  mcdlplugin.cpp
  mcdlreader.cpp
)

set(MOC_HEADERS
  mcdlplugin.h
)

qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

add_chemkit_plugin(mcdl ${SOURCES} ${MOC_SOURCES})
target_link_libraries(mcdl chemkit ${QT_LIBRARIES})
