find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

set(SOURCES
  mcdllineformat.cpp
  mcdlplugin.cpp
  mcdlreader.cpp
)

add_chemkit_plugin(mcdl ${SOURCES})
target_link_libraries(mcdl ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})
