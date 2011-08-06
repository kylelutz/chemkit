find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES
  mcdllineformat.cpp
  mcdlplugin.cpp
  mcdlreader.cpp
)

add_chemkit_plugin(mcdl ${SOURCES})
target_link_libraries(mcdl chemkit ${QT_LIBRARIES})
