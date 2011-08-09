find_package(Qt4 4.6 COMPONENTS QtCore QtXml REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
set(QT_USE_QTXML TRUE)
include(${QT_USE_FILE})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES
  cmlfileformat.cpp
  cmlplugin.cpp
)

add_chemkit_plugin(cml ${SOURCES})
target_link_libraries(cml chemkit chemkit-io ${QT_LIBRARIES})
