find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES
  fhzfileformat.cpp
  fhzplugin.cpp
)

add_chemkit_plugin(fhz ${SOURCES})
target_link_libraries(fhz chemkit ${QT_LIBRARIES})
